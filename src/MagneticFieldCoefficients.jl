import Base.getindex, Base.setindex!,
    Base.isapprox,
    Base.==,
    Base.+, Base.-, Base.*, Base./,
    Base.write,
    Base.size, Base.length,
    Base.iterate, Base.hash,
    Base.eltype

# Spherical harmonic coefficients describing a magnetic field
mutable struct MagneticFieldCoefficients
    coeffs::Array{SphericalHarmonicCoefficients,2} # coefficients
    radius::Float64 # radius of measured sphere
    center::Array{Float64,2} # center of measured sphere (depends on expansion point)
    ffp::Union{Array{Float64,2},Nothing} # field-free-point (if available, depends on expansion point)

    function MagneticFieldCoefficients(
            coeffs::Array{SphericalHarmonicCoefficients,2},
            radius::Float64,
            center::Array{Float64,2},
            ffp::Union{Array{Float64,2},Nothing})

        # test sizes of the arrays
        if size(coeffs, 1) != 3
            throw(DimensionMismatch("The coefficient matrix needs 3 entries (x,y,z) in the first dimension, not $(size(coeffs,1))"))
        elseif size(coeffs, 2) != size(center, 2)
            throw(DimensionMismatch("The number of patches of the coefficients and center does not match: $(size(coeffs,2)) != $(size(center,2))"))
        elseif !isnothing(ffp)
            if size(ffp, 1) != 3
                throw(DimensionMismatch("The FFP matrix needs 3 entries (x,y,z) in the first dimension, not $(size(coeffs,1))"))
            elseif size(coeffs, 2) != size(ffp, 2)
                throw(DimensionMismatch("The number of patches of the coefficients and FFPs does not match: $(size(coeffs,2)) != $(size(ffp,2))"))
            end
        end

        return new(coeffs, radius, center, ffp)
    end

end

# some other constructors
MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}) =
    MagneticFieldCoefficients(coeffs, 0.0)
MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}, radius::Float64) =
    MagneticFieldCoefficients(coeffs, radius, [0.0, 0.0, 0.0])
MagneticFieldCoefficients(
    coeffs::Array{SphericalHarmonicCoefficients,2},
    radius::Float64,
    center::VecOrMat{Float64}) = 
        MagneticFieldCoefficients(coeffs, radius, center, nothing)
MagneticFieldCoefficients(
    coeffs::Array{SphericalHarmonicCoefficients,2},
    radius::Float64,
    center::Vector{Float64},
    ffp::Union{Matrix{Float64},Nothing}) = 
        MagneticFieldCoefficients(
            coeffs,
            radius,
            hcat([center for p = 1:size(coeffs, 2)]...),
            ffp)

function MagneticFieldCoefficients(
            L::Int,
            R::Float64 = 1.0,
            solid::Bool = true;
            # optional
            radius::Float64 = 0.0,
            center::VecOrMat{Float64} = [0.0, 0.0, 0.0],
            numPatches::Int = 1)

    if L < 0
        throw(DomainError(L, "Input vector needs to be of size (L+1)², where L ∈ ℕ₀."))
    end

    return MagneticFieldCoefficients(
        [SphericalHarmonicCoefficients(L, R, solid) for j = 1:3, p = 1:numPatches],
        radius,
        center)
end

# constructor using t-design
MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}, tDesign::SphericalTDesign, ffp = nothing) = 
    MagneticFieldCoefficients(coeffs, ustrip(Unitful.m(tDesign.radius)), ustrip.(Unitful.m.(tDesign.center)), ffp)


# read coefficients/measurement data from an HDF5-file
function MagneticFieldCoefficients(path::String)

    coeffsMF = h5open(path, "r") do file
        if haskey(HDF5.root(file), "/coeffs")
            # coefficients available    
            return loadCoefficients(path)

        elseif haskey(HDF5.root(file), "/fields") &&
               haskey(HDF5.root(file), "/positions/tDesign")
            # calculate coefficients
            return loadTDesignCoefficients(path)

        else
            throw(ErrorException("Unknown file structure."))
        end
    end

    return coeffsMF
end

# read coefficients from an HDF5-file
function loadCoefficients(path::String)

    # load spherical harmonic coefficients
    shcoeffs = SphericalHarmonicCoefficients(path)

    coeffsMF = h5open(path, "r") do file
        if !haskey(HDF5.root(file), "/radius") || !haskey(HDF5.root(file), "/center")
            @warn "The file does not provide all necessary measurement informations."
        end

        # set all informations not given in the file to 0 or nothing
        radius = haskey(HDF5.root(file), "/radius") ? read(file, "/radius") : 0.0
        center = haskey(HDF5.root(file), "/center") ? read(file, "/center") : zeros(size(shcoeffs))
        ffp = haskey(HDF5.root(file), "/ffp") ? read(file, "/ffp") : nothing

        return MagneticFieldCoefficients(shcoeffs, radius, center, ffp)
    end

    return coeffsMF
end

# load coefficients from measurement data
# TODO: This should be merged with and moved to MPIFiles
function loadTDesignCoefficients(filename::String)
    field, radius, N, t, center, correction = h5open(filename, "r") do file
        field = read(file, "/fields")   # measured field (size: 3 x #points x #patches)
        radius = read(file, "/positions/tDesign/radius") # radius of the measured ball
        N = read(file, "/positions/tDesign/N")           # number of points of the t-design
        t = read(file, "/positions/tDesign/t")           # t of the t-design
        center = read(file, "/positions/tDesign/center") # center of the measured ball
        correction = read(file, "/sensor/correctionTranslation")
        return field, radius, N, t, center, correction
    end
    tDes = MPIFiles.loadTDesign(Int(t), N, radius * u"m", center .* u"m")
    coeffs = magneticField(tDes, field)
    for c = 1:size(coeffs, 2), j = 1:3
        coeffs[j, c] = SphericalHarmonicExpansions.translation(coeffs[j, c], correction[:, j])
    end

    coeffs_MF = MagneticFieldCoefficients(coeffs, radius, center)

    return coeffs_MF
end


# write coefficients to an HDF5-file
function write(path::String, coeffs::MagneticFieldCoefficients)

    # save SphericalHarmonicCoefficients
    write(path, coeffs.coeffs)

    # add field informations
    radius = coeffs.radius
    center = coeffs.center
    ffp = coeffs.ffp

    h5open(path, "cw") do file
        write(file, "/radius", radius)
        write(file, "/center", center)
        if ffp !== nothing
            write(file, "/ffp", ffp)
        end
    end
end

# Sizes
size(mfc::MagneticFieldCoefficients, kargs...) = size(mfc.coeffs, kargs...)
length(mfc::MagneticFieldCoefficients) = size(mfc,2)

# indexing
getindex(mfc::MagneticFieldCoefficients, i) = 
    MagneticFieldCoefficients(
        reshape(getindex(mfc.coeffs, :, i), 3, length(i)),
        mfc.radius, reshape(getindex(mfc.center, :, i), 3, length(i)),
        isnothing(mfc.ffp) ? nothing : reshape(getindex(mfc.ffp, :, i), 3, length(i)))

function setindex!(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients, i)
    # setindex not possible in some cases 
    if mfc1.radius != mfc2.radius
        throw(DomainError([mfc1.radius, mfc2.radius], "Coefficients do not have the same measurement radius."))
    elseif !isnothing(mfc1.ffp) && isnothing(mfc2.ffp)
        throw(DomainError([mfc1.ffp, mfc2.ffp], "Coefficients do not provide an FFP."))
    end

    # set coefficients and center
    setindex!(mfc1.coeffs, mfc2.coeffs, :, i)
    setindex!(mfc1.center, mfc2.center, :, i)
    # set FFP (ignore mfc2.ffp if mfc1.ffp === nothing)
    if !isnothing(mfc1.ffp) && !isnothing(mfc2.ffp)
        setindex!(mfc1.ffp, mfc2.ffp, :, i)
    end
    return nothing
end

# Iterator
iterate(mfc::MagneticFieldCoefficients, state = 1) = state <= length(mfc) ? (mfc[state], state + 1) : nothing

eltype(mfc::MagneticFieldCoefficients) = MagneticFieldCoefficients # element type

# Operations on MagneticFieldCoefficients
function isapprox(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients;
                  kargs...)

    val = all(isapprox.(mfc1.coeffs, mfc2.coeffs; kargs...)) &&
        isapprox(mfc1.radius, mfc2.radius; kargs...) &&
        isapprox(mfc1.center, mfc2.center; kargs...)

    if isnothing(mfc1.ffp) && isnothing(mfc2.ffp)
        return val
    elseif !isnothing(mfc1.ffp) && !isnothing(mfc2.ffp)
        return val && isapprox(mfc1.ffp, mfc2.ffp; kargs...)
    else
        @info "Only one of the coefficients has FFPs. Applying isapprox to the other values yields $val."
        return false
    end
end
==(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients) =
    mfc1.coeffs == mfc2.coeffs &&
    mfc1.radius == mfc2.radius &&
    mfc1.center == mfc2.center &&
    mfc1.ffp == mfc2.ffp

function hash(mfc::MagneticFieldCoefficients, h::UInt64)
    h = hash(:MagneticFieldCoefficients, h)
    h = hash(mfc.coeffs, h)
    h = hash(mfc.center, h)
    h = hash(mfc.radius, h)
    h = hash(mfc.ffp, h)
    return h
end

"""
    +(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; force::Bool=false)

`force = true` adds the coefficients even if the radius or center are not equal (set to values of the first coefficients).
"""
function +(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients;
           force::Bool = false)

    if force
        return MagneticFieldCoefficients(mfc1.coeffs .+ mfc2.coeffs,
                                         mfc1.radius, mfc1.center)
    end
    if mfc1.radius != mfc2.radius
        throw(DomainError([mfc1.radius, mfc2.radius], "Coefficients do not have the same measurement radius. (Use `force = true`.)"))
    end
    if mfc1.center != mfc2.center
        throw(DomainError([mfc1.center, mfc2.center], "Coefficients do not have the same measurement center. (Use `force = true`.)"))
    end
    return MagneticFieldCoefficients(mfc1.coeffs .+ mfc2.coeffs, mfc1.radius, mfc1.center)
end

"""
    -(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; force::Bool=false)

`force = true` subtracts the coefficients even if the radius or center are not equal (set to values of the first coefficients).
"""
function -(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients;
           force::Bool = false)
    if force
        return MagneticFieldCoefficients(
            mfc1.coeffs .- mfc2.coeffs,
            mfc1.radius,
            mfc1.center,
        )
    end
    if mfc1.radius != mfc2.radius
        throw(DomainError([mfc1.radius, mfc2.radius], "Coefficients do not have the same measurement radius."))
    end
    if mfc1.center != mfc2.center
        throw(DomainError([mfc1.center, mfc2.center], "Coefficients do not have the same measurement center."))
    end
    return MagneticFieldCoefficients(mfc1.coeffs .- mfc2.coeffs, mfc1.radius, mfc1.center)
end

+(mfc::MagneticFieldCoefficients, value::Number) =
    MagneticFieldCoefficients(mfc.coeffs .+ value, mfc.radius, mfc.center, mfc.ffp)
-(mfc::MagneticFieldCoefficients, value::Number) =
    MagneticFieldCoefficients(mfc.coeffs .- value, mfc.radius, mfc.center, mfc.ffp)
*(mfc::MagneticFieldCoefficients, value::Number) =
    MagneticFieldCoefficients(mfc.coeffs .* value, mfc.radius, mfc.center, mfc.ffp)
/(mfc::MagneticFieldCoefficients, value::Number) =
    MagneticFieldCoefficients(mfc.coeffs ./ value, mfc.radius, mfc.center, mfc.ffp)
+(value::Number, mfc::MagneticFieldCoefficients) = +(mfc::MagneticFieldCoefficients, value)
-(value::Number, mfc::MagneticFieldCoefficients) =
    MagneticFieldCoefficients(value .- mfc.coeffs, mfc.radius, mfc.center, mfc.ffp)
*(value::Number, mfc::MagneticFieldCoefficients) = *(mfc::MagneticFieldCoefficients, value)

## get some field-specific values
"""
    getOffset(mfc::MagneticFieldCoefficients, idx::AbstractUnitRange{Int64}=axes(mfc.coeffs,2))

Get the offset of the field described by mfc[idx].
"""
getOffset(mfc::MagneticFieldCoefficients, idx::AbstractUnitRange{Int64} = axes(mfc.coeffs, 2)) = [c[0, 0] for c in mfc.coeffs[:,idx]]
"""
    getGradient(mfc::MagneticFieldCoefficients, idx::AbstractUnitRange{Int64}=axes(mfc.coeffs,2))

Get the gradient of the field described by mfc[idx].
"""
function getGradient(
            mfc::MagneticFieldCoefficients,
            idx::AbstractUnitRange{Int64} = axes(mfc.coeffs, 2))
    # ideal gradient would be [[mfc.coeffs[1,p][1,1], mfc.coeffs[2,p][1,-1], mfc.coeffs[3,p][1,0]] for p in idx]

    G = Vector{Float64}[]
    for p in idx
        # calculate jacobian matrix
        @polyvar x y z
        expansion = sphericalHarmonicsExpansion.(mfc.coeffs[:, p], [x], [y], [z])
        jexp = differentiate(expansion, [x, y, z])
        push!(G, [jexp[i, i]((x, y, z) => [0.0, 0.0, 0.0]) for i = 1:3]) # gradient = diag(jacobian matrix)
    end

    return G
end
"""
    getJacobian(mfc::MagneticFieldCoefficients, idx::AbstractUnitRange{Int64}=axes(mfc.coeffs,2))

Get the Jacobian matrix of the field described by mfc[idx].
"""
function getJacobian(
            mfc::MagneticFieldCoefficients,
            idx::AbstractUnitRange{Int64} = axes(mfc.coeffs, 2))
    # ideal Jacobian matrix would be [vcat([[c[1,1], c[1,-1], c[1,0]]' for c in mfc.coeffs[:,p]]...) for p in idx]

    J = Matrix{Float64}[]
    for p in idx
        # calculate jacobian matrix
        @polyvar x y z
        expansion = sphericalHarmonicsExpansion.(mfc.coeffs[:, p], [x], [y], [z])
        jexp = differentiate(expansion, [x, y, z])
        push!(J, [jexp[i, j]((x, y, z) => [0.0, 0.0, 0.0]) for i = 1:3, j = 1:3]) # jacobian matrix
    end

    return J
end

## Load coefficients from t-design measurement ##
"""
    magneticField(tDesign::SphericalTDesign, field::Union{AbstractArray{T,2},AbstractArray{T,3}};
		       L::Int=Int(tDesign.T/2),
		       calcSolid::Bool=true) where T <: Real
*Description:*  Calculation of the spherical harmonic coefficients based on the measured t-design\\
 \\
*Input:*
- `tDesign`	- Measured t-design (type: SphericalTDesign)
- `field`       - Measured field (size = (J,N,C)) with J <= 3
**kwargs:**
- `L`           - Order up to which the coeffs be calculated (default: t/2)
- `calcSolid`   - Boolean (default: true)\\
    false -> spherical coefficients\\
    true -> solid coefficients
*Output:*
- `coeffs`    - spherical/solid coefficients, type: Array{SphericalHarmonicCoefficients}(3,C)
"""
function magneticField(
                tDesign::SphericalTDesign,
                field::Union{AbstractArray{T,2},AbstractArray{T,3}},
                L::Int = floor(Int, tDesign.T / 2),
                calcSolid::Bool = true) where {T<:Real}

    # get tDesign positions [m] and removing the unit
    # coordinates
    coords = Float64.(ustrip.(Unitful.m.(hcat([p for p in tDesign]...))))

    # radius
    R = Float64(ustrip(Unitful.m(tDesign.radius)))

    # center
    center = Float64.(ustrip.(Unitful.m.(tDesign.center)))

    return magneticField(coords, field, R, center, L, calcSolid)
end

function magneticField(
                coords::AbstractArray{T,2},
                field::Union{AbstractArray{T,2},AbstractArray{T,3}},
                R::T,
                center::Vector{T},
                L::Int,
                calcSolid::Bool = true) where {T<:Real}

    # transpose coords if its dimensions do not fit
    if size(coords, 1) != 3
        coords = coords'
    end

    # test dimensions of field array
    if size(field, 1) > 3
        throw(DimensionMismatch("The measured field has more than 3 entries in the first dimension: $(size(field,1))"))
    elseif size(field, 2) != size(coords, 2)
        throw(DimensionMismatch("The field vector does not match the size of the tdesign: $(size(field,2)) != $(size(coords,2))"))
    end

    coeffs = Array{SphericalHarmonicCoefficients}(undef, size(field, 1), size(field, 3))

    # rescale coordinates to t-design on unit sphere
    coords = coords .- center
    coords *= 1 / R
    for c in axes(field, 3)
        # calculation of the coefficients
        for j in axes(field, 1)
            coeffs[j, c] = SphericalHarmonicExpansions.sphericalQuadrature(field[j, :, c], coords', L)
            coeffs[j, c].R = R

            normalize!(coeffs[j, c], R)

            # convert spherical into solid coefficients
            if calcSolid
                solid!(coeffs[j, c])
            end
        end
    end

    return coeffs
end

# Shift coefficients into new expansion point
function shift!(coeffsMF::MagneticFieldCoefficients, v::Matrix{T}) where {T<:Real}

    # test dimensions of shifting vector/array
    if size(v, 1) != 3
        throw(DimensionMismatch("The shifting vector/matrix needs 3 entries but it has $(size(field,1))"))
    elseif size(v, 2) != size(coeffsMF, 2)
        throw(DimensionMismatch("The number of patches do not coincide: $(size(v,2)) != $(size(coeffsMF,2))"))
    end

    # shift coefficients by v
    for p = 1:size(coeffsMF, 2)
        coeffsMF.coeffs[:, p] = SphericalHarmonicExpansions.translation.(coeffsMF.coeffs[:, p], [v[:, p]])
    end

    # change center and FFP position due to new expansion point
    coeffsMF.center -= v
    if !isnothing(coeffsMF.ffp)
        coeffsMF.ffp -= v
    end

    return nothing
end

function shift!(coeffsMF::MagneticFieldCoefficients, v::Vector{T}) where {T<:Real}

    # create matrix from vector
    v = hcat([v for i = 1:size(coeffsMF, 2)]...)

    # shift coefficients
    shift!(coeffsMF, v)

    return nothing
end

function shift(coeffsMF::MagneticFieldCoefficients, v::Union{Vector{T},Matrix{T}}) where {T<:Real}

    coeffsMF_shifted = deepcopy(coeffsMF)
    shift!(coeffsMF_shifted, v)

    return coeffsMF_shifted
end

# shift coefficients into FFP
"""
    shiftFFP!(coeffsMF::MagneticFieldCoefficients)

Shift magnetic-field coefficients into FFP (and calculate it if not available).
"""
function shiftFFP!(coeffsMF::MagneticFieldCoefficients)

    # calculate FFP when not available
    if isnothing(coeffsMF.ffp)
        @info "Calculate FFP."
        findFFP!(coeffsMF)
    end

    # shift coefficients into FFP
    shift!(coeffsMF, coeffsMF.ffp)

    return nothing
end

## Calculate the FFP of the magnetic field ##

# Volumes where the FFP should be calculated
@enum Volume xyz xy xz yz

# get index set of the volume
function getVolumeIndex(vol::Symbol)  
    if vol == :xyz
        indVol = 1:3
    elseif vol == :xy
        indVol = 1:2
    elseif vol == :xz
        indVol = [1,3]
    elseif vol == :yz
        indVol = 2:3
    else # default/error
        indVol = 1:3
        @warn "Volume $vol not defined. Default 3D volume used."
    end

    return indVol
end
getVolumeIndex(vol::Union{String, Volume}) = getVolumeIndex(Symbol(vol))

"""
    findFFP(coeffsMF::MagneticFieldCoefficients; 
            vol::Union{Symbol, String, Volume} = xyz, 
            start::Union{Vector{Vector{Float64}},Vector{Float64}} = [[0.0;0.0;0.0]], 
            returnasmatrix::Bool=true)

*Description:*  Newton method to find the FFPs of the magnetic fields\\
 \\
*Input:*
- `coeffsMF`   - MagneticFieldCoefficients
**kwargs:**
- `vol` - provide the volume where the FFP search should be done\\
            options: $(instances(Volume))\\
            default: xyz
- `start` - start vector (important for 2D FFP search for the third coordinate)\\
            default: [0; 0; 0]
- `returnasmatrix` - Boolean\\
        true  -> return FFPs as Matrix with size (3,#Patches) (default)\\
        false -> return FFPs as Array of NLsolve.SolverResults with size #Patches
*Output:*
- `ffp` - FFPs of the magnetic field
"""
function findFFP(coeffsMF::MagneticFieldCoefficients;  
                 vol::Union{Symbol, String, Volume} = xyz, 
                 start::Union{Vector{Vector{Float64}},Vector{Float64}} = [[0.0;0.0;0.0]], # start value
                 returnasmatrix::Bool = true)

    # get spherical harmonic expansion
    @polyvar x y z
    expansion = sphericalHarmonicsExpansion.(coeffsMF.coeffs, x, y, z)

    # return all FFPs in a matrix or as array with the solver results
    ffp = returnasmatrix ? zeros(size(expansion)) : Array{NLsolve.SolverResults{Float64}}(undef, size(expansion, 2))

    # correct size of start vector
    if typeof(start) <: Vector{Float64}
        start = [start]
    end
    if length.(start) != 3 .* ones(length(start))
        throw(DimensionMismatch("Length of start vector needs to be 3, not $(length.(start))."))
    end
    if length(start) == 1 && size(coeffsMF,2) !== 1
        start = repeat(start,size(coeffsMF,2))
    elseif length(start) != size(coeffsMF,2)
        throw(DimensionMismatch("The numbers of start vectors ($(length(start))) has to be equal to the number of patches ($(size(coeffsMF,2)))."))
    end

    # get index set of the volume
    indVol = getVolumeIndex(vol)  

    # find FFPs
    for c in axes(expansion, 2)

        pol = [expansion[j,c] for j in indVol]
        dpol = vcat([hcat(differentiate.(p,(x,y,z))[indVol]...) for p in pol]...)

        function f!(fvec, xx)
            # for the planes use one value of start for the third dimension
            val = start[c]
            val[indVol] = xx
            # calculate field
            fvec[:] = [p((x,y,z) => val) for p in pol]
        end

        function g!(fjac, xx)
            # for the planes use one value of start for the third dimension
            val = start[c]
            val[indVol] = xx
            # calculate field derivative
            fjac[:,:] = [dp((x,y,z) => val) for dp in dpol]
        end

        if returnasmatrix
            ffp[:,c] = start[c] # start value needed for dimensions not used for FFP calculation
            ffp[indVol, c] = nlsolve(f!, g!, start[c][indVol], method = :newton, ftol = 1e-16).zero
        else
            ffp[c] = nlsolve(f!, g!, start[c][indVol], method = :newton, ftol = 1e-16)
        end
    end

    return ffp
end

"""
    ffp = findFFP!(coeffsMF::MagneticFieldCoefficients)

Find FFP and set it as coeffsMF.ffp.
"""
function findFFP!(coeffsMF::MagneticFieldCoefficients; kargs...)
    # find FFP
    ffp = findFFP(coeffsMF; returnasmatrix = true, kargs...)
    # set FFP
    coeffsMF.ffp = ffp

    return ffp
end
