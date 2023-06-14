import Base.isapprox, Base.==,
        Base.+, Base.-, Base.*, Base./, 
        Base.write, Base.size

# Spherical harmonic coefficients describing a magnetic field
mutable struct MagneticFieldCoefficients
    coeffs::Array{SphericalHarmonicCoefficients,2} # coefficients
    radius::Float64 # radius of measured sphere
    center::Vector{Float64} # center of measured sphere
    ffp::Union{Array{Float64,2},Nothing} # field-free-point (if available)

    function MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}, radius::Float64,
                                        center::Vector{Float64}, ffp::Union{Array{Float64,2},Nothing})
        # test sizes of the arrays
        if size(coeffs,1) != 3
            throw(DimensionMismatch("The coefficient matrix needs 3 entries (x,y,z) in the first dimension, not $(size(coeffs,1))"))
        elseif ffp !== nothing
            if size(ffp,1) != 3
                throw(DimensionMismatch("The FFP matrix needs 3 entries (x,y,z) in the first dimension, not $(size(coeffs,1))"))
            elseif size(coeffs,2) != size(ffp,2)
                throw(DimensionMismatch("The number of patches of the coefficients and FFPs does not match: $(size(coeffs,2)) != $(size(ffp,2))"))
            end
        end

        return new(coeffs,radius,center,ffp)
    end

end

# some other constructors
MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}) = MagneticFieldCoefficients(coeffs, 0.0)
MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}, radius::Float64) = MagneticFieldCoefficients(coeffs,radius,[0.0,0.0,0.0])
MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}, radius::Float64, center::Vector{Float64}) = MagneticFieldCoefficients(coeffs,radius,center,nothing)
MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}, radius::Float64, ffp::Array{Float64,2}) = MagneticFieldCoefficients(coeffs,radius,[0.0,0.0,0.0],ffp)

function MagneticFieldCoefficients(L::Int)
    if L<0
        throw(DomainError(L,"Input vector needs to be of size (L+1)², where L ∈ ℕ₀."))
    end
    return MagneticFieldCoefficients(fill(SphericalHarmonicCoefficients(L), (3,1)))
end

# constructor using t-design
MagneticFieldCoefficients(coeffs::Array{SphericalHarmonicCoefficients,2}, tDesign::SphericalTDesign, ffp=nothing) = MagneticFieldCoefficients(coeffs,ustrip(Unitful.m(tDesign.radius)), ustrip.(Unitful.m.(tDesign.center)), ffp)


# read coefficients from an HDF5-file
function MagneticFieldCoefficients(path::String)

    # load spherical harmonic coefficients
    shcoeffs = SphericalHarmonicCoefficients(path)

    coeffsMF = h5open(path,"r") do file
        if haskey(HDF5.root(file), "/radius") 
            # file contains all relevant information
            radius = read(file, "/radius")
            center = read(file, "/center")
            if haskey(HDF5.root(file), "/ffp")
                ffp = read(file, "/ffp")
                return MagneticFieldCoefficients(shcoeffs, radius, center, ffp)
            else
                # field has not FFP -> ffp = nothing
                return MagneticFieldCoefficients(shcoeffs, radius, center)
            end
        else
            # convert file of SphericalHarmonicCoefficients into MagneticFieldCoefficients
            # -> set all additional informations to 0 or nothing
            # use radius = 0.01 as default value
            return MagneticFieldCoefficients(shcoeffs, 0.01)
        end
    end

    return coeffsMF
end

# write coefficients to an HDF5-file
function write(path::String, coeffs::MagneticFieldCoefficients)

    # save SphericalHarmonicCoefficients
    write(path,coeffs.coeffs)

    # add field informations
    radius = coeffs.radius
    center = coeffs.center
    ffp = coeffs.ffp

    h5open(path,"cw") do file
        write(file, "/radius", radius)
        write(file, "/center", center)
        if ffp !== nothing
            write(file, "/ffp", ffp)
        end
    end
end

# Size
size(mfc::MagneticFieldCoefficients) = size(mfc.coeffs)

# Operations on MagneticFieldCoefficients
function isapprox(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; kargs...)
    val = all(isapprox.(mfc1.coeffs,mfc2.coeffs;kargs...)) && isapprox(mfc1.radius,mfc2.radius;kargs...) && isapprox(mfc1.center,mfc2.center;kargs...)
    if mfc1.ffp === nothing && mfc2.ffp === nothing
        return  val
    elseif mfc1.ffp !== nothing || mfc2.ffp !== nothing
        @info "Only one of the coefficients has FFPs. Applying isapprox to the other values yields $val."
        return false
    else
        return val && isapprox(mfc1.ffp,mfc2.ffp;kargs...)
    end
end
==(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients) = 
    mfc1.coeffs == mfc2.coeffs && mfc1.radius == mfc2.radius && mfc1.center == mfc2.center && mfc1.ffp == mfc2.ffp

"""
    +(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; force::Bool=false)

`force = true` adds the coefficients even if the radius or center are not equal (set to values of the first coefficients).
"""
function +(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; force::Bool=false)
    if force 
        return MagneticFieldCoefficients(mfc1.coeffs .+ mfc2.coeffs, mfc1.radius, mfc1.center)
    end
    if mfc1.radius != mfc2.radius
        throw(DomainError([mfc1.radius,mfc2.radius],"Coefficients do not have the same measurement radius."))
    end
    if mfc1.center != mfc2.center
        throw(DomainError([mfc1.center,mfc2.center],"Coefficients do not have the same measurement center."))
    end
    return MagneticFieldCoefficients(mfc1.coeffs .+ mfc2.coeffs, mfc1.radius, mfc1.center)
end

"""
    -(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; force::Bool=false)

`force = true` subtracts the coefficients even if the radius or center are not equal (set to values of the first coefficients).
"""
function -(mfc1::MagneticFieldCoefficients, mfc2::MagneticFieldCoefficients; force::Bool=false)
    if force 
        return MagneticFieldCoefficients(mfc1.coeffs .- mfc2.coeffs, mfc1.radius, mfc1.center)
    end
    if mfc1.radius != mfc2.radius
        throw(DomainError([mfc1.radius,mfc2.radius],"Coefficients do not have the same measurement radius."))
    end
    if mfc1.center != mfc2.center
        throw(DomainError([mfc1.center,mfc2.center],"Coefficients do not have the same measurement center."))
    end
    return MagneticFieldCoefficients(mfc1.coeffs .- mfc2.coeffs, mfc1.radius, mfc1.center)
end

+(mfc::MagneticFieldCoefficients, value::Number) = MagneticFieldCoefficients(mfc.coeffs .+ value, mfc.radius, mfc.center, mfc.ffp)
-(mfc::MagneticFieldCoefficients, value::Number) = MagneticFieldCoefficients(mfc.coeffs .- value, mfc.radius, mfc.center, mfc.ffp)
*(mfc::MagneticFieldCoefficients, value::Number) = MagneticFieldCoefficients(mfc.coeffs .* value, mfc.radius, mfc.center, mfc.ffp)
/(mfc::MagneticFieldCoefficients, value::Number) = MagneticFieldCoefficients(mfc.coeffs ./ value, mfc.radius, mfc.center, mfc.ffp)
+(value::Number, mfc::MagneticFieldCoefficients) = +(mfc::MagneticFieldCoefficients, value)
-(value::Number, mfc::MagneticFieldCoefficients) = MagneticFieldCoefficients(value .- mfc.coeffs, mfc.radius, mfc.center, mfc.ffp)
*(value::Number, mfc::MagneticFieldCoefficients) = *(mfc::MagneticFieldCoefficients, value)


## Load coefficients from t-design measurement ##
"""
    magneticField(tDesign::SphericalTDesign, field::Union{AbstractArray{T,2},AbstractArray{T,3}}, 
		       x::Variable, y::Variable, z::Variable;
		       L::Int=Int(tDesign.T/2),
		       calcSolid::Bool=true) where T <: Real
*Description:*  Calculation of the spherical harmonic coefficients and expansion based on the measured t-design\\
 \\
*Input:*
- `tDesign`	- Measured t-design (type: SphericalTDesign)
- `field`       - Measured field (size = (J,N,C)) with J <= 3
- `x, y, z`     - Cartesian coordinates
**kwargs:**
- `L`           - Order up to which the coeffs be calculated (default: t/2)
- `calcSolid`   - Boolean (default: true)\\
    false -> spherical coefficients\\
    true -> solid coefficients
*Output:*
- `coeffs`    - spherical/solid coefficients, type: Array{SphericalHarmonicCoefficients}(3,C)
- `expansion` - related expansion (Cartesian polynomial), type: Array{AbstractPolynomialLike}(3,C)
- `func`      - expansion converted to a function, type: Array{Function}(3,C)
"""
function magneticField(tDesign::SphericalTDesign, field::Union{AbstractArray{T,2},AbstractArray{T,3}}, 
		       x::Variable, y::Variable, z::Variable;
		       L::Int=floor(Int,tDesign.T/2),
		       calcSolid::Bool=true) where T <: Real

  # get tDesign positions [m] and removing the unit
  # coordinates
  coords = Float64.(ustrip.(Unitful.m.(hcat([p for p in tDesign]...))))

  # radius
  R = Float64(ustrip(Unitful.m(tDesign.radius)))

  # center
  center = Float64.(ustrip.(Unitful.m.(tDesign.center)))
  
  return magneticField(coords, field, 
		       R, center, L,
		       x, y, z,
		       calcSolid)
end

function magneticField(coords::AbstractArray{T, 2}, field::Union{AbstractArray{T, 2}, AbstractArray{T, 3}}, 
                       R::T, center::Vector{T}, L::Int,
                       x::Variable, y::Variable, z::Variable, 
                       calcSolid::Bool=true) where T <: Real

  # transpose coords if its dimensions do not fit
  if size(coords,1) != 3
    coords = coords'
  end

  # test dimensions of field array
  if size(field,1) > 3
    throw(DimensionMismatch("The measured field has more than 3 entries in the first dimension: $(size(field,1))"))
  elseif size(field,2) != size(coords,2)
    throw(DimensionMismatch("The field vector does not match the size of the tdesign: $(size(field,2)) != $(size(coords,2))"))
  end

  func= Array{Function}(undef,size(field,1),size(field,3))
  expansion = Array{Polynomial}(undef,size(field,1),size(field,3))
  coeffs = Array{SphericalHarmonicCoefficients}(undef,size(field,1),size(field,3))

  # rescale coordinates to t-design on unit sphere
  coords = coords .- center
  coords *= 1/R
  for c in axes(field,3)
    # calculation of the coefficients
    for j in axes(field,1)
      coeffs[j,c] = SphericalHarmonicExpansions.sphericalQuadrature(field[j,:,c],coords',L);
      coeffs[j,c].R = R

      normalize!(coeffs[j,c],R)

      # convert spherical into solid coefficients
      if calcSolid
          solid!(coeffs[j,c])
      end

      # calculation of the expansion
      expansion[j,c] = sphericalHarmonicsExpansion(coeffs[j,c],x,y,z) + 0*x;
      func[j,c] = @fastfunc expansion[j,c]+0*x+0*y+0*z
    end
  end

  return coeffs, expansion, func
end

#TODO: This should be merged with and moved to MPIFiles
function loadTDesignCoefficients(filename::String)
  field, radius, N, t, center, correction  = h5open(filename, "r") do file
    field = read(file,"/fields") 		# measured field (size: 3 x #points x #patches)
    radius = read(file,"/positions/tDesign/radius")	# radius of the measured ball
    N = read(file,"/positions/tDesign/N")		# number of points of the t-design
    t = read(file,"/positions/tDesign/t")		# t of the t-design
    center = read(file,"/positions/tDesign/center")	# center of the measured ball
    correction = read(file, "/sensor/correctionTranslation")
    return field, radius, N, t, center, correction
  end
  @polyvar x y z
  tDes = MPIFiles.loadTDesign(Int(t),N,radius*u"m", center.*u"m")
  coeffs, = magneticField(tDes, field, x,y,z)
  for c=1:size(coeffs,2), j = 1:3
    coeffs[j, c] = SphericalHarmonicExpansions.translation(coeffs[j, c], correction[:, j])
  end

  coeffs_MF = MagneticFieldCoefficients(coeffs, radius, center)

  return coeffs_MF
end