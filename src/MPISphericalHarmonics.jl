module MPISphericalHarmonics

using DocStringExtensions
using LinearAlgebra
using Unitful

using MPIMagneticFields
using SphericalHarmonicExpansions
using MPIFiles

mutable struct MagneticFieldCoefficients
  coeffs::Array{SphericalHarmonicCoefficients,2}
  radius::Float64
  center::Vector{Float64}
  ffp::Union{Array{Float64,2},Nothing}

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
  return MagneticFieldCoefficients(reshape([SphericalHarmonicCoefficients(L),SphericalHarmonicCoefficients(L),SphericalHarmonicCoefficients(L)],3,1))
end
  
# read coefficients from an HDF5-file
function MagneticFieldCoefficients(path::String)
  file = h5open(path,"r")

  # load spherical harmonic coefficients
  shcoeffs = SphericalHarmonicCoefficients(path)

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
    # use radius = 0.042 as default value
    return MagneticFieldCoefficients(shcoeffs, 0.042)
  end
end

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
		       L::Int=Int(floor(tDesign.T/2)),
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
function loadTDesign(filename::String)
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
  coeffs, expansion, func = magneticField(tDes, field, x,y,z)
  for c=1:size(coeffs,2), j = 1:3
    coeffs[j, c] = SphericalHarmonicExpansions.translation(coeffs[j, c], correction[:, j])
    expansion[j, c] = sphericalHarmonicsExpansion(coeffs[j, c], x, y, z);
    func[j, c] = @fastfunc expansion[j, c]
  end

  coeffs_MF = MagneticFieldCoefficients(coeffs, radius, center)

  return coeffs_MF, expansion, func
end

export SphericalHarmonicsDefinedField
Base.@kwdef mutable struct SphericalHarmonicsDefinedField <: AbstractMagneticField
  func::Array{Function, 2}
  patch::Integer = 1
end

function SphericalHarmonicsDefinedField(filename::String)
  coeffs_MF, expansion, func = loadTDesign(filename)
  return SphericalHarmonicsDefinedField(func=func)
end

MPIMagneticFields.fieldType(::SphericalHarmonicsDefinedField) = OtherField()
MPIMagneticFields.definitionType(::SphericalHarmonicsDefinedField) = SphericalHarmonicsDataBasedFieldDefinition()
MPIMagneticFields.timeDependencyType(::SphericalHarmonicsDefinedField) = TimeConstant()

MPIMagneticFields.value(field::SphericalHarmonicsDefinedField, r::PT) where {T <: Number, PT <: AbstractVector{T}} = [field.func[i, field.patch].(r...) for i=1:3]

selectPatch(field::SphericalHarmonicsDefinedField, patchNum) = field.patch = patchNum

end