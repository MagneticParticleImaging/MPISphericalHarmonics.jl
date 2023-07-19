module MPISphericalHarmonics

using DocStringExtensions
using LinearAlgebra
using Unitful

using MPIMagneticFields
using SphericalHarmonicExpansions
using MPIFiles
using NLsolve # for Newton solver

import Base.length

# load MagneticFieldCoefficients
include("MagneticFieldCoefficients.jl")
export MagneticFieldCoefficients
export shift, shift!, shiftFFP!

## SphericalHarmonicsDefinedField ##
export SphericalHarmonicsDefinedField
export selectPatch, length

Base.@kwdef mutable struct SphericalHarmonicsDefinedField <: AbstractMagneticField
  func::Array{SphericalHarmonicExpansions.StaticPolynomials.Polynomial, 2}
  patch::Integer = 1
end

function SphericalHarmonicsDefinedField(filename::String)

  # load coefficients
  mfc = MagneticFieldCoefficients(filename)
  # get field function
  func = fastfunc.(mfc.coeffs)

  return SphericalHarmonicsDefinedField(func=func)
end

# constructors for coefficients
SphericalHarmonicsDefinedField(coeffs::Array{SphericalHarmonicCoefficients}) = SphericalHarmonicsDefinedField(func = fastfunc.(coeffs))
SphericalHarmonicsDefinedField(coeffs_MF::MagneticFieldCoefficients) = SphericalHarmonicsDefinedField(coeffs_MF.coeffs)

MPIMagneticFields.FieldStyle(::SphericalHarmonicsDefinedField) = OtherField()
MPIMagneticFields.FieldDefinitionStyle(::SphericalHarmonicsDefinedField) = SphericalHarmonicsDataBasedFieldDefinition()
MPIMagneticFields.FieldTimeDependencyStyle(::SphericalHarmonicsDefinedField) = TimeConstant()

# get field values
MPIMagneticFields.value_(field::SphericalHarmonicsDefinedField, r) = [field.func[i, field.patch](r) for i=1:3]

# patches
length(field::SphericalHarmonicsDefinedField) = size(field.func,2) # get number of patches
function selectPatch(field::SphericalHarmonicsDefinedField, patchNum::Int) 
  if 1 <= patchNum <= length(field) # test if patch exists
    field.patch = patchNum
  else
    throw(DimensionMismatch("The field contains only $(length(field)) patches."))
  end
end


end
