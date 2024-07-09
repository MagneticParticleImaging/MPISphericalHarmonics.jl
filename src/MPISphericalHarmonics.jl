module MPISphericalHarmonics

using DocStringExtensions
using LinearAlgebra
using Unitful

using MPIMagneticFields
using SphericalHarmonicExpansions
using MPIFiles
using NLsolve # for Newton solver
using HDF5

import Base.length

# load MagneticFieldCoefficients
include("MagneticFieldCoefficients.jl")
export MagneticFieldCoefficients
export getOffset, getGradient, getJacobian
export shift, shift!, shiftFFP!
export findFFP, findFFP!

## SphericalHarmonicsDefinedField ##
export SphericalHarmonicsDefinedField
export setPatch!, getPatch, length

Base.@kwdef mutable struct SphericalHarmonicsDefinedField <: AbstractMagneticField
    func::Array{Union{Function,SphericalHarmonicExpansions.StaticPolynomials.Polynomial},2}
    patch::Integer = 1
end

function SphericalHarmonicsDefinedField(filename::String)

    # load coefficients
    mfc = MagneticFieldCoefficients(filename)
    # get field function
    func = fastfunc.(mfc.coeffs)

    return SphericalHarmonicsDefinedField(func = func)
end

# constructors for coefficients
SphericalHarmonicsDefinedField(coeffs::Array{SphericalHarmonicCoefficients}) =
    SphericalHarmonicsDefinedField(func = fastfunc.(coeffs))
SphericalHarmonicsDefinedField(coeffs_MF::MagneticFieldCoefficients) =
    SphericalHarmonicsDefinedField(coeffs_MF.coeffs)

MPIMagneticFields.FieldStyle(::SphericalHarmonicsDefinedField) = OtherField()
MPIMagneticFields.FieldDefinitionStyle(::SphericalHarmonicsDefinedField) =
    SphericalHarmonicsDataBasedFieldDefinition()
MPIMagneticFields.FieldTimeDependencyStyle(::SphericalHarmonicsDefinedField) =
    TimeConstant()

# get field values
MPIMagneticFields.value_(field::SphericalHarmonicsDefinedField, r) =
    [field.func[i, field.patch](r) for i = 1:3]

# patches
length(field::SphericalHarmonicsDefinedField) = size(field.func, 2) # get number of patches
function setPatch!(field::SphericalHarmonicsDefinedField, patch::Int) # set patch
    checkPatch(field, patch) # test if patch exists
    field.patch = patch # set patch

    return patch
end
getPatch(field::SphericalHarmonicsDefinedField) = field.patch # get current patch

# test if patch exists
function checkPatch(field::SphericalHarmonicsDefinedField, patch::Int)
    if patch > length(field) || patch < 1 
        throw(DimensionMismatch("Patch $patch requested, but the field contains only $(length(field)) patches."))
    end

    return nothing
end
end
