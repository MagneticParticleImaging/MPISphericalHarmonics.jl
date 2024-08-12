using MPISphericalHarmonics
using MPIMagneticFields
using MPIFiles
using SphericalHarmonicExpansions
using Test
using Aqua
using Scratch
using Unitful
using HDF5

const tmpdir = @get_scratch!("tmp")
@info "If you want to check the output of the tests, please head to $tmpdir."

## Calculate some field values
# load a spherical t-design
tDes = loadTDesign(6, 26, 42.0u"mm")
# ideal gradient field with FFP
idealField = IdealFFP([-1, -1, 2])
# get field values
fieldValues = hcat([idealField[ustrip.(Unitful.m.(pos))...] for pos in tDes]...)

## Calculate coefficients
coeffs = MPISphericalHarmonics.magneticField(tDes, fieldValues)

## Run tests
@testset "MPISphericalHarmonics.jl" begin
    @testset "Aqua" begin
        Aqua.test_all(MPISphericalHarmonics, ambiguities = false)
    end

    @testset "Ideal Coefficents" begin

        ## Test coefficients
        numCoeffs = length(coeffs[1, 1].c)
        # ideal coefficients
        ideal = [zeros(numCoeffs) for j = 1:3]
        ideal[1][4] = ideal[2][2] = idealField.gradient[1] # gradient (x,y)
        ideal[3][3] = idealField.gradient[3] # gradient (z)
        for j = 1:3
            @test isapprox(coeffs[j, 1].c, ideal[j], atol = 1e-10)
        end

        # test magneticField() and MagneticFieldCoefficients()
        include("testBasics.jl")

        # test MagneticFieldCoefficient operations (e.g., +, *, length, indexing, ...)
        include("testOperations.jl")

        # test load/write data from/to file
        include("testWriteLoad.jl")

        # test FFP stuff and shift operations
        include("testFFPandShift.jl")

        # test multiple patches 
        include("testMultiPatch.jl")
    end

end
