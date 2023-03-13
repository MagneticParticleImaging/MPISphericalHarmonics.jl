using MPISphericalHarmonics
using MPIMagneticFields
using MPIFiles 
using Test
using Aqua

@testset "MPISphericalHarmonics.jl" begin
  @testset "Aqua" begin
    Aqua.test_all(MPISphericalHarmonics, ambiguities=false)
  end

 
  @testset "Ideal Coefficents" begin

    ## Calculate some field values
    # load a spherical t-design
    tDes = loadTDesign(6,26,42u"mm")
    # ideal gradient field with FFP
    idealField = IdealFFP([-1,-1,2])
    # get field values
    fieldValues = hcat([idealField[ustrip.(Unitful.m.(pos))...] for pos in tDes]...)

    ## Calculate coefficients
    MPISphericalHarmonics.@polyvar x y z
    coeffs, = MPISphericalHarmonics.magneticField(tDes, fieldValues, x,y,z)

    ## Test coefficients
    numCoeffs = length(coeffs[1,1].c)
    # ideal coefficients
    ideal = [zeros(numCoeffs) for j=1:3]
    ideal[1][4] = ideal[2][2] = idealField.gradient[1] # gradient (x,y)
    ideal[3][3] = idealField.gradient[3] # gradient (z)
    for j=1:3
      @test isapprox(coeffs[j,1].c, ideal[j], atol = 1e-10)
    end

  end

end
