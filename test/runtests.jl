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
    tDes = loadTDesign(6,26,42.0u"mm")
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

    @testset "MagneticFieldCoefficients" begin

      ## Test Constructor
      # Errors
      @test_throws DomainError MagneticFieldCoefficients(-2)
      @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs[1:2,:])

      # standard constructors
      coeffsMF = MagneticFieldCoefficients(coeffs)
      @test coeffsMF.coeffs == coeffs
      @test coeffsMF.radius == 0.0
      @test coeffsMF.center == zeros(Float64,3)
      @test coeffsMF.ffp === nothing

      # constructor with t-design
      coeffsMF = MagneticFieldCoefficients(coeffs, tDes, zeros(Float64,3,1))
      @test coeffsMF.radius == 0.042
      @test coeffsMF.center == zeros(Float64,3)
      @test coeffsMF.ffp == zeros(Float64,3,1) 
    
      # constructor with wrong FFP sizes
      @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs, tDes, zeros(2,1))
      @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs, tDes, zeros(3,2))
    end

    @testset "Load data from file" begin
      ɛ = eps(Float64)

      ## measurement data (without coefficients)
      filename = "idealGradientField.h5"
      func = SphericalHarmonicsDefinedField(filename)
      @test isapprox(func[0.01,0.01,0.01], [-0.01,-0.01,0.02], atol=ε)

      # get coefficients
      coeffsMF, = MPISphericalHarmonics.loadTDesignCoefficients(filename)
      @test isapprox(coeffsMF.radius, 0.042, atol=ε) # radius
      @test isapprox(coeffsMF.coeffs[1][1,1], -1.0, atol=1e-10) # gradient (x)
      @test isapprox(coeffsMF.coeffs[2][1,-1], -1.0, atol=1e-10) # gradient (y)
      @test isapprox(coeffsMF.coeffs[3][1,0], 2.0, atol=1e-10) # gradient (z)

      ## load coefficients from file
      # write coefficients in file
      filename2 = "idealGradientFieldwithCoeffs.h5"
      filename3 = "idealGradientFieldwithCoeffs2.h5"
      filename4 = "idealGradientFieldwithCoeffs3.h5"
      cp(filename, filename2)
      write(filename2, coeffsMF.coeffs)
      # add radius and center
      cp(filename2, filename3)
      h5open(filename3, "cw") do file
        write(file, "/radius", coeffsMF.radius)
        write(file, "/center", coeffsMF.center)
      end
      # add FFP
      cp(filename3, filename4)
      h5open(filename4, "cw") do file
        write(file, "/ffp", zeros(3,1))
      end

      # only coefficients (no further informations given)
      coeffsTest = MagneticFieldCoefficients(filename2)
      @test isapprox(coeffsTest.radius, 0.01, atol=ε) # radius
 
      # with given radius & center
      coeffsTest = MagneticFieldCoefficients(filename3)
      @test isapprox(coeffsTest.radius, 0.042, atol=ε) # radius
      @test isapprox(coeffsTest.center, zeros(3), atol=ε) # center
      @test coeffsTest.ffp === nothing # FFP
      
      # with given FFP
      coeffsTest = MagneticFieldCoefficients(filename4)
      @test coeffsTest.ffp == zeros(3,1) # FFP

      # remove test files
      rm(filename2)
      rm(filename3)
      rm(filename4)
    end

    @testset "Multiple patches" begin
      #TODO
    end
  end

end
