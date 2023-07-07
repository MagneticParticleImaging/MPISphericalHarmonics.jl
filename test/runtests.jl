using MPISphericalHarmonics
using MPIMagneticFields
using MPIFiles 
using SphericalHarmonicExpansions
using Test
using Aqua

@testset "MPISphericalHarmonics.jl" begin
  @testset "Aqua" begin
    Aqua.test_all(MPISphericalHarmonics, ambiguities=false)
  end

 
  @testset "Ideal Coefficents" begin

    ## Calculate some field values
    # load a spherical t-design
    tDes = loadTDesign(6, 26, 42.0u"mm")
    # ideal gradient field with FFP
    idealField = IdealFFP([-1, -1, 2])
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

    @testset "magneticField" begin
      # further tests for function magneticField()
      coords = Float64.(ustrip.(Unitful.m.(hcat([p for p in tDes]...))))
      R = Float64(ustrip(Unitful.m(tDes.radius)))
      center = Float64.(ustrip.(Unitful.m.(tDes.center)))
      L = floor(Int,tDes.T/2)
  
      # transposed positions
      coeffsTest, = MPISphericalHarmonics.magneticField(coords', fieldValues, R, center, L, x, y, z)
      for j=1:3
        @test isapprox(coeffs[j,1].c, coeffsTest[j,1].c, atol = 1e-10)
      end

      # Errors
      @test_throws DimensionMismatch MPISphericalHarmonics.magneticField(coords, zeros(4,length(tDes)), R, center, L, x, y, z) # >3 field values in the first dimension 
      @test_throws DimensionMismatch MPISphericalHarmonics.magneticField(coords, fieldValues[:,1:end-1], R, center, L, x, y, z) # number of field values != number of measured positions
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
      @test coeffsMF.center == zeros(Float64,3,1)
      @test isnothing(coeffsMF.ffp)

      L = 1
      coeffsMF = MagneticFieldCoefficients(L)
      @test size(coeffsMF.coeffs) == (3,1)
      [@test coeffsMF.coeffs[j,1].L == L for j=1:3]

      coeffsMF = MagneticFieldCoefficients(coeffs,0.042,zeros(3,1))
      @test coeffsMF.coeffs == coeffs
      @test coeffsMF.radius == 0.042
      @test coeffsMF.center == zeros(Float64,3,1)
      @test isnothing(coeffsMF.ffp)

      # constructor with t-design
      coeffsMF = MagneticFieldCoefficients(coeffs, tDes, zeros(Float64,3,1))
      @test coeffsMF.radius == 0.042
      @test coeffsMF.center == zeros(Float64,3,1)
      @test coeffsMF.ffp == zeros(Float64,3,1) 
    
      # constructor with wrong FFP sizes
      @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs, tDes, zeros(2,1))
      @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs, tDes, zeros(3,2))

      # test size
      @test size(coeffsMF) == (3,1)
    end

    @testset "MagneticFieldCoefficient operations" begin
      # coefficients for the tests
      shc0 = fill(SphericalHarmonicCoefficients(zeros(4)), (3,1))
      shc1 = fill(SphericalHarmonicCoefficients(ones(4)), (3,1))
      shc2 = 2 .* shc1
      c0 = MagneticFieldCoefficients(shc0,0.042,zeros(3)) 
      c1 = MagneticFieldCoefficients(shc1,0.042,zeros(3)) 
      c2 = MagneticFieldCoefficients(shc2,0.042,zeros(3)) 
      c1R = MagneticFieldCoefficients(shc1,0.01,zeros(3)) # different radius
      c1C = MagneticFieldCoefficients(shc1,0.042,ones(3)) # different center
      c1F = MagneticFieldCoefficients(shc1,0.042,zeros(3),zeros(3,1)) # FFP 1
      c1F2 = MagneticFieldCoefficients(shc1,0.042,zeros(3),ones(3,1)) # FFP 2

      # isapprox 
      @test !isapprox(c1, c2) # wrong coefficients
      @test !isapprox(c1, c1R) # wrong radius
      @test !isapprox(c1, c1C) # wrong center
      @test !isapprox(c1, c1F) # one FFP
      @test !isapprox(c1F, c1F2) # wrong FFP

      # addition/subtraction # @test_throws DomainError mfc1+mfc2
      @test_throws DomainError c1 + c1R # DomainError radius
      @test_throws DomainError c2 - c1R # DomainError radius
      @test_throws DomainError c1 + c1C # DomainError center
      @test_throws DomainError c2 - c1C # DomainError center
      @test isapprox(+(c1,c1R,force=true), c2) # force = true 
      @test isapprox(-(c2,c1R,force=true), c1) # force = true 
      @test isapprox(c1 + c1, c2) # correct
      @test isapprox(c2 - c1, c1) # correct

      # operations MFC and value
      @test isapprox(1 + c1, c2)
      @test isapprox(c0 + 2, c2)
      @test isapprox(1 - c1, c0)
      @test isapprox(c2 - 2, c0)
      @test isapprox(0 * c1, c0)
      @test isapprox(c1 * 2, c2)
      @test isapprox(c2 / 2, c1)

    end

    @testset "Load/write data from/to file" begin
      ɛ = eps(Float64)

      ## measurement data (without coefficients)
      filename = "idealGradientField.h5"
      field = SphericalHarmonicsDefinedField(filename)
      @test isapprox(field[0.01,0.01,0.01], [-0.01,-0.01,0.02], atol=ε)

      # get coefficients
      coeffsMF = MPISphericalHarmonics.loadTDesignCoefficients(filename)
      @test isapprox(coeffsMF.radius, 0.042, atol=ε) # radius
      @test isapprox(coeffsMF.coeffs[1][1,1], -1.0, atol=1e-10) # gradient (x)
      @test isapprox(coeffsMF.coeffs[2][1,-1], -1.0, atol=1e-10) # gradient (y)
      @test isapprox(coeffsMF.coeffs[3][1,0], 2.0, atol=1e-10) # gradient (z)

      ## load coefficients from file
      # write coefficients in file
      filename2 = "Coeffs.h5"
      filename3 = "Coeffs2.h5"
      filename4 = "Coeffs3.h5"
      filenameW = "CoeffsW.h5"
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
      field = SphericalHarmonicsDefinedField(filename4)
      @test isapprox(field[-0.02,0.02,-0.03], [0.02,-0.02,-0.06], atol=ε)

      # test write
      MPISphericalHarmonics.write(filenameW, coeffsTest)
      coeffsW = MagneticFieldCoefficients(filenameW)
      @test isapprox(coeffsW.radius, 0.042, atol=ε) # radius
      @test isapprox(coeffsW.center, zeros(3), atol=ε) # center
      @test coeffsW.ffp == zeros(3,1) # FFP
      coeffsW = nothing # This maybe masks an implementation error
      GC.gc()

      # remove test files
      rm(filename2)
      rm(filename3)
      rm(filename4)
      rm(filenameW)
    end

    @testset "SphericalHarmonicsDefinedField (multiple patches)" begin
      ## Multi-patch setting: Second field with offset
      coeffsPatch = hcat(deepcopy(coeffs),deepcopy(coeffs)) # two patches
      for j=1:3 coeffsPatch[j,2][0,0] = 0.01 end # set offset
      coeffsPatch_MF = MagneticFieldCoefficients(coeffsPatch) # create MagneticFieldCoefficients
      field = SphericalHarmonicsDefinedField(coeffsPatch_MF) # test constructor on coefficients

      # Test field types
      @test FieldStyle(field) isa OtherField
      @test FieldDefinitionStyle(field) isa SphericalHarmonicsDataBasedFieldDefinition
      @test FieldTimeDependencyStyle(field) isa TimeConstant

      # Test number of patches
      @test length(field) == 2

      ## Test FFPs (for both patches)
      # First patch
      @test field.patch == 1
      ffp = zeros(3)
      @test isapprox(field[ffp...], zeros(3), atol=1e-10)

      # Second patch
      selectPatch(field,2)   
      @test field.patch == 2
      # test offset in (0, 0, 0)
      offset = ones(3) .* 0.01
      @test isapprox(field[0, 0, 0], offset, atol=1e-10)
      # test FFP
      ffp = [0.01, 0.01, -0.005]
      @test isapprox(field[ffp...], zeros(3), atol=1e-10)

      # Test wrong patch number
      @test_throws DimensionMismatch selectPatch(field,3)
    end
  end

end
