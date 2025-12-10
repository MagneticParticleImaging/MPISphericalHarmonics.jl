@testset "magneticField" begin
  # further tests for function magneticField()
  coords = Float64.(ustrip.(Unitful.m.(hcat([p for p in tDes]...))))
  R = Float64(ustrip(Unitful.m(tDes.radius)))
  center = Float64.(ustrip.(Unitful.m.(tDes.center)))
  L = floor(Int, tDes.T / 2)

  # transposed positions
  coeffsTest = MPISphericalHarmonics.magneticField(coords', fieldValues, R, center, L)
  for j = 1:3
      @test isapprox(coeffs[j, 1].c, coeffsTest[j, 1].c, atol = 1e-10)
  end

  # Errors
  # >3 field values in the first dimension:
  @test_throws DimensionMismatch MPISphericalHarmonics.magneticField(coords, zeros(4, length(tDes)), R, center, L)  
  # number of field values != number of measured positions:
  @test_throws DimensionMismatch MPISphericalHarmonics.magneticField(coords, fieldValues[:, 1:end-1], R, center, L) 
end

@testset "MagneticFieldCoefficients" begin

  ## Test Constructor
  # Errors
  @test_throws DomainError MagneticFieldCoefficients(-2)
  @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs[1:2, :])

  # standard constructors
  coeffsMF = MagneticFieldCoefficients(coeffs)
  @test coeffsMF.coeffs == coeffs
  @test coeffsMF.radius == 0.0
  @test coeffsMF.center == zeros(Float64, 3, 1)
  @test isnothing(coeffsMF.ffp)

  L = 1
  coeffsMF = MagneticFieldCoefficients(L)
  @test size(coeffsMF) == (3, 1)
  [@test coeffsMF.coeffs[j, 1].L == L for j = 1:3]

  # multiple patches
  coeffsMF = MagneticFieldCoefficients(L, numPatches = 2)
  @test size(coeffsMF) == (3, 2)
  @test length(coeffsMF) == 2 # number of patches

  coeffsMF = MagneticFieldCoefficients(coeffs, 0.042, zeros(3, 1))
  @test coeffsMF.coeffs == coeffs
  @test coeffsMF.radius == 0.042
  @test coeffsMF.center == zeros(Float64, 3, 1)
  @test isnothing(coeffsMF.ffp)

  # constructor with wrong sizes of the center
  @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs, 0.042, zeros(2, 1))
  @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs, 0.042, zeros(3, 2))

  # constructor with t-design
  coeffsMF = MagneticFieldCoefficients(coeffs, tDes, zeros(Float64, 3, 1))
  @test coeffsMF.radius == 0.042
  @test coeffsMF.center == zeros(Float64, 3, 1)
  @test coeffsMF.ffp == zeros(Float64, 3, 1)

  # constructor with wrong FFP sizes
  @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs, tDes, zeros(2, 1))
  @test_throws DimensionMismatch MagneticFieldCoefficients(coeffs, tDes, zeros(3, 2))

  # test offset/gradient/Jacobian
  @test isapprox(getOffset(coeffsMF), [0.0, 0.0, 0.0], atol = 1e-10)
  @test isapprox(getGradient(coeffsMF), [[-1.0, -1.0, 2.0]], atol = 1e-10)
  @test isapprox(getJacobian(coeffsMF), [[-1.0 0.0 0.0; 0.0 -1.0 0.0; 0.0 0.0 2.0]], atol = 1e-10)
end