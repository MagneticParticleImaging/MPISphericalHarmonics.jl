@testset "SphericalHarmonicsDefinedField (multiple patches)" begin
  ## Multi-patch setting: Second field with offset
  coeffsPatch = hcat(deepcopy(coeffs), deepcopy(coeffs)) # two patches
  for j = 1:3
      coeffsPatch[j, 2][0, 0] = 0.01
  end # set offset
  coeffsPatch_MF = MagneticFieldCoefficients(coeffsPatch) # create MagneticFieldCoefficients
  field = SphericalHarmonicsDefinedField(coeffsPatch_MF) # test constructor on coefficients
  fieldSP = SphericalHarmonicsDefinedField(coeffsPatch_MF[1]) # test constructor on coefficients of a single patch

  # Test field types
  @test FieldStyle(field) isa OtherField
  @test FieldDefinitionStyle(field) isa SphericalHarmonicsDataBasedFieldDefinition
  @test FieldTimeDependencyStyle(field) isa TimeConstant

  # Test number of patches
  @test length(field) == 2
  @test length(fieldSP) == 1

  ## Test FFPs (for both patches)
  # First patch
  @test field.patch == 1
  ffp = zeros(3)
  @test isapprox(field[ffp...], zeros(3), atol = 1e-10)

  # Second patch
  setPatch!(field, 2)
  @test field.patch == 2 # tests setPatch!
  @test getPatch(field) == 2 # tests getPatch
  # test offset in (0, 0, 0)
  offset = ones(3) .* 0.01
  @test isapprox(field[0, 0, 0], offset, atol = 1e-10)
  # test FFP
  ffp = [0.01, 0.01, -0.005]
  @test isapprox(field[ffp...], zeros(3), atol = 1e-10)

  # Test wrong patch number
  @test_throws DimensionMismatch setPatch!(field, 3)
end