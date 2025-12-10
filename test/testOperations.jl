@testset "MagneticFieldCoefficient operations" begin
  # coefficients for the tests
  shc0 = fill(SphericalHarmonicCoefficients(zeros(4)), (3, 1))
  shc1 = fill(SphericalHarmonicCoefficients(ones(4)), (3, 1))
  shc2 = 2 .* shc1
  c0 = MagneticFieldCoefficients(shc0, 0.042, zeros(3))
  c1 = MagneticFieldCoefficients(shc1, 0.042, zeros(3))
  c2 = MagneticFieldCoefficients(shc2, 0.042, zeros(3))
  c1R = MagneticFieldCoefficients(shc1, 0.01, zeros(3)) # different radius
  c1C = MagneticFieldCoefficients(shc1, 0.042, ones(3)) # different center
  c1F = MagneticFieldCoefficients(shc1, 0.042, zeros(3), zeros(3, 1)) # FFP 1
  c1F2 = MagneticFieldCoefficients(shc1, 0.042, zeros(3), ones(3, 1)) # FFP 2

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
  @test isapprox(+(c1, c1R, force = true), c2) # force = true 
  @test isapprox(-(c2, c1R, force = true), c1) # force = true 
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

  # indexing
  # setup MagneticFieldCoefficients with multiple patches
  csh = reshape(hcat(SphericalHarmonicCoefficients.([rand(16) for i = 1:12])...), 3, 4)
  center = rand(3, 4)
  ffp = rand(3, 4)
  coeffsMF = MagneticFieldCoefficients(csh, 0.01, center, ffp)
  c1i = MagneticFieldCoefficients(csh[:, 2:2], 0.01, center[:, 2:2], ffp[:, 2:2])
  c2i = MagneticFieldCoefficients(csh[:, 3:4], 0.01, center[:, 3:4], ffp[:, 3:4])

  # test getindex
  @test isapprox(coeffsMF[2], c1i)
  @test isapprox(coeffsMF[3:4], c2i)

  # test setindex
  coeffsMF[1] = coeffsMF[2]
  @test isapprox(coeffsMF[1], c1i)
  coeffsMF[1:2] = coeffsMF[3:4]
  @test isapprox(coeffsMF[1:2], c2i)

  # test iterator
  for (i,c) in enumerate(coeffsMF)
      @test coeffsMF[i] == c
  end

  # test element type
  @test eltype(coeffsMF) == MagneticFieldCoefficients

  # test hash
  @test hash(coeffsMF[3:4]) == hash(c2i) # same hash
  @test hash(coeffsMF) !== hash(c1i) # different hash

  # test hcat
  ccat = hcat(c1F, c1F2) # coefficients with FFPs
      @test isapprox(ccat.ffp, [0 1; 0 1; 0 1]) # test FFP
      @test isapprox(ccat.center, zeros(3,2)) # test center
      @test isapprox(ccat.radius, 0.042) # test radius
      @test all(isapprox.(ccat.coeffs, hcat(shc1, shc1))) # test coefficients
  ccat = hcat(c0, c1, c2) # coefficients without FFPs
      @test isnothing(ccat.ffp) # test FFP
      @test length(ccat) == 3 # test number of patches
  # errors
      @test_throws DomainError hcat(c1, c1R) # different radius
      @test_throws DomainError hcat(c1, c1F) # one with FFP, one without FFP
  # force concatenation
  ccat = hcat(c1, c1R, force = true) # different radius
      @test length(ccat) == 2 # test number of patches
      @test isapprox(ccat.radius, c1.radius) # radius of first coefficients used
  ccat = hcat(c1F, c1, force = true) # one with FFP, one without FFP
      @test isapprox(ccat[2], c1) # test coefficients without FFP
      @test isnothing(ccat.ffp) # no FFP

  # test errors 
  c3 = MagneticFieldCoefficients(csh[:, 2:2], 0.02, center[:, 2:2], ffp[:, 2:2]) # different radius
  c4 = MagneticFieldCoefficients(csh[:, 2:2], 0.01, center[:, 2:2]) # no FFP

  @test_throws DomainError coeffsMF[1] = c3 # different radius
  @test_throws DomainError coeffsMF[1] = c4 # no FFP
end