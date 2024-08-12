@testset "Load/write data from/to file" begin
  ɛ = eps(Float64)

  ## Create measurement files
  filename = joinpath(tmpdir, "idealGradientField.h5")
  filenameV1 = joinpath(tmpdir, "idealGradientFieldV1.h5") # old file structure

  # write data (v2)
  h5open(filename, "w") do file
      write(file,"/fields", fieldValues) 		# measured field (size: 3 x #points x #patches)
      write(file,"/positions/tDesign/radius", ustrip(u"m", tDes.radius))	# radius of the measured ball
      write(file,"/positions/tDesign/N", size(tDes.positions,2))		# number of points of the t-design
      write(file,"/positions/tDesign/t", tDes.T)		# t of the t-design
      write(file,"/positions/tDesign/center", ustrip.(u"m", tDes.center))	# center of the measured ball
      write(file, "/sensor/correctionTranslation", zeros(3,3)) 
  end
  # write data (v1)
  h5open(filenameV1, "w") do file
      write(file,"fields", fieldValues) 		# measured field (size: 3 x #points x #patches)
      write(file,"fieldsError", zeros(Float64, size(fieldValues))) # error of the measured field
      write(file,"positions", tDes.positions)	# measured positions
      write(file,"positionsTDesignRadius", ustrip(u"m", tDes.radius))	# radius of the measured ball
      write(file,"positionsTDesignN", size(tDes.positions,2))		# number of points of the t-design
      write(file,"positionsTDesignT", tDes.T)		# t of the t-design
      write(file,"positionsCenter", ustrip.(u"m", tDes.center))	# center of the measured ball
      write(file,"currents",ones(3, 1)) # currents used for the measurement
  end

  ## measurement data (without coefficients, v1)
  field = SphericalHarmonicsDefinedField(filenameV1)
  @test isapprox(field[0.01, 0.01, 0.01], [-0.01, -0.01, 0.02], atol = ε)

  ## measurement data (without coefficients)
  field = SphericalHarmonicsDefinedField(filename)
  @test isapprox(field[0.01, 0.01, 0.01], [-0.01, -0.01, 0.02], atol = ε)

  # get coefficients
  coeffsMF = MagneticFieldCoefficients(filename)
  @test isapprox(coeffsMF.radius, 0.042, atol = ε) # radius
  @test isapprox(coeffsMF.coeffs[1][1, 1], -1.0, atol = 1e-10) # gradient (x)
  @test isapprox(coeffsMF.coeffs[2][1, -1], -1.0, atol = 1e-10) # gradient (y)
  @test isapprox(coeffsMF.coeffs[3][1, 0], 2.0, atol = 1e-10) # gradient (z)

  ## load coefficients from file
  # write coefficients in file
  filename2 = joinpath(tmpdir, "Coeffs.h5")
  filename3 = joinpath(tmpdir, "Coeffs2.h5")
  filename4 = joinpath(tmpdir, "Coeffs3.h5")
  filenameW = joinpath(tmpdir, "CoeffsW.h5")
  filenameE = joinpath(tmpdir, "CoeffsE.h5")
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
      write(file, "/ffp", zeros(3, 1))
  end

  # only coefficients (no further informations given)
  coeffsTest = MagneticFieldCoefficients(filename2)
  @test isapprox(coeffsTest.radius, 0.0, atol = ε) # radius

  # with given radius & center
  coeffsTest = MagneticFieldCoefficients(filename3)
  @test isapprox(coeffsTest.radius, 0.042, atol = ε) # radius
  @test isapprox(coeffsTest.center, zeros(3), atol = ε) # center
  @test coeffsTest.ffp === nothing # FFP

  # with given FFP
  coeffsTest = MagneticFieldCoefficients(filename4)
  @test coeffsTest.ffp == zeros(3, 1) # FFP
  field = SphericalHarmonicsDefinedField(filename4)
  @test isapprox(field[-0.02, 0.02, -0.03], [0.02, -0.02, -0.06], atol = ε)

  # test write
  MPISphericalHarmonics.write(filenameW, coeffsTest)
  coeffsW = MagneticFieldCoefficients(filenameW)
  @test isapprox(coeffsW.radius, 0.042, atol = ε) # radius
  @test isapprox(coeffsW.center, zeros(3), atol = ε) # center
  @test coeffsW.ffp == zeros(3, 1) # FFP
  coeffsW = nothing # This maybe masks an implementation error

  # test error
  h5write(filenameE, "test", 0)
  @test_throws ErrorException MagneticFieldCoefficients(filenameE)

  GC.gc()

  # remove test files
  rm(filename2)
  rm(filename3)
  rm(filename4)
  rm(filenameW)
  rm(filenameE)

  # remove measurement files
  rm(filename)
  rm(filenameV1)
end