@testset "Find FFP" begin

  ## 3D with FFP ##
  # setup an ideal gradient field with FFP
  coeffsMF = MagneticFieldCoefficients( 
              reshape(
                  SphericalHarmonicCoefficients(
                      [[0.2, 0.0, 0.0, -1.0], # x-field
                      [0.2, -1.0, 0.0, 0.0], # y-field
                      [0.2, 0.0, 2.0, 0.0]], # z-field
                      ones(3), # normalized
                      BitVector([1,1,1])), # solid coefficients
                  3,1), # get a matrix
              0.5) # radius

  # find FFP (and set it as coeffsMF.ffp)
  findFFP!(coeffsMF)
  @test isapprox(coeffsMF.ffp, [0.2,0.2,-0.1]) # correct FFP

  # returnasmatrix = false
  ffp = findFFP(coeffsMF, returnasmatrix = false)
  @test isapprox(getproperty.(ffp, :zero), [[0.2,0.2,-0.1]])


  ## 2D with different FFLs ##
  # setup ideal FFLs
  coeffsFFL = MagneticFieldCoefficients(1,1.0,true, radius=0.5, numPatches=4)
  # FFL in z-direction (patch 1)
      coeffsFFL.coeffs[1,1].c = [0.2, 0.0, 0.0, -1.0]
      coeffsFFL.coeffs[2,1].c = [0.2, 1.0, 0.0, 0.0]
      coeffsFFL.coeffs[3,1].c = [0.0, 0.0, 0.0, 0.0]
  # FFL in x-direction (patch 2)
      coeffsFFL.coeffs[1,2].c = [0.0, 0.0, 0.0, 0.0]
      coeffsFFL.coeffs[2,2].c = [0.2, 1.0, 0.0, 0.0]
      coeffsFFL.coeffs[3,2].c = [0.2, 0.0, -1.0, 0.0]
  # FFL in y-direction (patch 3)
      coeffsFFL.coeffs[1,3].c = [0.2, 0.0, 1.0, 0.0]
      coeffsFFL.coeffs[2,3].c = [0.0, 0.0, 0.0, 0.0]
      coeffsFFL.coeffs[3,3].c = [0.2, 0.0, 0.0, 1.0]
  # FFL in y-direction (patch 4)
      coeffsFFL.coeffs[1,4].c = [0.2, 0.0, -2.0, 0.0]
      coeffsFFL.coeffs[2,4].c = [0.0, 0.0, 0.0, 0.0]
      coeffsFFL.coeffs[3,4].c = [0.2, 0.0, 0.0, -2.0]

  # correct FFPs for start=[0,0,0]
  ffpCorrect = [[0.2,-0.2,0.0],
                [0.0,-0.2,0.2],
                [-0.2,0.0,-0.2],
                [0.1,0.0,0.1]]
  # correct FFPs for patch 3 & 4 for start = [0,0.1,0]
  ffpCorrectStart = [-0.2 0.1 -0.2;
                      0.1 0.1 0.1]'

  # find FFP
      # FFL in z-direction
      ffp = findFFP(coeffsFFL[1], vol=:xy) # use Symbol
      @test isapprox(ffp, ffpCorrect[1])

      # FFL in x-direction
      ffp = findFFP(coeffsFFL[2], vol=MPISphericalHarmonics.yz) # use Enum
      @test isapprox(ffp, ffpCorrect[2])

      # FFL in y-direction
          # find FFP for patch 3
          ffp = findFFP(coeffsFFL[3], vol="xz") # use String
          @test isapprox(ffp, ffpCorrect[3])
          # find FFP for patch 3 & 4
          # with start vector (with length(start) != numPatches)
          ffp = findFFP(coeffsFFL[3:4], vol=:xz, start=[0.0,0.1,0.0])
          @test isapprox(ffp, ffpCorrectStart)

  ## Errors
  @test_throws DimensionMismatch findFFP(coeffsFFL[3:4], vol=:xz, start=[0.0,0.1]) # length(start) != 3
  @test_throws DimensionMismatch findFFP(coeffsFFL[3:4], vol=:xz, start=[zeros(3) for i=1:3]) # number of start vectors != numPatches
  @test_logs (:warn,"Volume IBI not defined. Default 3D volume used.")
          findFFP(coeffsMF, vol = "IBI") # volume not defined

  @testset "Shift Coefficients" begin
      # shift coeffficients by any vector or matrix
      coeffsMF2 = hcat(coeffsMF, coeffsMF) # two patches
      # shift coeffs by vector
          shift!(coeffsMF2, [0.1,-0.1,0.1])
          @test isapprox(coeffsMF2.center, -[0.1,-0.1,0.1] .* ones(3,2)) # center should be shifted by the negative shift
          @test isapprox(getOffset(coeffsMF2), [0.1,0.3,0.4] .* ones(3,2)) # new offset due to shift
      # shift coeffs by matrix
          shift!(coeffsMF2, [-0.2 0.1; -0.2 -0.1; -0.2 0.1])
          @test isapprox(coeffsMF2.center, -[-0.1 0.2; -0.3 -0.2; -0.1 0.2]) # center should be shifted
          @test isapprox(getOffset(coeffsMF2), [0.3 0.0; 0.5 0.4; 0.0 0.6]) # new offset due to shift
      # shift coeffs by vector and do not overwrite
          coeffsShift = shift(coeffsMF,[0.1,-0.1,0.1])
          # test shifted coefficients
          @test isapprox(coeffsShift.center, -[0.1,-0.1,0.1]) # center should be shifted
          @test isapprox(getOffset(coeffsShift), [0.1,0.3,0.4]) # new offset due to shift
          # test original coefficients (no change)
          @test isapprox(coeffsMF.center, zeros(3,1)) # center
          @test isapprox(getOffset(coeffsMF), 0.2 .* ones(3,1)) # offset

      # shift coefficients into their FFP
          coeffsMF.ffp = nothing # remove already calculated FFP
          shiftFFP!(coeffsMF) # shift into the FFP including the calculation
          @test isapprox(coeffsMF.ffp, [0.0,0.0,0.0]) # FFP should be in the origin
          @test isapprox(coeffsMF.center, -[0.2,0.2,-0.1]) # center should be -FFP
          @test isapprox(getOffset(coeffsMF), [0.0,0.0,0.0]) # there should be no offset in the FFP

      # test error (shift matrix with wrong size)
          @test_throws DimensionMismatch shift!(coeffsMF, [0.1,-0.1]) # wrong length of vector
          @test_throws DimensionMismatch shift!(coeffsMF, ones(3,2)) # wrong size of matrix (coeffs have only one patch)
  end
end