###################################################
# Example:
# 1. Create some ideal simulated data (2 patches)
# 2. Calculate coefficients from the data
#  2.1 Save coefficients
# 3. Calculate FFP and shift the coefficients
# 4. Plot the fields
################################################### 
# Julia: 1.10

# Activate local environment
using Pkg
Pkg.activate(joinpath(@__DIR__,"."))
Pkg.instantiate()

# Packages
using MPISphericalHarmonics
using MPISphericalHarmonics.MPIFiles, MPISphericalHarmonics.Unitful # load t-design
using MPISphericalHarmonics.MPIMagneticFields # simulation of ideal fields

# Plotting
using PyPlot
using LinearAlgebra

## 1. Simulate measurement data of an ideal field ##
tDes = loadTDesign(6, 26, 40.0u"mm") # spherical t-design (measurement points)
idealField = [IdealFFP([-1, -1, 2]),   # create an ideal field with field-free point
	IdealFFP([-1, -1, 2]) + IdealHomogeneousField([0.01,0.02,0.0])] # ideal field with shifted FFP
fieldValues = cat([hcat([idealField[p][ustrip.(Unitful.m.(pos))...] for pos in tDes]...) for p=1:2]...,dims=3) # get the field values (measurement data)

## 2. Calculate the coefficients ##
coeffs_ = MPISphericalHarmonics.magneticField(tDes, fieldValues) # SphericalHarmonicCoefficients
coeffs = MagneticFieldCoefficients(coeffs_, ustrip(Unitful.m(tDes.radius)))

## 2.1 optional: write coefficients in a file ##
filename = "coeffsExample.h5"
write(filename, coeffs)
# coeffs = MagneticFieldCoefficients(filename) # load coefficients from file

## 3. Calculate FFP and shift coefficients into the FFP ##
@info "Initial center:" coeffs.center
shiftFFP!(coeffs)
@info "Shifted center:" coeffs.center

## 4. Plotting ##
# load field function
field = SphericalHarmonicsDefinedField(coeffs)
# field = SphericalHarmonicsDefinedField(filename) # load field from coefficients-file

# range of the plots
plotRange = -coeffs.radius:0.001:coeffs.radius
quiverRange = -coeffs.radius:0.005:coeffs.radius

# plot
figure(figsize=(12,6))
for p=1:length(coeffs) # plot both patches
  subplot(1,2,p)
  selectPatch(field,p) # switch patches

  # plot norm
  imshow(norm.(field[plotRange, plotRange, 0]), origin="lower",
		extent = (plotRange[1], plotRange[end], 
			  plotRange[1], plotRange[end])),
  # plot arrows
  quiver(quiverRange, quiverRange, 
	 getindex.(field[quiverRange, quiverRange, 0],2), 
	 getindex.(field[quiverRange, quiverRange, 0],1))

  # draw circle of the measurement 
  plot(coeffs.radius .* sin.(0:0.01:2*pi) .+ coeffs.center[1,p], 
       coeffs.radius .* cos.(0:0.01:2*pi) .+ coeffs.center[2,p], 
       color="white")

  # labels,...
  xlim((-1,1).*coeffs.radius), ylim((-1,1).*coeffs.radius)
  xlabel("x / m"), ylabel("y / m")
  title("Patch $p")
end
suptitle("Fields shifted into their FFP")
