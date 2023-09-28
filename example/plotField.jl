using PyPlot
using LinearAlgebra
using Artifacts

using MPISphericalHarmonics

filename = joinpath(artifact"example", "FFP000COMSOLOptimized_tDesgn_t=12_r=45.h5")
field = SphericalHarmonicsDefinedField(filename)

figure()
imshow(norm.(field[-0.02:0.001:0.02, -0.02:0.001:0.02, 0]))