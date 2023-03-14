using PyPlot
using LinearAlgebra

using MPISphericalHarmonics

field = SphericalHarmonicsDefinedField(".\\MPISphericalHarmonics.jl\\FFP000COMSOLOptimized_tDesgn_t=12_r=45.h5")
imshow(norm.(field[-0.02:0.001:0.02, -0.02:0.001:0.02, 0]))