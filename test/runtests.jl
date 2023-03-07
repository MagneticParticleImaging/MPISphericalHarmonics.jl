using MPISphericalHarmonics
using Test
using Aqua

@testset "MPISphericalHarmonics.jl" begin
  @testset "Aqua" begin
    Aqua.test_all(MPISphericalHarmonics, ambiguities=false)
  end


end
