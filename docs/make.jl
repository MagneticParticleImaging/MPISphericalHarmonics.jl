using MPISphericalHarmonics
using Documenter

DocMeta.setdocmeta!(MPISphericalHarmonics, :DocTestSetup, :(using MPISphericalHarmonics); recursive=true)

makedocs(;
    modules=[MPISphericalHarmonics],
    authors="Marija Boberg <m.boberg@uke.de> and contributors",
    repo="https://github.com/MagneticParticleImaging/MPISphericalHarmonics.jl/blob/{commit}{path}#{line}",
    sitename="MPISphericalHarmonics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://github.com/MagneticParticleImaging/MPISphericalHarmonics.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(repo   = "github.com/MagneticParticleImaging/MPISphericalHarmonics.jl.git",
           target = "build")
