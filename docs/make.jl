using MagneticModels
using Documenter

DocMeta.setdocmeta!(MagneticModels, :DocTestSetup, :(using MagneticModels); recursive=true)

makedocs(;
    modules=[MagneticModels],
    authors="Beforerr <zzj956959688@gmail.com> and contributors",
    sitename="MagneticModels.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaSpacePhysics.github.io/MagneticModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaSpacePhysics/MagneticModels.jl",
    devbranch="main",
)
