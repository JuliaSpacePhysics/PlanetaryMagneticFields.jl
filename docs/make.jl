using PlanetaryMagneticFields
using Documenter

DocMeta.setdocmeta!(PlanetaryMagneticFields, :DocTestSetup, :(using PlanetaryMagneticFields); recursive = true)

makedocs(;
    modules = [PlanetaryMagneticFields],
    authors = "Beforerr <zzj956959688@gmail.com> and contributors",
    sitename = "PlanetaryMagneticFields.jl",
    format = Documenter.HTML(;
        canonical = "https://JuliaSpacePhysics.github.io/PlanetaryMagneticFields.jl",
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Examples" => "examples.md",
    ],
)

deploydocs(;
    repo = "github.com/JuliaSpacePhysics/PlanetaryMagneticFields.jl",
    devbranch = "main",
    push_preview = true,
)
