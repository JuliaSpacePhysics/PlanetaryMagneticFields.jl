using Test
using PlanetaryMagneticFields
using LinearAlgebra

@testset "load_model" begin
    # Test loading by unique model name (Symbol)
    model_sym = load_model(:JRM09)
    @test isa(model_sym, PlanetaryMagneticFields.MagneticModel)
    @test model_sym.obj.name == :jupiter

    # Test with keyword arguments
    model_truncated = load_model(:jupiter, "jrm33"; max_degree = 13)
    @test degree(model_truncated) == 13

    # Test convenience accessor functions
    model_jrm33 = JRM33(max_degree = 13)
    @test degree(model_jrm33) == 13
end

@testset "callable model" begin
    using Unitful
    model = JRM09()

    # Position at 1.5 RJ
    r_RJ = 1.5
    θ, φ = π / 4, 0.0

    # Evaluate in planetary radii
    B_RJ = model(r_RJ, θ, φ)
    # Same position in km
    r_km = r_RJ * 71492.0u"km"  # Jupiter radius
    @test B_RJ == model(r_km, θ, φ)
    𝐫 = [r_km * sin(θ) * cos(φ), r_km * sin(θ) * sin(φ), r_km * cos(θ)]
    @test B_RJ == model(𝐫)
end

@testset "available_models" begin
    models = available_models(:jupiter)
    @test "jrm09" in models
    @test "jrm33" in models

    models = available_models()
    @test length(models) == length(unique(models))
end

@testset "all planets" begin
    # Test that all planets are supported
    planets = [:mercury, :earth, :mars, :jupiter, :saturn, :uranus, :neptune, :ganymede]
    for planet in planets
        models = available_models(planet)
        @test length(models) > 0
    end

    # Test loading a model for each planet
    test_models = [
        (:mercury, "anderson2012"),
        (:earth, "igrf2020"),
        (:mars, "langlais2019"),
        (:jupiter, "jrm09"),
        (:saturn, "cassini11"),
        (:uranus, "ah5"),
        (:neptune, "gsfco8"),
        (:ganymede, "kivelson2002a"),
    ]

    for (planet, model_name) in test_models
        model = load_model(planet, model_name)
        @test isa(model, PlanetaryMagneticFields.MagneticModel)

        # Test that we can evaluate the field
        B = model(1.0, π / 4, 0.0)
        @test length(B) == 3
        @test all(isfinite, B)
    end
end
