using Test
using PlanetaryMagneticFields
using PlanetaryMagneticFields: SphericalHarmonicModel, GaussCoefficients
using LinearAlgebra

@testset "SphericalHarmonicModel construction" begin
    g = zeros(3, 3)
    h = zeros(3, 3)
    g[2, 1] = 410244.7
    coeffs = GaussCoefficients(g, h)
    model = SphericalHarmonicModel("test", coeffs)
    @test degree(model) == 2
    @test order(model) == 2
end

@testset "Cartesian vs Spherical evaluation consistency" begin
    # Load a real model
    g = zeros(3, 3)
    h = zeros(3, 3)
    g[2, 1] = 410244.7
    g[2, 2] = -71498.3
    h[2, 2] = 21330.5
    coeffs = GaussCoefficients(g, h)
    model = SphericalHarmonicModel("test", coeffs)

    # Test position
    x, y, z = 1.5, 0.5, 1.0  # Cartesian (in planetary radii)
    (r, θ, φ) = PlanetaryMagneticFields.car2sph(x, y, z)

    # Evaluate both ways
    B_from_cart = model(x, y, z; in = :cartesian)
    B_sph = PlanetaryMagneticFields.evalsph(model.coeffs, r, θ, φ)
    B_from_sph = PlanetaryMagneticFields.sph2car(
        B_sph[1], B_sph[2], B_sph[3], 1, θ, φ
    )
    @test B_from_cart ≈ B_from_sph rtol = 1.0e-10
end

@testset "Pure dipole field (analytical test)" begin
    # Create a pure dipole model (only g_1^0 coefficient)
    g = zeros(2, 2)
    h = zeros(2, 2)
    g[2, 1] = 428000.0  # g_1^0 in nT (approximate Jupiter dipole)

    coeffs = GaussCoefficients(g, h)
    model = SphericalHarmonicModel("dipole", coeffs)

    # Test at equator (θ = π/2)
    r, θ, φ = 2.0, π / 2, 0.0  # 2 planetary radii, equator
    B = model(r, θ, φ)

    # For a pure dipole at the equator:
    # B_r = 2 * g_1^0 * (a/r)^3 * cos(θ) = 0 (since cos(π/2) = 0)
    # B_θ = g_1^0 * (a/r)^3 * sin(θ) = g_1^0 / r^3 (since sin(π/2) = 1)
    # B_φ = 0

    expected_Bθ = g[2, 1] / (r^3)
    @test abs(B[1]) < 1.0e-10  # B_r ≈ 0
    @test B[2] ≈ expected_Bθ rtol = 1.0e-6
    @test abs(B[3]) < 1.0e-10  # B_φ ≈ 0

    # Test at north pole (θ = 0)
    r, θ, φ = 2.0, 0.0, 0.0
    B = model(r, θ, φ)

    # At the pole:
    # B_r = 2 * g_1^0 * (a/r)^3
    # B_θ = 0
    # B_φ = 0
    expected_Br = 2.0 * g[2, 1] / (r^3)
    @test B[1] ≈ expected_Br rtol = 1.0e-6
    @test abs(B[2]) < 1.0e-6  # B_θ ≈ 0
    @test abs(B[3]) < 1.0e-10  # B_φ ≈ 0
end
