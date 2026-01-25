using Test
using Dates
using StaticArrays
using GeoCotrans
using TsyganenkoModels
using PlanetaryMagneticFields

struct DummyModel <: PlanetaryMagneticFields.MagneticFieldModel
    field::SVector{3, Float64}
end

(m::DummyModel)(r, θ, φ) = m.field
PlanetaryMagneticFields.evaluate_field_spherical(m::DummyModel, r, θ, φ) = m.field

@testset "CompositeModel summation" begin
    base = PlanetaryMagneticFields.Earth
    m1 = MagneticModel(DummyModel(SVector(1.0, 0.0, 0.0)), base)
    m2 = MagneticModel(DummyModel(SVector(0.0, 2.0, 0.0)), base)
    comp = compose_models(m1, m2)

    @test comp(1.0, π / 2, 0.0) == SVector(1.0, 2.0, 0.0)
    cart = comp(1.0, π / 2, 0.0; out = :cartesian)
    @test all(isfinite, cart)
end

@testset "Earth IGRF + Tsyganenko" begin
    drivers = TsyganenkoDrivers(Date(2020, 1, 1); pdyn = 2.0, dst = 0.0, byimf = 0.0, bzimf = 0.0)
    igrf = "igrf2020"
    composite = load_model(:earth, "t96"; drivers = drivers, igrf_model = igrf, max_degree = 4)
    r, θ, φ = 5.0, π / 2, 0.0

    internal = load_model(:earth, igrf; max_degree = 4)
    B_internal = internal(r, θ, φ)

    pos_geo = GeoCotrans.GEO(PlanetaryMagneticFields.spherical_to_cartesian(r, θ, φ)...)
    pos_gsm = GeoCotrans.geo2gsm(pos_geo, drivers.time)
    B_gsm = TsyganenkoModels.t96(pos_gsm, drivers.time, drivers.pdyn, drivers.dst, drivers.byimf, drivers.bzimf)
    B_geo = GeoCotrans.gsm2geo(B_gsm, drivers.time)
    B_external = PlanetaryMagneticFields.cartesian_field_to_spherical(B_geo..., θ, φ)

    @test composite(r, θ, φ) ≈ B_internal + B_external
end
