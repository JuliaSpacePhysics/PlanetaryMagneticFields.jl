using Dates: AbstractTime, Date
using GeoCotrans
using TsyganenkoModels

const TSYGANENKO_MODELS = Set([:t89, :t96, :t01, :ts04])

"""
    TsyganenkoDrivers(time = Date(2025, 1, 1); pdyn = 2.0, dst = 0.0, byimf = 0.0, bzimf = 0.0)

Configuration for Tsyganenko external field models.
"""
struct TsyganenkoDrivers
    time::AbstractTime
    pdyn::Float64
    dst::Float64
    byimf::Float64
    bzimf::Float64
end

TsyganenkoDrivers(
    time::AbstractTime = Date(2025, 1, 1);
    pdyn = 2.0,
    dst = 0.0,
    byimf = 0.0,
    bzimf = 0.0,
) = TsyganenkoDrivers(time, pdyn, dst, byimf, bzimf)

"""
    TsyganenkoModel(variant, drivers)

Wrap a Tsyganenko external model (T89, T96, T01, TS04) so it can be composed
with other `MagneticModel`s.
"""
struct TsyganenkoModel <: ExternalFieldModel
    variant::Symbol
    drivers::TsyganenkoDrivers

    function TsyganenkoModel(variant::Symbol, drivers::TsyganenkoDrivers = TsyganenkoDrivers())
        variant in TSYGANENKO_MODELS || error("Unsupported Tsyganenko model: $variant. Supported: $(collect(TSYGANENKO_MODELS))")
        return new(variant, drivers)
    end
end

function _tsyg_fn(variant::Symbol)
    return getfield(TsyganenkoModels, variant)
end

function (m::TsyganenkoModel)(r, θ, φ)
    pos_geo = GeoCotrans.GEO(spherical_to_cartesian(r, θ, φ)...)
    pos_gsm = GeoCotrans.geo2gsm(pos_geo, m.drivers.time)
    fn = _tsyg_fn(m.variant)
    B_gsm = fn(pos_gsm[1], pos_gsm[2], pos_gsm[3], m.drivers.time, m.drivers.pdyn, m.drivers.dst, m.drivers.byimf, m.drivers.bzimf)
    B_vec = B_gsm isa GeoCotrans.CoordinateVector ? B_gsm : GeoCotrans.GSM(B_gsm)
    B_geo = GeoCotrans.gsm2geo(B_vec, m.drivers.time)
    return cartesian_field_to_spherical(B_geo..., θ, φ)
end

tsyg_model(planet::Planet, variant::Symbol, drivers::TsyganenkoDrivers) =
    MagneticModel(TsyganenkoModel(variant, drivers), planet)
