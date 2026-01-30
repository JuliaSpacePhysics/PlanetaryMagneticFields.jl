using LinearAlgebra
using AxisKeys

_field_func(idx::Int) = idx <= 3 ? x -> getindex(x, idx) : norm
_field_func(::Nothing) = identity
_field_func(f) = f

function fieldmap(model, r, nlat, nlon; idx = identity)
    func = _field_func(idx)
    # Create latitude and longitude grids
    # Latitude: -90 to 90 degrees (convert to colatitude for evaluation)
    # Longitude: -180 to 180 degrees
    lats = range(-90, 90, length = nlat)
    lons = range(-180, 180, length = nlon)

    # Preallocate field array (longitude × latitude for GeoMakie surface!)
    field = zeros(Float64, nlon, nlat)

    θs = @. deg2rad(90 - lats)'
    φs = @. deg2rad(mod(lons, 360))
    field = func.(model.(r, θs, φs; in = :spherical, out = :spherical))

    return KeyedArray(field; lon = lons, lat = lats)
end

function fieldmap(model, r = 1.0; nlat = 180, nlon = 360, kw...)
    return fieldmap(model, r, nlat, nlon; kw...)
end

function evalmodel(m::SphericalHarmonicModel{<:GaussCoefficients}, r, θ, φ, _)
    return evalsph(m.coeffs, r, θ, φ, m.degree, m.order)
end

# Time-varying evaluation
function evalmodel(m::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients}, r, θ, φ, t)
    @assert r > 0
    return evalsph(m.coeffs(t), r, θ, φ, m.degree, m.order)
end
