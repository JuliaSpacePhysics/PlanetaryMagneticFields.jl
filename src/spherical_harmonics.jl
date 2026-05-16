"""
    LazyFieldMap{T,M,R,Θ,Φ} <: AbstractMatrix{T}

Lazy `(nlon, nlat)` grid of magnetic field vectors at radial distance `r`.
Each `field[i, j]` calls `model(r, θs[j], φs[i]; in = :spherical, out = :spherical)`
on access, so broadcasting a projection (`getindex.(field, 1)`, `norm.(field)`)
materializes only the requested scalar array — no intermediate `Matrix{SVector{3}}`.

# Properties
- `model`, `r`, `θs`, `φs`: stored fields (radians).
- `lons`, `lats`: derived axes in degrees, suitable for plotting:
  `lons = range(-180, 180, length = nlon)`, `lats = range(-90, 90, length = nlat)`.
"""
struct LazyFieldMap{T, M, R, Θ, Φ} <: AbstractMatrix{T}
    model::M
    r::R
    θs::Θ
    φs::Φ
end

function Base.getproperty(m::LazyFieldMap, name::Symbol)
    name === :lons && return range(-180, 180, length = length(getfield(m, :φs)))
    name === :lats && return range(-90, 90, length = length(getfield(m, :θs)))
    return getfield(m, name)
end

Base.propertynames(::LazyFieldMap) = (:model, :r, :θs, :φs, :lons, :lats)

function LazyFieldMap(model, r, θs::AbstractVector, φs::AbstractVector)
    T = Base.promote_op(model, typeof(r), eltype(θs), eltype(φs))
    return LazyFieldMap{T, typeof(model), typeof(r), typeof(θs), typeof(φs)}(model, r, θs, φs)
end

Base.size(m::LazyFieldMap) = (length(m.φs), length(m.θs))
Base.IndexStyle(::Type{<:LazyFieldMap}) = IndexCartesian()
@inline Base.getindex(m::LazyFieldMap, i::Int, j::Int) =
    @inbounds m.model(m.r, m.θs[j], m.φs[i]; in = :spherical, out = :spherical)

function fieldmap(model, r, nlat, nlon)
    # Latitude spans -90 to 90 (south to north); θ = colatitude = π/2 - lat.
    # Longitude spans -180 to 180; sin/cos are 2π-periodic so no wrap needed.
    θs = range(π, 0, length = nlat)
    φs = deg2rad.(range(-180, 180, length = nlon))
    return LazyFieldMap(model, r, θs, φs)
end

fieldmap(model, r = 1.0; nlat = 180, nlon = 360) = fieldmap(model, r, nlat, nlon)

function evalmodel(m::SphericalHarmonicModel{<:GaussCoefficients}, r, θ, φ, _)
    return evalsph(m.coeffs, r, θ, φ, m.degree, m.order)
end

# Time-varying evaluation
function evalmodel(m::SphericalHarmonicModel{<:TimeVaryingGaussCoefficients}, r, θ, φ, t)
    @assert r > 0
    return evalsph(m.coeffs(t), r, θ, φ, m.degree, m.order)
end
