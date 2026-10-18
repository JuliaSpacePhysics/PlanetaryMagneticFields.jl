struct Planet
    name::Symbol
    radius::Float64  # Mean radius in km
end

struct PlanetFrame <: AbstractReferenceFrame end

SphericalHarmonicModel(name, coeffs, p::Planet; kw...) =
    SphericalHarmonicModel(name, coeffs; radius = p.radius, frame = PlanetFrame(), kw...)

degree(m::SphericalHarmonicModel) = m.degree
order(m::SphericalHarmonicModel) = m.order
