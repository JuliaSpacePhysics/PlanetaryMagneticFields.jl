module PlanetaryMagneticFieldsUnitfulExt

using PlanetaryMagneticFields
using PlanetaryMagneticFields: SphericalHarmonicModel
using Unitful

function (m::SphericalHarmonicModel)(r::Unitful.Length, θ, φ; kw...)
    # Convert position to planetary radii
    r_radii = ustrip(u"km", r) / m.obj.radius
    return m(r_radii, θ, φ; kw...)
end

# Cartesian coordinates with units
function (m::SphericalHarmonicModel)(x::Unitful.Length, y::Unitful.Length, z::Unitful.Length; kw...)
    # Convert positions to planetary radii
    x_radii = ustrip(u"km", x) / m.obj.radius
    y_radii = ustrip(u"km", y) / m.obj.radius
    z_radii = ustrip(u"km", z) / m.obj.radius
    return m(x_radii, y_radii, z_radii; in = :cartesian, kw...)
end

end # module
