"""
    cartesian_to_spherical(x, y, z)

Convert Cartesian coordinates to spherical coordinates.

# Arguments
- `x, y, z`: Cartesian coordinates (same units, e.g., km or planetary radii)

# Returns
- `(r, θ, φ)`: Spherical coordinates
  - `r`: radial distance (same units as input)
  - `θ`: colatitude in radians [0, π]
  - `φ`: east longitude in radians [0, 2π]

# Coordinate Convention
- x-axis: points to (θ=π/2, φ=0)
- y-axis: points to (θ=π/2, φ=π/2)
- z-axis: points to north pole (θ=0)
"""
function cartesian_to_spherical(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)

    if r < 1.0e-10
        # At origin
        return (0.0, 0.0, 0.0)
    end

    θ = acos(clamp(z / r, -1.0, 1.0))  # Colatitude [0, π]
    φ = atan(y, x)  # Longitude [-π, π]

    # Normalize φ to [0, 2π]
    if φ < 0
        φ += 2π
    end

    return (r, θ, φ)
end

"""
    spherical_to_cartesian(r, θ, φ)

Convert spherical coordinates to Cartesian coordinates.

# Arguments
- `r`: radial distance
- `θ`: colatitude in radians [0, π]
- `φ`: east longitude in radians [0, 2π]

# Returns
- `(x, y, z)`: Cartesian coordinates (same units as r)
"""
function spherical_to_cartesian(r, θ, φ)
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    cosφ = cos(φ)

    x = r * sinθ * cosφ
    y = r * sinθ * sinφ
    z = r * cosθ

    return (x, y, z)
end

"""
    spherical_field_to_cartesian(Br, Bθ, Bφ,
                                   θ, φ)

Convert magnetic field vector from spherical to Cartesian components.

# Arguments
- `Br, Bθ, Bφ`: Field components in spherical coordinates
- `θ`: colatitude in radians [0, π]
- `φ`: east longitude in radians [0, 2π]

# Returns
- `SVector{3,Float64}`: Field vector [Bx, By, Bz] in Cartesian coordinates

# Transformation
```math
\\begin{pmatrix} B_x \\\\ B_y \\\\ B_z \\end{pmatrix} =
\\begin{pmatrix}
  \\sin θ \\cos φ & \\cos θ \\cos φ & -\\sin φ \\\\
  \\sin θ \\sin φ & \\cos θ \\sin φ &  \\cos φ \\\\
  \\cos θ         & -\\sin θ          &  0
\\end{pmatrix}
\\begin{pmatrix} B_r \\\\ B_θ \\\\ B_φ \\end{pmatrix}
```
"""
function spherical_field_to_cartesian(
        Br, Bθ, Bφ,
        θ, φ
    )
    sinθ = sin(θ)
    cosθ = cos(θ)
    sinφ = sin(φ)
    cosφ = cos(φ)

    Bx = Br * sinθ * cosφ + Bθ * cosθ * cosφ - Bφ * sinφ
    By = Br * sinθ * sinφ + Bθ * cosθ * sinφ + Bφ * cosφ
    Bz = Br * cosθ - Bθ * sinθ

    return SVector{3, Float64}(Bx, By, Bz)
end
