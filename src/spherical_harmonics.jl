"""
Spherical harmonic field evaluation for PlanetaryMagneticFields.jl

This module implements Schmidt semi-normalized associated Legendre polynomials
and spherical harmonic magnetic field evaluation.
"""

using StaticArrays
using LinearAlgebra

"""
    schmidt_legendre!(P::Matrix{Float64}, dP::Matrix{Float64}, cosθ, sinθ, max_degree::Int)

Compute Schmidt semi-normalized associated Legendre polynomials and their derivatives.

This function computes P_n^m(cos θ) and dP_n^m/dθ for n = 0 to max_degree and m = 0 to n,
using the Schmidt semi-normalization convention standard in geomagnetism.

# Arguments
- `P::Matrix{Float64}`: Output matrix for P_n^m values (size: (max_degree+1) × (max_degree+1))
- `dP::Matrix{Float64}`: Output matrix for dP_n^m/dθ values (same size as P)
- `cosθ`: cos(colatitude)
- `sinθ`: sin(colatitude)
- `max_degree::Int`: Maximum degree to compute

# Schmidt Semi-Normalization
The Schmidt semi-normalized Legendre polynomials satisfy:
```math
\\int_0^π [P_n^m(\\cos θ)]² \\sin θ dθ = \\frac{2(n+m)!}{(2n+1)(n-m)!}
```
"""
function schmidt_legendre!(
        P::Matrix{Float64}, dP::Matrix{Float64},
        cosθ, sinθ, max_degree::Int
    )
    # Initialize
    fill!(P, 0.0)
    fill!(dP, 0.0)

    # P_0^0 = 1 (degree 0, order 0)
    P[1, 1] = 1.0
    dP[1, 1] = 0.0

    # Handle the case when we're at the pole
    if abs(sinθ) < 1.0e-10
        # At the poles, only m=0 terms are non-zero
        for n in 1:max_degree
            # Use recursion for m=0 only
            if n == 1
                P[2, 1] = cosθ
                dP[2, 1] = -sinθ
            else
                k = Float64(n)
                P[n + 1, 1] = ((2k - 1) / k) * cosθ * P[n, 1] - ((k - 1) / k) * P[n - 1, 1]
                dP[n + 1, 1] = ((2k - 1) / k) * (cosθ * dP[n, 1] - sinθ * P[n, 1]) -
                    ((k - 1) / k) * dP[n - 1, 1]
            end
        end
        return
    end

    # General case: not at the poles
    # Using Schmidt quasi-normalized recursion formulas (geomagnetic convention)

    # P_1^1 = sinθ is an initial condition for Schmidt quasi-normalization
    P[2, 2] = sinθ
    dP[2, 2] = cosθ

    # Diagonal terms for n >= 2: P_n^n = sqrt((2n-1)/(2n)) * sinθ * P_{n-1}^{n-1}
    for m in 2:max_degree
        fm = Float64(m)
        scale = sqrt((2fm - 1) / (2fm))
        P[m + 1, m + 1] = sinθ * P[m, m] * scale
        dP[m + 1, m + 1] = (sinθ * dP[m, m] + cosθ * P[m, m]) * scale
    end

    # Off-diagonal terms: P_n^m for n > m
    # General recursion: P_n^m = a_nm * cosθ * P_{n-1}^m - b_nm * P_{n-2}^m
    # where a_nm = (2n-1) / sqrt(n² - m²)
    #       b_nm = sqrt((n-1)² - m²) / sqrt(n² - m²)
    for m in 0:max_degree
        fm = Float64(m)
        for n in (m + 1):max_degree
            fn = Float64(n)

            # Compute coefficients
            n2_m2 = fn * fn - fm * fm  # n² - m²
            a_nm = (2fn - 1) / sqrt(n2_m2)

            if n == m + 1
                # First off-diagonal: b_nm term is zero (P_{n-2}^m doesn't exist)
                P[n + 1, m + 1] = a_nm * cosθ * P[n, m + 1]
                dP[n + 1, m + 1] = a_nm * (cosθ * dP[n, m + 1] - sinθ * P[n, m + 1])
            else
                # General case with both terms
                nm1_2_m2 = (fn - 1) * (fn - 1) - fm * fm  # (n-1)² - m²
                b_nm = sqrt(nm1_2_m2 / n2_m2)

                P[n + 1, m + 1] = a_nm * cosθ * P[n, m + 1] - b_nm * P[n - 1, m + 1]
                dP[n + 1, m + 1] = a_nm * (cosθ * dP[n, m + 1] - sinθ * P[n, m + 1]) -
                    b_nm * dP[n - 1, m + 1]
            end
        end
    end
    return
end

"""
    evaluate_field_spherical(model::SphericalHarmonicModel, r, θ, φ)

Evaluate the magnetic field at a point in spherical coordinates.

# Arguments
- `model::SphericalHarmonicModel`: The magnetic field model
- `r`: Radial distance in planetary radii (dimensionless)
- `θ`: Colatitude in radians [0, π]
- `φ`: East longitude in radians [0, 2π]

# Returns
- `SVector{3,Float64}`: Magnetic field vector [B_r, B_θ, B_φ] in nanoTesla

# Coordinate System
- Spherical coordinates (r, θ, φ) where:
  - r: radial distance from planet center
  - θ: colatitude (0 at north pole, π at south pole)
  - φ: east longitude (0 at prime meridian, increases eastward)

# Mathematical Formulation
The magnetic field is derived from the scalar potential:
```math
B_r = -\\frac{∂V}{∂r}, \\quad
B_θ = -\\frac{1}{r}\\frac{∂V}{∂θ}, \\quad
B_φ = -\\frac{1}{r\\sin θ}\\frac{∂V}{∂φ}
```
"""
function evaluate_field_spherical(model::SphericalHarmonicModel, r, θ, φ)
    r > 0 || error("Radius must be positive")
    0 <= θ <= π || error("Colatitude θ must be in [0, π]")

    coeffs = model.coeffs
    max_degree = coeffs.degree

    # Precompute trigonometric functions
    cosθ = cos(θ)
    sinθ = sin(θ)

    # Allocate arrays for Legendre polynomials
    P = zeros(Float64, max_degree + 1, max_degree + 1)
    dP = zeros(Float64, max_degree + 1, max_degree + 1)

    # Compute Schmidt-normalized associated Legendre polynomials
    schmidt_legendre!(P, dP, cosθ, sinθ, max_degree)

    # Initialize field components
    Br = 0.0
    Bθ = 0.0
    Bφ = 0.0

    # Precompute powers of (1/r)
    # r is in planetary radii (dimensionless), so (1/r)^k gives the radial dependence
    # For internal field: B components use (1/r)^(n+2)
    ratio = 1.0 / r
    ratio_powers = zeros(Float64, max_degree + 3)
    ratio_powers[1] = ratio
    for k in 2:(max_degree + 3)
        ratio_powers[k] = ratio_powers[k - 1] * ratio
    end

    # Sum over degrees and orders
    for n in 1:max_degree
        fn = Float64(n)
        ratio_pow = ratio_powers[n + 2]  # (1/r)^(n+2) for magnetic field

        for m in 0:min(n, coeffs.order)
            # Get coefficients
            g = coeffs.g[n + 1, m + 1]
            h = coeffs.h[n + 1, m + 1]

            # Compute trigonometric terms
            sin_mφ, cos_mφ = sincos(m * φ)

            # Spherical harmonic term
            Y = g * cos_mφ + h * sin_mφ

            # Derivatives
            dY_dφ = m * (-g * sin_mφ + h * cos_mφ)

            # Contribution to field components
            # B_r = -∂V/∂r: factor of (n+1) from derivative of (a/r)^(n+1)
            Br += (fn + 1) * ratio_pow * Y * P[n + 1, m + 1]

            # B_θ = -(1/r)∂V/∂θ
            Bθ -= ratio_pow * Y * dP[n + 1, m + 1]

            # B_φ = -(1/(r sin θ))∂V/∂φ
            if abs(sinθ) > 1.0e-10
                Bφ -= ratio_pow * dY_dφ * P[n + 1, m + 1] / sinθ
            end
        end
    end

    return SVector{3, Float64}(Br, Bθ, Bφ)
end


"""
    evaluate_field_cartesian(model::SphericalHarmonicModel, x, y, z)

Evaluate the magnetic field at a point in Cartesian coordinates.

# Arguments
- `model::SphericalHarmonicModel`: The magnetic field model
- `x, y, z`: Cartesian coordinates in planetary radii

# Returns
- `SVector{3,Float64}`: Magnetic field vector [Bx, By, Bz] in nanoTesla
"""
function evaluate_field_cartesian(
        model::SphericalHarmonicModel,
        x, y, z
    )
    # Convert to spherical coordinates
    (r, θ, φ) = cartesian_to_spherical(x, y, z)

    # Evaluate field in spherical coordinates
    B_sph = evaluate_field_spherical(model, r, θ, φ)

    # Convert field to Cartesian coordinates
    B_cart = spherical_field_to_cartesian(B_sph[1], B_sph[2], B_sph[3], θ, φ)

    return B_cart
end
