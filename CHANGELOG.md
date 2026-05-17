# Changelog

## [Unreleased]

## [0.2.0] - 2026-05-16

### Changed

- **Breaking**: `fieldmap` now returns a lazy `LazyFieldMap <: AbstractMatrix{SVector{3,Float64}}` instead of an `AxisKeys.KeyedArray` and no longer accepts the `idx` kwarg. Each element is `[Br, Bθ, Bφ]`, computed on access — projecting via `getindex.(field, 1)` or `norm.(field)` materializes only the requested component. Recover axes via `field.lons` and `field.lats` (degrees).
- Removed `AxisKeys` and `StaticArrays` as hard dependencies.


[unreleased]: https://github.com/JuliaSpacePhysics/PlanetaryMagneticFields.jl/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/JuliaSpacePhysics/PlanetaryMagneticFields.jl/releases/tag/v0.2.0
