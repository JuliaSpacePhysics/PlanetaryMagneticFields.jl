# Changelog

## [Unreleased]

### Breaking

- `fieldmap` now returns a lazy `LazyFieldMap <: AbstractMatrix{SVector{3,Float64}}`
  instead of an `AxisKeys.KeyedArray` and no longer accepts the `idx` kwarg.
  Each element is `[Bᵣ, Bθ, Bφ]`, computed on access — projecting via
  `getindex.(field, 1)` or `norm.(field)` materializes only the requested component.
  Recover axes via `field.lons` and `field.lats` (degrees).
- Removed `AxisKeys` and `StaticArrays` as hard dependencies.

### Changed