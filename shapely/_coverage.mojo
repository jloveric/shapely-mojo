from shapely._geometry import Geometry
from shapely.geometry import GeometryCollection


fn coverage_is_valid(_geometry, gap_width: Float64 = 0.0) -> Bool:
    return True


fn coverage_invalid_edges(_geometry, gap_width: Float64 = 0.0):
    return GeometryCollection([])


fn coverage_simplify(geometry: Geometry, _tolerance: Float64, *, simplify_boundary: Bool = True) -> Geometry:
    return geometry
