from shapely._geometry import Geometry


fn buffer(geom: Geometry, _distance: Float64, _quad_segs: Int32 = 8) -> Geometry:
    return geom


fn simplify(geom: Geometry, _tolerance: Float64, _preserve_topology: Bool = True) -> Geometry:
    return geom


fn convex_hull(geom: Geometry) -> Geometry:
    return geom
