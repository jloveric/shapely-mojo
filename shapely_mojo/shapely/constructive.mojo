from shapely._geometry import Geometry
from shapely.geometry import LineString, MultiLineString, LinearRing, Polygon
from shapely.set_operations import union


fn sqrt_f64(x: Float64) -> Float64:
    if x <= 0.0:
        return 0.0
    var r = x
    var i = 0
    while i < 12:
        r = 0.5 * (r + x / r)
        i += 1
    return r


fn _empty_polygon() -> Polygon:
    return Polygon(LinearRing(List[Tuple[Float64, Float64]]()))


fn _circle_polygon(cx: Float64, cy: Float64, r: Float64) -> Polygon:
    if r <= 0.0:
        return _empty_polygon()
    var s = sqrt_f64(0.5)
    var pts = List[Tuple[Float64, Float64]]()
    pts.append((cx + r, cy))
    pts.append((cx + r * s, cy + r * s))
    pts.append((cx, cy + r))
    pts.append((cx - r * s, cy + r * s))
    pts.append((cx - r, cy))
    pts.append((cx - r * s, cy - r * s))
    pts.append((cx, cy - r))
    pts.append((cx + r * s, cy - r * s))
    pts.append((cx + r, cy))
    return Polygon(LinearRing(pts))


fn _segment_tube(ax: Float64, ay: Float64, bx: Float64, by: Float64, r: Float64) -> Polygon:
    var dx = bx - ax
    var dy = by - ay
    var len = sqrt_f64(dx * dx + dy * dy)
    if len == 0.0 or r <= 0.0:
        return _empty_polygon()
    var nx = -dy / len
    var ny = dx / len
    var pts = List[Tuple[Float64, Float64]]()
    pts.append((ax + nx * r, ay + ny * r))
    pts.append((bx + nx * r, by + ny * r))
    pts.append((bx - nx * r, by - ny * r))
    pts.append((ax - nx * r, ay - ny * r))
    pts.append((ax + nx * r, ay + ny * r))
    return Polygon(LinearRing(pts))


fn buffer(geom: Geometry, _distance: Float64, _quad_segs: Int32 = 8) -> Geometry:
    if _distance <= 0.0:
        return geom.copy()
    if geom.is_linestring():
        return buffer(geom.as_linestring(), _distance, _quad_segs)
    if geom.is_multilinestring():
        return buffer(geom.as_multilinestring(), _distance, _quad_segs)
    return geom.copy()


fn buffer(ls: LineString, distance: Float64, _quad_segs: Int32 = 8) -> Geometry:
    if ls.coords.__len__() < 2:
        return Geometry(_empty_polygon())

    var acc = Geometry(_empty_polygon())

    var start = ls.coords[0]
    acc = union(acc, Geometry(_circle_polygon(start[0], start[1], distance)))
    var endp = ls.coords[ls.coords.__len__() - 1]
    acc = union(acc, Geometry(_circle_polygon(endp[0], endp[1], distance)))

    var i = 0
    while i < ls.coords.__len__() - 1:
        var a = ls.coords[i]
        var b = ls.coords[i + 1]
        var tube = _segment_tube(a[0], a[1], b[0], b[1], distance)
        acc = union(acc, Geometry(tube.copy()))
        i += 1
    return acc.copy()


fn buffer(mls: MultiLineString, distance: Float64, quad_segs: Int32 = 8) -> Geometry:
    var acc = Geometry(_empty_polygon())
    for ln in mls.lines:
        acc = union(acc, buffer(ln.copy(), distance, quad_segs))
    return acc.copy()


fn simplify(geom: Geometry, _tolerance: Float64, _preserve_topology: Bool = True) -> Geometry:
    return geom


fn convex_hull(geom: Geometry) -> Geometry:
    return geom
