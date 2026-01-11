from shapely._geometry import Geometry
from shapely.geometry import LineString, MultiLineString, LinearRing, Polygon
from shapely.set_operations import union


alias CapStyle = Int32
alias JoinStyle = Int32

alias CAP_ROUND = 1
alias CAP_FLAT = 2
alias CAP_SQUARE = 3

alias JOIN_ROUND = 1
alias JOIN_BEVEL = 2
alias JOIN_MITRE = 3


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


fn _unit_tangent(ax: Float64, ay: Float64, bx: Float64, by: Float64) -> (Float64, Float64):
    var dx = bx - ax
    var dy = by - ay
    var len = sqrt_f64(dx * dx + dy * dy)
    if len == 0.0:
        return (0.0, 0.0)
    return (dx / len, dy / len)


fn _unit_normal_left(tx: Float64, ty: Float64) -> (Float64, Float64):
    return (-ty, tx)


fn _dot(ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
    return ax * bx + ay * by


fn _cross(ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
    return ax * by - ay * bx


fn _line_intersection(
    p1x: Float64,
    p1y: Float64,
    d1x: Float64,
    d1y: Float64,
    p2x: Float64,
    p2y: Float64,
    d2x: Float64,
    d2y: Float64,
) -> (Tuple[Float64, Float64], Bool):
    # Solve p1 + t*d1 = p2 + u*d2
    var denom = _cross(d1x, d1y, d2x, d2y)
    if denom == 0.0:
        return ((0.0, 0.0), False)
    var rx = p2x - p1x
    var ry = p2y - p1y
    var t = _cross(rx, ry, d2x, d2y) / denom
    return ((p1x + t * d1x, p1y + t * d1y), True)


fn _append_arc_8(
    mut out: List[Tuple[Float64, Float64]],
    cx: Float64,
    cy: Float64,
    r: Float64,
    sx: Float64,
    sy: Float64,
    ex: Float64,
    ey: Float64,
    ccw: Bool,
):
    # Approximate an arc using the 8-direction unit circle points.
    # We choose the closest direction indices to (s) and (e), then step.
    var s = sqrt_f64(0.5)
    var dirs = List[Tuple[Float64, Float64]]()
    dirs.append((1.0, 0.0))
    dirs.append((s, s))
    dirs.append((0.0, 1.0))
    dirs.append((-s, s))
    dirs.append((-1.0, 0.0))
    dirs.append((-s, -s))
    dirs.append((0.0, -1.0))
    dirs.append((s, -s))

    fn best_idx(vx: Float64, vy: Float64, dirs: List[Tuple[Float64, Float64]]) -> Int32:
        var best = -1.0e308
        var besti: Int32 = 0
        var i = 0
        while i < dirs.__len__():
            var d = _dot(vx, vy, dirs[i][0], dirs[i][1])
            if d > best:
                best = d
                besti = Int32(i)
            i += 1
        return besti

    var si = Int(best_idx(sx, sy, dirs))
    var ei = Int(best_idx(ex, ey, dirs))

    var i = si
    # include start direction
    out.append((cx + dirs[i][0] * r, cy + dirs[i][1] * r))
    if ccw:
        while i != ei:
            i = (i + 1) % 8
            out.append((cx + dirs[i][0] * r, cy + dirs[i][1] * r))
    else:
        while i != ei:
            i = (i + 7) % 8
            out.append((cx + dirs[i][0] * r, cy + dirs[i][1] * r))


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


fn buffer(
    geom: Geometry,
    _distance: Float64,
    _quad_segs: Int32,
    cap_style: CapStyle,
    join_style: JoinStyle,
    mitre_limit: Float64 = 5.0,
) -> Geometry:
    if _distance <= 0.0:
        return geom.copy()
    if geom.is_linestring():
        return buffer(geom.as_linestring(), _distance, _quad_segs, cap_style, join_style, mitre_limit)
    if geom.is_multilinestring():
        return buffer(geom.as_multilinestring(), _distance, _quad_segs, cap_style, join_style, mitre_limit)
    return geom.copy()


fn buffer(ls: LineString, distance: Float64, _quad_segs: Int32 = 8) -> Geometry:
    return buffer(ls, distance, _quad_segs, CAP_ROUND, JOIN_ROUND, 5.0)


fn buffer(
    ls: LineString,
    distance: Float64,
    _quad_segs: Int32,
    cap_style: CapStyle,
    join_style: JoinStyle,
    _mitre_limit: Float64 = 5.0,
) -> Geometry:
    if ls.coords.__len__() < 2:
        return Geometry(_empty_polygon())

    # Build a single polygon ring for the polyline buffer (no unions).
    var n = ls.coords.__len__()

    # Optionally extend endpoints for square caps
    var pts = ls.coords.copy()
    if cap_style == CAP_SQUARE:
        var t0 = _unit_tangent(pts[0][0], pts[0][1], pts[1][0], pts[1][1])
        pts[0] = (pts[0][0] - t0[0] * distance, pts[0][1] - t0[1] * distance)
        var t1 = _unit_tangent(pts[n - 2][0], pts[n - 2][1], pts[n - 1][0], pts[n - 1][1])
        pts[n - 1] = (pts[n - 1][0] + t1[0] * distance, pts[n - 1][1] + t1[1] * distance)

    # Precompute tangents and normals for each segment
    var txs = List[Float64]()
    var tys = List[Float64]()
    var nxs = List[Float64]()
    var nys = List[Float64]()
    var i = 0
    while i < n - 1:
        var t = _unit_tangent(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1])
        txs.append(t[0])
        tys.append(t[1])
        var nn = _unit_normal_left(t[0], t[1])
        nxs.append(nn[0])
        nys.append(nn[1])
        i += 1

    # Left and right offset polylines (joined at vertices)
    var left = List[Tuple[Float64, Float64]]()
    var right = List[Tuple[Float64, Float64]]()

    # Start vertex
    left.append((pts[0][0] + nxs[0] * distance, pts[0][1] + nys[0] * distance))
    right.append((pts[0][0] - nxs[0] * distance, pts[0][1] - nys[0] * distance))

    # Internal vertices
    var k = 1
    while k < n - 1:
        var px = pts[k][0]
        var py = pts[k][1]

        var n0x = nxs[k - 1]
        var n0y = nys[k - 1]
        var n1x = nxs[k]
        var n1y = nys[k]
        var t0x = txs[k - 1]
        var t0y = tys[k - 1]
        var t1x = txs[k]
        var t1y = tys[k]

        # Left side intersection
        var p0x = px + n0x * distance
        var p0y = py + n0y * distance
        var p1x = px + n1x * distance
        var p1y = py + n1y * distance
        var li = _line_intersection(p0x, p0y, t0x, t0y, p1x, p1y, t1x, t1y)
        var lpt = li[0]
        var lok = li[1]

        # Right side intersection (use -normals)
        var q0x = px - n0x * distance
        var q0y = py - n0y * distance
        var q1x = px - n1x * distance
        var q1y = py - n1y * distance
        var ri = _line_intersection(q0x, q0y, t0x, t0y, q1x, q1y, t1x, t1y)
        var rpt = ri[0]
        var rok = ri[1]

        if join_style == JOIN_ROUND:
            # Only add the round arc on the *outer* (convex) side of the turn.
            # On the inner (concave) side, connect with the offset-line intersection.
            var cr = _cross(t0x, t0y, t1x, t1y)
            if cr > 0.0:
                # Left turn: left side is outer.
                _append_arc_8(left, px, py, distance, n0x, n0y, n1x, n1y, True)
                if rok:
                    right.append(rpt)
                else:
                    right.append((q0x, q0y))
                    right.append((q1x, q1y))
            elif cr < 0.0:
                # Right turn: right side is outer.
                _append_arc_8(right, px, py, distance, -n0x, -n0y, -n1x, -n1y, False)
                if lok:
                    left.append(lpt)
                else:
                    left.append((p0x, p0y))
                    left.append((p1x, p1y))
            else:
                # Straight: just join with intersections if available.
                if lok:
                    left.append(lpt)
                else:
                    left.append((p0x, p0y))
                    left.append((p1x, p1y))
                if rok:
                    right.append(rpt)
                else:
                    right.append((q0x, q0y))
                    right.append((q1x, q1y))
        elif join_style == JOIN_MITRE and lok and rok:
            # Apply a simple mitre limit: fallback to bevel if too far
            var dlx = lpt[0] - px
            var dly = lpt[1] - py
            var drx = rpt[0] - px
            var dry = rpt[1] - py
            var ml = sqrt_f64(dlx * dlx + dly * dly)
            var mr = sqrt_f64(drx * drx + dry * dry)
            if ml <= _mitre_limit * distance:
                left.append(lpt)
            else:
                left.append((p0x, p0y))
                left.append((p1x, p1y))
            if mr <= _mitre_limit * distance:
                right.append(rpt)
            else:
                right.append((q0x, q0y))
                right.append((q1x, q1y))
        else:
            # bevel (or parallel fallback)
            left.append((p0x, p0y))
            left.append((p1x, p1y))
            right.append((q0x, q0y))
            right.append((q1x, q1y))

        k += 1

    # End vertex
    left.append((pts[n - 1][0] + nxs[n - 2] * distance, pts[n - 1][1] + nys[n - 2] * distance))
    right.append((pts[n - 1][0] - nxs[n - 2] * distance, pts[n - 1][1] - nys[n - 2] * distance))

    # Assemble ring: left forward + cap + right reversed + cap
    var ring = List[Tuple[Float64, Float64]]()
    for p in left:
        ring.append(p)

    # end cap
    if cap_style == CAP_ROUND:
        var nx = nxs[n - 2]
        var ny = nys[n - 2]
        # arc from +n to -n around end, ccw
        _append_arc_8(ring, pts[n - 1][0], pts[n - 1][1], distance, nx, ny, -nx, -ny, True)

    var rr = right.__len__() - 1
    while rr >= 0:
        ring.append(right[rr])
        rr -= 1

    # start cap
    if cap_style == CAP_ROUND:
        var nx0 = nxs[0]
        var ny0 = nys[0]
        # arc from -n to +n around start, ccw
        _append_arc_8(ring, pts[0][0], pts[0][1], distance, -nx0, -ny0, nx0, ny0, True)

    if ring.__len__() > 0:
        var first = ring[0]
        var last = ring[ring.__len__() - 1]
        if first[0] != last[0] or first[1] != last[1]:
            ring.append(first)

    return Geometry(Polygon(LinearRing(ring)))


fn buffer(mls: MultiLineString, distance: Float64, quad_segs: Int32 = 8) -> Geometry:
    var acc = Geometry(_empty_polygon())
    for ln in mls.lines:
        acc = union(acc, buffer(ln.copy(), distance, quad_segs))
    return acc.copy()


fn buffer(
    mls: MultiLineString,
    distance: Float64,
    quad_segs: Int32,
    cap_style: CapStyle,
    join_style: JoinStyle,
    mitre_limit: Float64 = 5.0,
) -> Geometry:
    # Avoid expensive dissolving/union: return a MultiPolygon of per-line buffers.
    var polys = List[Polygon]()
    for ln in mls.lines:
        var g = buffer(ln.copy(), distance, quad_segs, cap_style, join_style, mitre_limit)
        if g.is_polygon():
            polys.append(g.as_polygon())
        elif g.is_multipolygon():
            for p in g.as_multipolygon().polys:
                polys.append(p.copy())
    return Geometry(MultiPolygon(polys))


fn simplify(geom: Geometry, _tolerance: Float64, _preserve_topology: Bool = True) -> Geometry:
    return geom


fn convex_hull(geom: Geometry) -> Geometry:
    return geom
