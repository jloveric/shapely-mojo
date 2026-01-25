from shapely._geometry import Geometry
from shapely.geometry import (
    Point,
    LineString,
    MultiLineString,
    LinearRing,
    Polygon,
    GeometryCollection,
    MultiPoint,
    MultiPolygon,
)
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


fn _collect_coords(geom: Geometry, mut out: List[Tuple[Float64, Float64]]):
    if geom.is_point():
        var p = geom.as_point()
        out.append((p.x, p.y))
        return
    if geom.is_linestring():
        var ls = geom.as_linestring()
        for c in ls.coords:
            out.append(c)
        return
    if geom.is_polygon():
        var p = geom.as_polygon()
        for c in p.shell.coords:
            out.append(c)
        for h in p.holes:
            for c in h.coords:
                out.append(c)
        return
    if geom.is_multipoint():
        var mp = geom.as_multipoint()
        for p in mp.points:
            out.append((p.x, p.y))
        return
    if geom.is_multilinestring():
        var mls = geom.as_multilinestring()
        for ln in mls.lines:
            for c in ln.coords:
                out.append(c)
        return
    if geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        for poly in mp.polys:
            for c in poly.shell.coords:
                out.append(c)
            for h in poly.holes:
                for c in h.coords:
                    out.append(c)
        return
    if geom.is_geometrycollection():
        var gc = geom.as_geometrycollection()
        for g in gc.geoms:
            _collect_coords(g.copy(), out)
        return


fn _sort_points_lex(mut pts: List[Tuple[Float64, Float64]]):
    # insertion sort
    var i = 1
    while i < pts.__len__():
        var v = pts[i]
        var j = i - 1
        while j >= 0:
            var pj = pts[j]
            if pj[0] < v[0] or (pj[0] == v[0] and pj[1] <= v[1]):
                break
            pts[j + 1] = pj
            j -= 1
        pts[j + 1] = v
        i += 1


fn _unique_sorted_points(pts: List[Tuple[Float64, Float64]]) -> List[Tuple[Float64, Float64]]:
    if pts.__len__() == 0:
        return pts.copy()
    var out = List[Tuple[Float64, Float64]]()
    out.append(pts[0])
    var i = 1
    while i < pts.__len__():
        var p = pts[i]
        var last = out[out.__len__() - 1]
        if p[0] != last[0] or p[1] != last[1]:
            out.append(p)
        i += 1
    return out.copy()


fn _hull_ring(points: List[Tuple[Float64, Float64]]) -> List[Tuple[Float64, Float64]]:
    if points.__len__() <= 1:
        return points.copy()
    var pts = points.copy()
    _sort_points_lex(pts)
    pts = _unique_sorted_points(pts)
    if pts.__len__() <= 1:
        return pts.copy()

    var lower = List[Tuple[Float64, Float64]]()
    for p in pts:
        while lower.__len__() >= 2:
            var b = lower[lower.__len__() - 1]
            var a = lower[lower.__len__() - 2]
            var cr = _cross(b[0] - a[0], b[1] - a[1], p[0] - b[0], p[1] - b[1])
            if cr > 0.0:
                break
            lower = lower[: lower.__len__() - 1]
        lower.append(p)

    var upper = List[Tuple[Float64, Float64]]()
    var i = pts.__len__() - 1
    while i >= 0:
        var p = pts[i]
        while upper.__len__() >= 2:
            var b = upper[upper.__len__() - 1]
            var a = upper[upper.__len__() - 2]
            var cr = _cross(b[0] - a[0], b[1] - a[1], p[0] - b[0], p[1] - b[1])
            if cr > 0.0:
                break
            upper = upper[: upper.__len__() - 1]
        upper.append(p)
        i -= 1

    # concatenate without duplicating endpoints
    var ring = List[Tuple[Float64, Float64]]()
    var j = 0
    while j < lower.__len__():
        ring.append(lower[j])
        j += 1
    j = 1
    while j + 1 < upper.__len__():
        ring.append(upper[j])
        j += 1
    return ring.copy()


fn _point_seg_dist2(px: Float64, py: Float64, ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
    var vx = bx - ax
    var vy = by - ay
    var wx = px - ax
    var wy = py - ay
    var c1 = vx * wx + vy * wy
    if c1 <= 0.0:
        var dx = px - ax
        var dy = py - ay
        return dx * dx + dy * dy
    var c2 = vx * vx + vy * vy
    if c2 <= c1:
        var dx2 = px - bx
        var dy2 = py - by
        return dx2 * dx2 + dy2 * dy2
    var t = c1 / c2
    var projx = ax + t * vx
    var projy = ay + t * vy
    var dx3 = px - projx
    var dy3 = py - projy
    return dx3 * dx3 + dy3 * dy3


fn _simplify_linestring(ls: LineString, tol: Float64) -> LineString:
    if ls.coords.__len__() <= 2:
        return ls.copy()
    if tol <= 0.0:
        return ls.copy()
    var tol2 = tol * tol
    var n = ls.coords.__len__()
    var keep = List[Bool]()
    var i = 0
    while i < n:
        keep.append(False)
        i += 1
    keep[0] = True
    keep[n - 1] = True

    var stack_a = List[Int]()
    var stack_b = List[Int]()
    stack_a.append(0)
    stack_b.append(n - 1)
    while stack_a.__len__() > 0:
        var a = stack_a[stack_a.__len__() - 1]
        var b = stack_b[stack_b.__len__() - 1]
        stack_a = stack_a[: stack_a.__len__() - 1]
        stack_b = stack_b[: stack_b.__len__() - 1]

        var ax = ls.coords[a][0]
        var ay = ls.coords[a][1]
        var bx = ls.coords[b][0]
        var by = ls.coords[b][1]
        var maxd = -1.0
        var maxi = -1
        var j = a + 1
        while j < b:
            var px = ls.coords[j][0]
            var py = ls.coords[j][1]
            var d2 = _point_seg_dist2(px, py, ax, ay, bx, by)
            if d2 > maxd:
                maxd = d2
                maxi = j
            j += 1
        if maxd > tol2 and maxi != -1:
            keep[maxi] = True
            stack_a.append(a)
            stack_b.append(maxi)
            stack_a.append(maxi)
            stack_b.append(b)

    var out = List[Tuple[Float64, Float64]]()
    i = 0
    while i < n:
        if keep[i]:
            out.append(ls.coords[i])
        i += 1
    return LineString(out)


fn _simplify_ring_coords(coords: List[Tuple[Float64, Float64]], tol: Float64) -> List[Tuple[Float64, Float64]]:
    if coords.__len__() < 4:
        return coords.copy()
    var open = coords.copy()
    var first = open[0]
    var last = open[open.__len__() - 1]
    if first[0] == last[0] and first[1] == last[1]:
        open = open[: open.__len__() - 1]
    if open.__len__() < 3:
        return List[Tuple[Float64, Float64]]()
    var simplified = _simplify_linestring(LineString(open), tol).coords.copy()
    if simplified.__len__() < 3:
        return List[Tuple[Float64, Float64]]()
    var out = simplified.copy()
    var f = out[0]
    var l = out[out.__len__() - 1]
    if f[0] != l[0] or f[1] != l[1]:
        out.append(f)
    return out.copy()


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


fn _line_intersection_tu(
    p1x: Float64,
    p1y: Float64,
    d1x: Float64,
    d1y: Float64,
    p2x: Float64,
    p2y: Float64,
    d2x: Float64,
    d2y: Float64,
) -> (Tuple[Float64, Float64], Float64, Float64, Bool):
    # Solve p1 + t*d1 = p2 + u*d2
    var denom = _cross(d1x, d1y, d2x, d2y)
    if denom == 0.0:
        return ((0.0, 0.0), 0.0, 0.0, False)
    var rx = p2x - p1x
    var ry = p2y - p1y
    var t = _cross(rx, ry, d2x, d2y) / denom
    var u = _cross(rx, ry, d1x, d1y) / denom
    return ((p1x + t * d1x, p1y + t * d1y), t, u, True)


fn _append_arc(
    mut out: List[Tuple[Float64, Float64]],
    cx: Float64,
    cy: Float64,
    r: Float64,
    quad_segs: Int32,
    sx: Float64,
    sy: Float64,
    ex: Float64,
    ey: Float64,
    ccw: Bool,
    force_short: Bool,
):
    # Approximate an arc using a direction table derived from the 8-direction unit
    # circle, subdividing each 45-degree octant using normalized interpolation.
    # quad_segs is the number of segments per 90-degree quarter circle.
    var s = sqrt_f64(0.5)
    var base = List[Tuple[Float64, Float64]]()
    base.append((1.0, 0.0))
    base.append((s, s))
    base.append((0.0, 1.0))
    base.append((-s, s))
    base.append((-1.0, 0.0))
    base.append((-s, -s))
    base.append((0.0, -1.0))
    base.append((s, -s))

    var segs: Int = Int(quad_segs)
    if segs < 1:
        segs = 1
    # Number of subdivisions per 45-degree octant.
    var per_oct: Int = Int(segs / 2)
    if per_oct < 1:
        per_oct = 1

    var dirs = List[Tuple[Float64, Float64]]()
    var bi = 0
    while bi < 8:
        var a = base[bi]
        var b = base[(bi + 1) % 8]
        var t: Int = 0
        while t <= per_oct:
            if bi != 0 and t == 0:
                t += 1
                continue
            var tt = Float64(t) / Float64(per_oct)
            var vx = a[0] * (1.0 - tt) + b[0] * tt
            var vy = a[1] * (1.0 - tt) + b[1] * tt
            var vl = sqrt_f64(vx * vx + vy * vy)
            if vl != 0.0:
                vx /= vl
                vy /= vl
            dirs.append((vx, vy))
            t += 1
        bi += 1

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

    var n = dirs.__len__()
    var i = si
    # Decide sweep. For joins we want the caller's exterior sweep direction.
    # For caps we want the short arc (and a stable semicircle for opposite vectors).
    var do_ccw = ccw
    var half: Int = Int(n / 2)

    # Detect opposite vectors (used for round caps): force exactly a semicircle.
    var sdot = sx * ex + sy * ey
    if force_short and sdot < -0.90:
        ei = (si + half) % n
        do_ccw = ccw
    elif force_short:
        var steps_ccw = (ei - si) % n
        if steps_ccw < 0:
            steps_ccw += n
        var steps_cw = (si - ei) % n
        if steps_cw < 0:
            steps_cw += n
        if ccw and steps_ccw > half:
            do_ccw = False
        if (not ccw) and steps_cw > half:
            do_ccw = True
    # include start direction
    out.append((cx + dirs[i][0] * r, cy + dirs[i][1] * r))
    if do_ccw:
        while i != ei:
            i = (i + 1) % n
            out.append((cx + dirs[i][0] * r, cy + dirs[i][1] * r))
    else:
        while i != ei:
            i = (i + n - 1) % n
            out.append((cx + dirs[i][0] * r, cy + dirs[i][1] * r))


fn _append_arc_join(
    mut out: List[Tuple[Float64, Float64]],
    cx: Float64,
    cy: Float64,
    r: Float64,
    quad_segs: Int32,
    sx: Float64,
    sy: Float64,
    ex: Float64,
    ey: Float64,
):
    # Build a short arc between start/end unit vectors by iterative subdivision.
    # This avoids direction-table quantization collapsing joins to a mitre.
    var segs: Int = Int(quad_segs) * 4
    if segs < 1:
        segs = 1
    var iters: Int = 0
    var pieces: Int = 1
    while pieces < segs:
        pieces *= 2
        iters += 1

    var dirs = List[Tuple[Float64, Float64]]()
    dirs.append((sx, sy))
    dirs.append((ex, ey))

    var k = 0
    while k < iters:
        var nd = List[Tuple[Float64, Float64]]()
        var i = 0
        while i < dirs.__len__() - 1:
            var a = dirs[i]
            var b = dirs[i + 1]
            nd.append(a)
            var mx = a[0] + b[0]
            var my = a[1] + b[1]
            var ml = sqrt_f64(mx * mx + my * my)
            if ml != 0.0:
                mx /= ml
                my /= ml
                nd.append((mx, my))
            i += 1
        nd.append(dirs[dirs.__len__() - 1])
        dirs = nd.copy()
        k += 1

    var j = 0
    while j < dirs.__len__():
        var d = dirs[j]
        out.append((cx + d[0] * r, cy + d[1] * r))
        j += 1


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


fn _disk(cx: Float64, cy: Float64, r: Float64, quad_segs: Int32) -> Polygon:
    if r <= 0.0:
        return _empty_polygon()
    var segs: Int = Int(quad_segs)
    if segs < 1:
        segs = 1
    # Approximate full circle using 8-direction base subdivided by quad_segs.
    var s = sqrt_f64(0.5)
    var base = List[Tuple[Float64, Float64]]()
    base.append((1.0, 0.0))
    base.append((s, s))
    base.append((0.0, 1.0))
    base.append((-s, s))
    base.append((-1.0, 0.0))
    base.append((-s, -s))
    base.append((0.0, -1.0))
    base.append((s, -s))

    var per_oct: Int = Int(segs / 2)
    if per_oct < 1:
        per_oct = 1

    var ring = List[Tuple[Float64, Float64]]()
    var bi = 0
    while bi < 8:
        var a = base[bi]
        var b = base[(bi + 1) % 8]
        var t: Int = 0
        while t <= per_oct:
            if bi != 0 and t == 0:
                t += 1
                continue
            var tt = Float64(t) / Float64(per_oct)
            var vx = a[0] * (1.0 - tt) + b[0] * tt
            var vy = a[1] * (1.0 - tt) + b[1] * tt
            var vl = sqrt_f64(vx * vx + vy * vy)
            if vl != 0.0:
                vx /= vl
                vy /= vl
            ring.append((cx + vx * r, cy + vy * r))
            t += 1
        bi += 1
    if ring.__len__() > 0:
        ring.append(ring[0])
    return Polygon(LinearRing(ring))


fn circle(cx: Float64, cy: Float64, radius: Float64, quad_segs: Int32 = 16) -> Geometry:
    return Geometry(_disk(cx, cy, radius, quad_segs))


fn circle(center: Point, radius: Float64, quad_segs: Int32 = 16) -> Geometry:
    return Geometry(_disk(center.x, center.y, radius, quad_segs))


fn buffer(geom: Geometry, _distance: Float64, _quad_segs: Int32 = 16) -> Geometry:
    if _distance <= 0.0:
        return geom.copy()
    if geom.is_linestring():
        return buffer(geom.as_linestring(), _distance, _quad_segs)
    if geom.is_multilinestring():
        return buffer(geom.as_multilinestring(), _distance, _quad_segs)
    if geom.is_polygon():
        return buffer(geom.as_polygon(), _distance, _quad_segs)
    if geom.is_multipolygon():
        return buffer(geom.as_multipolygon(), _distance, _quad_segs)
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
    if geom.is_polygon():
        return buffer(geom.as_polygon(), _distance, _quad_segs, cap_style, join_style, mitre_limit)
    if geom.is_multipolygon():
        return buffer(geom.as_multipolygon(), _distance, _quad_segs, cap_style, join_style, mitre_limit)
    return geom.copy()


fn buffer(ls: LineString, distance: Float64, _quad_segs: Int32 = 16) -> Geometry:
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

    # Robust buffering for round joins: Minkowski-sum style union of segment tubes
    # and vertex disks. This naturally trims corners and clips join circles.
    if join_style == JOIN_ROUND:
        var pts = ls.coords.copy()
        var n = pts.__len__()

        # Square caps extend endpoints along tangents.
        if cap_style == CAP_SQUARE:
            var t0 = _unit_tangent(pts[0][0], pts[0][1], pts[1][0], pts[1][1])
            pts[0] = (pts[0][0] - t0[0] * distance, pts[0][1] - t0[1] * distance)
            var t1 = _unit_tangent(pts[n - 2][0], pts[n - 2][1], pts[n - 1][0], pts[n - 1][1])
            pts[n - 1] = (pts[n - 1][0] + t1[0] * distance, pts[n - 1][1] + t1[1] * distance)

        var acc = Geometry(_empty_polygon())

        # Segment tubes
        var i = 0
        while i < n - 1:
            var tube = Geometry(_segment_tube(pts[i][0], pts[i][1], pts[i + 1][0], pts[i + 1][1], distance))
            acc = union(acc, tube)
            i += 1

        # Round join disks at internal vertices
        var k = 1
        while k < n - 1:
            var d = Geometry(_disk(pts[k][0], pts[k][1], distance, _quad_segs))
            acc = union(acc, d)
            k += 1

        # Round caps: endpoint disks. Flat/square caps: omit disks at ends.
        if cap_style == CAP_ROUND:
            acc = union(acc, Geometry(_disk(pts[0][0], pts[0][1], distance, _quad_segs)))
            acc = union(acc, Geometry(_disk(pts[n - 1][0], pts[n - 1][1], distance, _quad_segs)))

        return acc.copy()

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
        var li = _line_intersection_tu(p0x, p0y, t0x, t0y, p1x, p1y, t1x, t1y)
        var lpt = li[0]
        var _ = li[1]
        var _ = li[2]
        var lok = li[3]

        # Right side intersection (use -normals)
        var q0x = px - n0x * distance
        var q0y = py - n0y * distance
        var q1x = px - n1x * distance
        var q1y = py - n1y * distance
        var ri = _line_intersection_tu(q0x, q0y, t0x, t0y, q1x, q1y, t1x, t1y)
        var rpt = ri[0]
        var _ = ri[1]
        var _ = ri[2]
        var rok = ri[3]

        if join_style == JOIN_ROUND:
            # Only add the round arc on the *outer* (convex) side of the turn.
            # On the inner (concave) side, connect with the offset-line intersection.
            var cr = _cross(t0x, t0y, t1x, t1y)
            if cr > 0.0:
                # Left turn: left side is outer.
                _append_arc_join(left, px, py, distance, _quad_segs, n0x, n0y, n1x, n1y)
                if rok:
                    right.append(rpt)
                else:
                    right.append((q0x, q0y))
                    right.append((q1x, q1y))
            elif cr < 0.0:
                # Right turn: right side is outer.
                _append_arc_join(right, px, py, distance, _quad_segs, -n0x, -n0y, -n1x, -n1y)
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
            if join_style == JOIN_BEVEL:
                var cr = _cross(t0x, t0y, t1x, t1y)
                if cr > 0.0:
                    # Left turn: right side is outer (convex): bevel with endpoints.
                    # Left side is inner (concave): use intersection only (no extra points).
                    right.append((q0x, q0y))
                    right.append((q1x, q1y))
                    if lok:
                        left.append(lpt)
                elif cr < 0.0:
                    # Right turn: left side is outer (convex): bevel with endpoints.
                    # Right side is inner (concave): use intersection only (no extra points).
                    left.append((p0x, p0y))
                    left.append((p1x, p1y))
                    if rok:
                        right.append(rpt)
                else:
                    left.append((p0x, p0y))
                    left.append((p1x, p1y))
                    right.append((q0x, q0y))
                    right.append((q1x, q1y))
            else:
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
        # arc from +n to -n around end (exterior)
        _append_arc(ring, pts[n - 1][0], pts[n - 1][1], distance, _quad_segs, nx, ny, -nx, -ny, False, True)

    var rr = right.__len__() - 1
    while rr >= 0:
        ring.append(right[rr])
        rr -= 1

    # start cap
    if cap_style == CAP_ROUND:
        var nx0 = nxs[0]
        var ny0 = nys[0]
        # arc from -n to +n around start (exterior)
        _append_arc(ring, pts[0][0], pts[0][1], distance, _quad_segs, -nx0, -ny0, nx0, ny0, False, True)

    if ring.__len__() > 0:
        var first = ring[0]
        var last = ring[ring.__len__() - 1]
        if first[0] != last[0] or first[1] != last[1]:
            ring.append(first)

    return Geometry(Polygon(LinearRing(ring)))


fn buffer(p: Polygon, distance: Float64, quad_segs: Int32 = 16) -> Geometry:
    return buffer(p, distance, quad_segs, CAP_ROUND, JOIN_ROUND, 5.0)


fn buffer(
    p: Polygon,
    distance: Float64,
    quad_segs: Int32,
    cap_style: CapStyle,
    join_style: JoinStyle,
    mitre_limit: Float64 = 5.0,
) -> Geometry:
    if distance <= 0.0:
        return Geometry(p.copy())
    if p.is_empty():
        return Geometry(_empty_polygon())

    # Best-effort polygon buffer using boundary buffering.
    # Expand exterior and shrink/fill holes by unioning with buffered boundary rings.
    var acc = Geometry(p.copy())

    var shell_ls = LineString(p.shell.coords.copy())
    acc = union(acc, buffer(shell_ls, distance, quad_segs, cap_style, join_style, mitre_limit))

    for h in p.holes:
        var hole_ls = LineString(h.coords.copy())
        acc = union(acc, buffer(hole_ls, distance, quad_segs, cap_style, join_style, mitre_limit))

    return acc.copy()


fn buffer(mp: MultiPolygon, distance: Float64, quad_segs: Int32 = 16) -> Geometry:
    return buffer(mp, distance, quad_segs, CAP_ROUND, JOIN_ROUND, 5.0)


fn buffer(
    mp: MultiPolygon,
    distance: Float64,
    quad_segs: Int32,
    cap_style: CapStyle,
    join_style: JoinStyle,
    mitre_limit: Float64 = 5.0,
) -> Geometry:
    if distance <= 0.0:
        return Geometry(mp.copy())
    var polys = List[Polygon]()
    for p in mp.polys:
        var g = buffer(p.copy(), distance, quad_segs, cap_style, join_style, mitre_limit)
        if g.is_polygon() and not g.is_empty():
            polys.append(g.as_polygon())
        elif g.is_multipolygon():
            for pp in g.as_multipolygon().polys:
                polys.append(pp.copy())
    return Geometry(MultiPolygon(polys))


fn buffer(mls: MultiLineString, distance: Float64, quad_segs: Int32 = 16) -> Geometry:
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
    if geom.is_linestring():
        return Geometry(_simplify_linestring(geom.as_linestring(), _tolerance))
    if geom.is_multilinestring():
        var mls = geom.as_multilinestring()
        var out = List[LineString]()
        for ln in mls.lines:
            out.append(_simplify_linestring(ln.copy(), _tolerance))
        return Geometry(MultiLineString(out))
    if geom.is_polygon():
        var poly = geom.as_polygon()
        var shell_coords = _simplify_ring_coords(poly.shell.coords, _tolerance)
        if shell_coords.__len__() < 4:
            return Geometry(_empty_polygon())
        var holes = List[LinearRing]()
        for h in poly.holes:
            var hc = _simplify_ring_coords(h.coords, _tolerance)
            if hc.__len__() >= 4:
                holes.append(LinearRing(hc))
        return Geometry(Polygon(LinearRing(shell_coords), holes))
    if geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        var polys = List[Polygon]()
        for p in mp.polys:
            var g = simplify(Geometry(p.copy()), _tolerance, _preserve_topology)
            if g.is_polygon() and not g.is_empty():
                polys.append(g.as_polygon())
        return Geometry(MultiPolygon(polys))
    return geom.copy()


fn convex_hull(geom: Geometry) -> Geometry:
    var pts = List[Tuple[Float64, Float64]]()
    _collect_coords(geom.copy(), pts)
    if pts.__len__() == 0:
        return Geometry(GeometryCollection([]))
    var ring = _hull_ring(pts)
    if ring.__len__() == 1:
        return Geometry(Point(ring[0][0], ring[0][1]))
    if ring.__len__() == 2:
        return Geometry(LineString([ring[0], ring[1]]))
    if ring.__len__() > 2:
        var out = ring.copy()
        var first = out[0]
        var last = out[out.__len__() - 1]
        if first[0] != last[0] or first[1] != last[1]:
            out.append(first)
        return Geometry(Polygon(LinearRing(out)))
    return Geometry(GeometryCollection([]))
