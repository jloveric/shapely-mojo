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
from shapely.set_operations import union, difference
from shapely.ops import polygonize_full
from shapely.algorithms import signed_area_coords, segments_intersect, point_in_ring


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


fn circle(cx: Float64, cy: Float64, r: Float64, quad_segs: Int32 = 8) -> Geometry:
    _ = quad_segs
    return Geometry(_circle_polygon(cx, cy, r))


fn _unit_tangent(ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Tuple[Float64, Float64]:
    var dx = bx - ax
    var dy = by - ay
    var len = sqrt_f64(dx * dx + dy * dy)
    if len == 0.0:
        return (0.0, 0.0)
    return (dx / len, dy / len)


fn _unit_normal_left(tx: Float64, ty: Float64) -> Tuple[Float64, Float64]:
    return (-ty, tx)


fn _dot(ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
    return ax * bx + ay * by


fn _cross(ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
    return ax * by - ay * bx


fn _close_ring(mut coords: List[Tuple[Float64, Float64]]):
    if coords.__len__() == 0:
        return
    var first = coords[0]
    var last = coords[coords.__len__() - 1]
    if first[0] != last[0] or first[1] != last[1]:
        coords.append(first)


fn _reverse_coords(coords: List[Tuple[Float64, Float64]]) -> List[Tuple[Float64, Float64]]:
    var out = List[Tuple[Float64, Float64]]()
    var i: Int = coords.__len__() - 1
    while i >= 0:
        out.append(coords[i])
        i -= 1
    _close_ring(out)
    return out.copy()


fn _offset_ring_inward_round(
    coords: List[Tuple[Float64, Float64]],
    distance: Float64,
    quad_segs: Int32,
) -> List[Tuple[Float64, Float64]]:
    var out = List[Tuple[Float64, Float64]]()
    if coords.__len__() < 4:
        return out.copy()

    var area = signed_area_coords(coords)
    if area == 0.0:
        return out.copy()
    var outward_right = area > 0.0

    var sign = 1.0
    if not outward_right:
        sign = -1.0

    var m = coords.__len__() - 1

    var txs = List[Float64]()
    var tys = List[Float64]()
    var nxs = List[Float64]()
    var nys = List[Float64]()
    var i = 0
    while i < m:
        var j = i + 1
        if j == m:
            j = 0
        var t = _unit_tangent(coords[i][0], coords[i][1], coords[j][0], coords[j][1])
        txs.append(t[0])
        tys.append(t[1])
        var nn = _unit_normal_left(t[0], t[1])
        nxs.append(nn[0])
        nys.append(nn[1])
        i += 1

    var k = 0
    while k < m:
        var px = coords[k][0]
        var py = coords[k][1]

        var prev = k - 1
        if prev < 0:
            prev = m - 1

        var n0x = nxs[prev]
        var n0y = nys[prev]
        var n1x = nxs[k]
        var n1y = nys[k]
        var t0x = txs[prev]
        var t0y = tys[prev]
        var t1x = txs[k]
        var t1y = tys[k]

        var o0x = px + n0x * distance * sign
        var o0y = py + n0y * distance * sign
        var o1x = px + n1x * distance * sign
        var o1y = py + n1y * distance * sign
        var sx = n0x * sign
        var sy = n0y * sign
        var ex = n1x * sign
        var ey = n1y * sign

        var li = _line_intersection_tu(o0x, o0y, t0x, t0y, o1x, o1y, t1x, t1y)
        var ipt = li[0]
        var lok = li[3]

        var cr = _cross(t0x, t0y, t1x, t1y)
        var want_arc = False
        if outward_right:
            if cr < 0.0:
                want_arc = True
        else:
            if cr > 0.0:
                want_arc = True

        if want_arc:
            _append_arc_join(out, px, py, distance, quad_segs, sx, sy, ex, ey)
        else:
            if lok:
                out.append(ipt)
            else:
                out.append((o0x, o0y))
                out.append((o1x, o1y))

        k += 1

    _close_ring(out)
    if signed_area_coords(out) < 0.0:
        out = _reverse_coords(out)
    return out.copy()


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
            var _ = lower.pop()
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
            var _ = upper.pop()
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
        var _ = stack_a.pop()
        var _2 = stack_b.pop()

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
        var _ = open.pop()
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
) -> Tuple[Tuple[Float64, Float64], Bool]:
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
) -> Tuple[Tuple[Float64, Float64], Float64, Float64, Bool]:
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
        var first = ring[0]
        var last = ring[ring.__len__() - 1]
        if first[0] == last[0] and first[1] == last[1]:
            var _ = ring.pop()
        ring.append(ring[0])
    return Polygon(LinearRing(ring))


fn _is_convex_closed_ring(coords: List[Tuple[Float64, Float64]]) -> Bool:
    if coords.__len__() < 4:
        return False

    var first = coords[0]
    var last = coords[coords.__len__() - 1]
    if first[0] != last[0] or first[1] != last[1]:
        return False
    var area = signed_area_coords(coords)
    if area == 0.0:
        return False
    var m = coords.__len__() - 1
    var i = 0
    while i < m:
        var im1 = i - 1
        if im1 < 0:
            im1 = m - 1
        var ip1 = i + 1
        if ip1 == m:
            ip1 = 0
        var ax = coords[i][0] - coords[im1][0]
        var ay = coords[i][1] - coords[im1][1]
        var bx = coords[ip1][0] - coords[i][0]
        var by = coords[ip1][1] - coords[i][1]
        var cr = _cross(ax, ay, bx, by)
        if area > 0.0:
            if cr < 0.0:
                return False
        else:
            if cr > 0.0:
                return False
        i += 1
    return True


fn _dedup_consecutive_closed(
    coords: List[Tuple[Float64, Float64]]
) -> List[Tuple[Float64, Float64]]:
    var out = List[Tuple[Float64, Float64]]()
    var i = 0
    while i < coords.__len__():
        if out.__len__() == 0:
            out.append(coords[i])
        else:
            var prev = out[out.__len__() - 1]
            var cur = coords[i]
            if prev[0] != cur[0] or prev[1] != cur[1]:
                out.append(cur)
        i += 1
    if out.__len__() >= 2:
        var f = out[0]
        var l = out[out.__len__() - 1]
        if f[0] != l[0] or f[1] != l[1]:
            out.append(f)
    return out.copy()


fn _coords_bbox(
    coords: List[Tuple[Float64, Float64]]
) -> Tuple[Float64, Float64, Float64, Float64]:
    if coords.__len__() == 0:
        return (0.0, 0.0, 0.0, 0.0)
    var minx = coords[0][0]
    var miny = coords[0][1]
    var maxx = minx
    var maxy = miny
    var i = 1
    while i < coords.__len__():
        var c = coords[i]
        if c[0] < minx:
            minx = c[0]
        if c[0] > maxx:
            maxx = c[0]
        if c[1] < miny:
            miny = c[1]
        if c[1] > maxy:
            maxy = c[1]
        i += 1
    return (minx, miny, maxx, maxy)


fn _rings_intersect(
    a: List[Tuple[Float64, Float64]],
    b: List[Tuple[Float64, Float64]],
) -> Bool:
    if a.__len__() < 2 or b.__len__() < 2:
        return False
    var ab = _coords_bbox(a)
    var bb = _coords_bbox(b)
    if ab[2] < bb[0] or bb[2] < ab[0] or ab[3] < bb[1] or bb[3] < ab[1]:
        return False

    var an = a.__len__() - 1
    var bn = b.__len__() - 1
    var i = 0
    while i < an:
        var a1 = a[i]
        var a2 = a[i + 1]
        var aminx = a1[0] if a1[0] < a2[0] else a2[0]
        var amaxx = a2[0] if a2[0] > a1[0] else a1[0]
        var aminy = a1[1] if a1[1] < a2[1] else a2[1]
        var amaxy = a2[1] if a2[1] > a1[1] else a1[1]
        var j = 0
        while j < bn:
            var b1 = b[j]
            var b2 = b[j + 1]
            var bminx = b1[0] if b1[0] < b2[0] else b2[0]
            var bmaxx = b2[0] if b2[0] > b1[0] else b1[0]
            var bminy = b1[1] if b1[1] < b2[1] else b2[1]
            var bmaxy = b2[1] if b2[1] > b1[1] else b1[1]
            if amaxx < bminx or bmaxx < aminx or amaxy < bminy or bmaxy < aminy:
                j += 1
                continue
            if segments_intersect(a1, a2, b1, b2):
                return True
            j += 1
        i += 1
    return False


fn _ring_has_self_intersection(coords: List[Tuple[Float64, Float64]]) -> Bool:
    if coords.__len__() < 4:
        return False
    var n = coords.__len__() - 1
    var i = 0
    while i < n:
        var a1 = coords[i]
        var a2 = coords[i + 1]
        var aminx = a1[0] if a1[0] < a2[0] else a2[0]
        var amaxx = a2[0] if a2[0] > a1[0] else a1[0]
        var aminy = a1[1] if a1[1] < a2[1] else a2[1]
        var amaxy = a2[1] if a2[1] > a1[1] else a1[1]
        var j = i + 1
        while j < n:
            if j == i:
                j += 1
                continue
            if j == i + 1:
                j += 1
                continue
            if i == 0 and j == n - 1:
                j += 1
                continue
            var b1 = coords[j]
            var b2 = coords[j + 1]
            var bminx = b1[0] if b1[0] < b2[0] else b2[0]
            var bmaxx = b2[0] if b2[0] > b1[0] else b1[0]
            var bminy = b1[1] if b1[1] < b2[1] else b2[1]
            var bmaxy = b2[1] if b2[1] > b1[1] else b1[1]
            if amaxx < bminx or bmaxx < aminx or amaxy < bminy or bmaxy < aminy:
                j += 1
                continue
            if segments_intersect(a1, a2, b1, b2):
                return True
            j += 1
        i += 1
    return False


fn _offset_ring_outward_round(
    coords: List[Tuple[Float64, Float64]],
    distance: Float64,
    quad_segs: Int32,
) -> List[Tuple[Float64, Float64]]:
    var out = List[Tuple[Float64, Float64]]()
    if coords.__len__() < 4:
        return out.copy()

    var area = signed_area_coords(coords)
    if area == 0.0:
        return out.copy()
    var outward_right = area > 0.0

    var m = coords.__len__() - 1

    var txs = List[Float64]()
    var tys = List[Float64]()
    var nxs = List[Float64]()
    var nys = List[Float64]()
    var i = 0
    while i < m:
        var j = i + 1
        if j == m:
            j = 0
        var t = _unit_tangent(coords[i][0], coords[i][1], coords[j][0], coords[j][1])
        txs.append(t[0])
        tys.append(t[1])
        var nn = _unit_normal_left(t[0], t[1])
        nxs.append(nn[0])
        nys.append(nn[1])
        i += 1

    var k = 0
    while k < m:
        var px = coords[k][0]
        var py = coords[k][1]

        var prev = k - 1
        if prev < 0:
            prev = m - 1

        var n0x = nxs[prev]
        var n0y = nys[prev]
        var n1x = nxs[k]
        var n1y = nys[k]
        var t0x = txs[prev]
        var t0y = tys[prev]
        var t1x = txs[k]
        var t1y = tys[k]

        var o0x = px + n0x * distance
        var o0y = py + n0y * distance
        var o1x = px + n1x * distance
        var o1y = py + n1y * distance
        var sx = n0x
        var sy = n0y
        var ex = n1x
        var ey = n1y
        if outward_right:
            o0x = px - n0x * distance
            o0y = py - n0y * distance
            o1x = px - n1x * distance
            o1y = py - n1y * distance
            sx = -n0x
            sy = -n0y
            ex = -n1x
            ey = -n1y

        var li = _line_intersection_tu(o0x, o0y, t0x, t0y, o1x, o1y, t1x, t1y)
        var ipt = li[0]
        var lok = li[3]

        var cr = _cross(t0x, t0y, t1x, t1y)
        var want_arc = False
        if outward_right:
            if cr > 0.0:
                want_arc = True
        else:
            if cr < 0.0:
                want_arc = True

        if want_arc:
            _append_arc_join(out, px, py, distance, quad_segs, sx, sy, ex, ey)
        else:
            if lok:
                out.append(ipt)
            else:
                out.append((o0x, o0y))
                out.append((o1x, o1y))

        k += 1

    _close_ring(out)
    if signed_area_coords(out) < 0.0:
        out = _reverse_coords(out)
    return out.copy()


fn _is_axis_aligned_box_shell(
    coords: List[Tuple[Float64, Float64]],
    xmin: Float64,
    ymin: Float64,
    xmax: Float64,
    ymax: Float64,
) -> Bool:
    if coords.__len__() != 5:
        return False
    var c0 = coords[0]
    var c1 = coords[1]
    var c2 = coords[2]
    var c3 = coords[3]
    var c4 = coords[4]
    if c0[0] != c4[0] or c0[1] != c4[1]:
        return False
    if (
        c0[0] == xmax
        and c0[1] == ymin
        and c1[0] == xmax
        and c1[1] == ymax
        and c2[0] == xmin
        and c2[1] == ymax
        and c3[0] == xmin
        and c3[1] == ymin
    ):
        return True
    if (
        c0[0] == xmin
        and c0[1] == ymin
        and c1[0] == xmin
        and c1[1] == ymax
        and c2[0] == xmax
        and c2[1] == ymax
        and c3[0] == xmax
        and c3[1] == ymin
    ):
        return True
    return False


fn _rounded_rect_polygon(
    xmin: Float64,
    ymin: Float64,
    xmax: Float64,
    ymax: Float64,
    r: Float64,
    quad_segs: Int32,
) -> Polygon:
    if r <= 0.0:
        return Polygon(LinearRing([(xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin), (xmax, ymin)]))

    var ring = List[Tuple[Float64, Float64]]()
    _append_arc(ring, xmax, ymax, r, quad_segs, 1.0, 0.0, 0.0, 1.0, True, True)
    ring.append((xmin, ymax + r))
    _append_arc(ring, xmin, ymax, r, quad_segs, 0.0, 1.0, -1.0, 0.0, True, True)
    ring.append((xmin - r, ymin))
    _append_arc(ring, xmin, ymin, r, quad_segs, -1.0, 0.0, 0.0, -1.0, True, True)
    ring.append((xmax, ymin - r))
    _append_arc(ring, xmax, ymin, r, quad_segs, 0.0, -1.0, 1.0, 0.0, True, True)
    ring.append((xmax + r, ymax))

    if ring.__len__() > 0:
        var first = ring[0]
        var last = ring[ring.__len__() - 1]
        if first[0] != last[0] or first[1] != last[1]:
            ring.append(first)
    return Polygon(LinearRing(ring))


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

    var is_closed = False
    if ls.coords.__len__() >= 4:
        var f = ls.coords[0]
        var l = ls.coords[ls.coords.__len__() - 1]
        if f[0] == l[0] and f[1] == l[1]:
            is_closed = True

    var effective_cap: CapStyle = cap_style
    if is_closed:
        # Closed rings don't have caps.
        effective_cap = CAP_FLAT

    # Robust buffering for round joins: Minkowski-sum style union of segment tubes
    # and vertex disks. This naturally trims corners and clips join circles.
    if join_style == JOIN_ROUND:
        var pts = ls.coords.copy()
        var n = pts.__len__()

        # Square caps extend endpoints along tangents.
        if effective_cap == CAP_SQUARE and not is_closed:
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

        # Closed rings need a join disk at the closure vertex as well.
        if is_closed:
            acc = union(acc, Geometry(_disk(pts[0][0], pts[0][1], distance, _quad_segs)))

        # Round caps: endpoint disks. Flat/square caps: omit disks at ends.
        if effective_cap == CAP_ROUND and not is_closed:
            acc = union(acc, Geometry(_disk(pts[0][0], pts[0][1], distance, _quad_segs)))
            acc = union(acc, Geometry(_disk(pts[n - 1][0], pts[n - 1][1], distance, _quad_segs)))

        return acc.copy()

    # Build a single polygon ring for the polyline buffer (no unions).
    var n = ls.coords.__len__()

    # Optionally extend endpoints for square caps (open polylines only)
    var pts = ls.coords.copy()
    if effective_cap == CAP_SQUARE and not is_closed:
        var t0 = _unit_tangent(pts[0][0], pts[0][1], pts[1][0], pts[1][1])
        pts[0] = (pts[0][0] - t0[0] * distance, pts[0][1] - t0[1] * distance)
        var t1 = _unit_tangent(pts[n - 2][0], pts[n - 2][1], pts[n - 1][0], pts[n - 1][1])
        pts[n - 1] = (pts[n - 1][0] + t1[0] * distance, pts[n - 1][1] + t1[1] * distance)

    if is_closed:
        # Closed ring: all vertices are internal and segment data wraps around.
        var m = n - 1

        var txs = List[Float64]()
        var tys = List[Float64]()
        var nxs = List[Float64]()
        var nys = List[Float64]()
        var i = 0
        while i < m:
            var j = i + 1
            if j == m:
                j = 0
            var t = _unit_tangent(pts[i][0], pts[i][1], pts[j][0], pts[j][1])
            txs.append(t[0])
            tys.append(t[1])
            var nn = _unit_normal_left(t[0], t[1])
            nxs.append(nn[0])
            nys.append(nn[1])
            i += 1

        var left = List[Tuple[Float64, Float64]]()
        var right = List[Tuple[Float64, Float64]]()

        var k = 0
        while k < m:
            var px = pts[k][0]
            var py = pts[k][1]

            var prev = k - 1
            if prev < 0:
                prev = m - 1

            var n0x = nxs[prev]
            var n0y = nys[prev]
            var n1x = nxs[k]
            var n1y = nys[k]
            var t0x = txs[prev]
            var t0y = tys[prev]
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

            var force_bevel = False
            if join_style == JOIN_MITRE:
                var can_mitre = lok and rok
                if can_mitre:
                    # Apply a simple mitre limit: fallback to bevel if too far
                    var dlx = lpt[0] - px
                    var dly = lpt[1] - py
                    var drx = rpt[0] - px
                    var dry = rpt[1] - py
                    var ml = sqrt_f64(dlx * dlx + dly * dly)
                    var mr = sqrt_f64(drx * drx + dry * dry)
                    if ml <= _mitre_limit * distance and mr <= _mitre_limit * distance:
                        left.append(lpt)
                        right.append(rpt)
                    else:
                        force_bevel = True
                else:
                    force_bevel = True

            if join_style == JOIN_BEVEL or force_bevel:
                # bevel (or mitre fallback)
                var cr = _cross(t0x, t0y, t1x, t1y)
                if cr > 0.0:
                    # Left turn: right side is outer (convex): bevel with endpoints.
                    # Left side is inner (concave): use forward intersection only.
                    right.append((q0x, q0y))
                    right.append((q1x, q1y))
                    if lok:
                        var dlx = lpt[0] - px
                        var dly = lpt[1] - py
                        var ml = sqrt_f64(dlx * dlx + dly * dly)
                        if ml <= _mitre_limit * distance:
                            left.append(lpt)
                        else:
                            left.append((p0x, p0y))
                            left.append((p1x, p1y))
                    else:
                        left.append((p0x, p0y))
                        left.append((p1x, p1y))
                elif cr < 0.0:
                    # Right turn: left side is outer (convex): bevel with endpoints.
                    # Right side is inner (concave): use forward intersection only.
                    left.append((p0x, p0y))
                    left.append((p1x, p1y))
                    if rok:
                        var drx = rpt[0] - px
                        var dry = rpt[1] - py
                        var mr = sqrt_f64(drx * drx + dry * dry)
                        if mr <= _mitre_limit * distance:
                            right.append(rpt)
                        else:
                            right.append((q0x, q0y))
                            right.append((q1x, q1y))
                    else:
                        right.append((q0x, q0y))
                        right.append((q1x, q1y))
                else:
                    left.append((p0x, p0y))
                    left.append((p1x, p1y))
                    right.append((q0x, q0y))
                    right.append((q1x, q1y))
            elif join_style != JOIN_MITRE:
                # parallel/unknown fallback
                left.append((p0x, p0y))
                left.append((p1x, p1y))
                right.append((q0x, q0y))
                right.append((q1x, q1y))

            k += 1

        _close_ring(left)
        _close_ring(right)

        var ring_area = signed_area_coords(ls.coords)

        var outer_coords = List[Tuple[Float64, Float64]]()
        var inner_coords = List[Tuple[Float64, Float64]]()
        if ring_area > 0.0:
            # CCW ring: outward is to the right.
            for p in right:
                outer_coords.append(p)
            for p in left:
                inner_coords.append(p)
        else:
            # CW ring: outward is to the left.
            for p in left:
                outer_coords.append(p)
            for p in right:
                inner_coords.append(p)

        _close_ring(outer_coords)
        _close_ring(inner_coords)

        if signed_area_coords(outer_coords) < 0.0:
            outer_coords = _reverse_coords(outer_coords)
        if signed_area_coords(inner_coords) > 0.0:
            inner_coords = _reverse_coords(inner_coords)

        # Build as a boolean difference for robustness (prevents occasional seams on
        # circle-like rings where the direct shell+hole polygon can be invalid).
        var outer_poly = Polygon(LinearRing(outer_coords))
        var inner_poly = Polygon(LinearRing(inner_coords))
        return difference(outer_poly, inner_poly)

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

    if join_style == JOIN_ROUND and p.holes.__len__() == 0:
        var b = p.bounds()
        if _is_axis_aligned_box_shell(p.shell.coords, b[0], b[1], b[2], b[3]):
            return Geometry(_rounded_rect_polygon(b[0], b[1], b[2], b[3], distance, quad_segs))
        if _is_convex_closed_ring(p.shell.coords):
            var ring = _offset_ring_outward_round(p.shell.coords, distance, quad_segs)
            if ring.__len__() > 0:
                return Geometry(Polygon(LinearRing(ring)))
        var ring2 = _offset_ring_outward_round(p.shell.coords, distance, quad_segs)
        if ring2.__len__() >= 4:
            var ring3 = List[Tuple[Float64, Float64]]()
            var k2 = 0
            while k2 < ring2.__len__():
                if ring3.__len__() == 0:
                    ring3.append(ring2[k2])
                else:
                    var prev = ring3[ring3.__len__() - 1]
                    var cur = ring2[k2]
                    if prev[0] != cur[0] or prev[1] != cur[1]:
                        ring3.append(cur)
                k2 += 1
            if ring3.__len__() >= 2:
                var f = ring3[0]
                var l = ring3[ring3.__len__() - 1]
                if f[0] != l[0] or f[1] != l[1]:
                    ring3.append(f)
            if ring3.__len__() >= 4 and not _ring_has_self_intersection(ring3):
                return Geometry(Polygon(LinearRing(ring3)))
            var segs = List[LineString]()
            var i2 = 0
            while i2 < ring2.__len__() - 1:
                var pts = List[Tuple[Float64, Float64]]()
                pts.append(ring2[i2])
                pts.append(ring2[i2 + 1])
                segs.append(LineString(pts))
                i2 += 1
            if segs.__len__() > 0:
                var mls = MultiLineString(segs)
                var res = polygonize_full(Geometry(mls.copy()))
                ref polys = res[0]
                var best_area = -1.0
                var best = _empty_polygon()
                for g in polys.geoms:
                    if g.is_polygon():
                        var pp = g.as_polygon()
                        var a = pp.area()
                        if a > best_area:
                            best_area = a
                            best = pp.copy()
                    elif g.is_multipolygon():
                        for pp in g.as_multipolygon().polys:
                            var a = pp.area()
                            if a > best_area:
                                best_area = a
                                best = pp.copy()
                if best_area > 0.0:
                    return Geometry(best.copy())

    if join_style == JOIN_ROUND and p.holes.__len__() > 0:
        var shell_ring = _offset_ring_outward_round(p.shell.coords, distance, quad_segs)
        var shell_ring2 = _dedup_consecutive_closed(shell_ring)
        if shell_ring2.__len__() >= 4 and not _ring_has_self_intersection(shell_ring2):
            var shell_lr = LinearRing(shell_ring2)
            var holes_out = List[LinearRing]()
            var hole_coords = List[List[Tuple[Float64, Float64]]]()

            var ok = True
            var hi = 0
            while hi < p.holes.__len__():
                ref h = p.holes[hi]
                var hr = _offset_ring_inward_round(h.coords, distance, quad_segs)
                var hr2 = _dedup_consecutive_closed(hr)
                if hr2.__len__() < 4:
                    hi += 1
                    continue
                if _ring_has_self_intersection(hr2):
                    ok = False
                    break
                var pt = hr2[0]
                var pin = point_in_ring(Point(pt[0], pt[1]), shell_lr)
                if pin != 1:
                    ok = False
                    break

                hole_coords.append(hr2.copy())
                holes_out.append(LinearRing(hr2))
                hi += 1

            if ok and hole_coords.__len__() > 1:
                var a = 0
                while a < hole_coords.__len__():
                    var b = a + 1
                    while b < hole_coords.__len__():
                        ref ra = hole_coords[a]
                        ref rb = hole_coords[b]
                        if _rings_intersect(ra, rb):
                            ok = False
                            break
                        var pta = ra[0]
                        var ptb = rb[0]
                        if point_in_ring(Point(pta[0], pta[1]), LinearRing(rb)) == 1:
                            ok = False
                            break
                        if point_in_ring(Point(ptb[0], ptb[1]), LinearRing(ra)) == 1:
                            ok = False
                            break
                        b += 1
                    if not ok:
                        break
                    a += 1

            if ok:
                return Geometry(Polygon(shell_lr, holes_out))

    # Best-effort polygon buffer using boundary buffering.
    # Expand exterior and shrink/fill holes by unioning with buffered boundary rings.
    var acc = Geometry(p.copy())

    # Polygon rings are closed; cap style does not apply. Force flat caps to avoid
    # spurious end-cap artifacts when buffering closed rings.
    var ring_cap: CapStyle = CAP_FLAT

    var shell_ls = LineString(p.shell.coords.copy())
    acc = union(acc, buffer(shell_ls, distance, quad_segs, ring_cap, join_style, mitre_limit))

    for h in p.holes:
        var hole_ls = LineString(h.coords.copy())
        acc = union(acc, buffer(hole_ls, distance, quad_segs, ring_cap, join_style, mitre_limit))

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
