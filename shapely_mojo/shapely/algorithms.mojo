from shapely.geometry import Point, LineString, LinearRing, Polygon


fn abs_f64(x: Float64) -> Float64:
    if x < 0.0:
        return -x
    return x


fn orientation(
    ax: Float64, ay: Float64, bx: Float64, by: Float64, cx: Float64, cy: Float64
) -> Int32:
    var v = (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
    if v > 0.0:
        return 1
    if v < 0.0:
        return -1
    return 0


fn on_segment(
    ax: Float64, ay: Float64, bx: Float64, by: Float64, cx: Float64, cy: Float64
) -> Bool:
    if orientation(ax, ay, bx, by, cx, cy) != 0:
        return False
    var minx = ax if ax < bx else bx
    var maxx = ax if ax > bx else bx
    var miny = ay if ay < by else by
    var maxy = ay if ay > by else by
    return cx >= minx and cx <= maxx and cy >= miny and cy <= maxy


fn segments_intersect(
    a1: Tuple[Float64, Float64],
    a2: Tuple[Float64, Float64],
    b1: Tuple[Float64, Float64],
    b2: Tuple[Float64, Float64],
) -> Bool:
    var o1 = orientation(a1[0], a1[1], a2[0], a2[1], b1[0], b1[1])
    var o2 = orientation(a1[0], a1[1], a2[0], a2[1], b2[0], b2[1])
    var o3 = orientation(b1[0], b1[1], b2[0], b2[1], a1[0], a1[1])
    var o4 = orientation(b1[0], b1[1], b2[0], b2[1], a2[0], a2[1])

    if o1 != o2 and o3 != o4:
        return True
    if o1 == 0 and on_segment(a1[0], a1[1], a2[0], a2[1], b1[0], b1[1]):
        return True
    if o2 == 0 and on_segment(a1[0], a1[1], a2[0], a2[1], b2[0], b2[1]):
        return True
    if o3 == 0 and on_segment(b1[0], b1[1], b2[0], b2[1], a1[0], a1[1]):
        return True
    if o4 == 0 and on_segment(b1[0], b1[1], b2[0], b2[1], a2[0], a2[1]):
        return True
    return False


fn point_on_linestring(p: Point, ls: LineString) -> Bool:
    if ls.coords.__len__() < 2:
        return False
    var px = p.x
    var py = p.y
    for i in range(0, ls.coords.__len__() - 1):
        var a = ls.coords[i]
        var b = ls.coords[i + 1]
        if on_segment(a[0], a[1], b[0], b[1], px, py):
            return True
    return False


fn point_in_ring(pt: Point, ring: LinearRing) -> Int32:
    # returns: 2 if on boundary, 1 if inside, 0 if outside
    var n = ring.coords.__len__()
    if n < 2:
        return 0
    var inside = False
    var j = n - 1
    var px = pt.x
    var py = pt.y
    for i in range(0, n):
        var xi = ring.coords[i][0]
        var yi = ring.coords[i][1]
        var xj = ring.coords[j][0]
        var yj = ring.coords[j][1]
        # boundary check
        if on_segment(xj, yj, xi, yi, px, py):
            return 2
        var intersect = ((yi > py) != (yj > py)) and (
            px
            < (xj - xi) * (py - yi) / ((yj - yi) if (yj - yi) != 0.0 else 1e-18)
            + xi
        )
        if intersect:
            inside = not inside
        j = i
    return 1 if inside else 0


fn point_in_polygon(pt: Point, poly: Polygon) -> Int32:
    ref shell = poly.shell
    var shell_res = point_in_ring(pt, shell)
    if shell_res == 0:
        return 0
    if shell_res == 2:
        return 2
    for h in poly.holes:
        var res = point_in_ring(pt, h)
        if res == 2:
            return 2
        if res == 1:
            return 0
    return 1


fn any_segment_intersection(a: LineString, b: LineString) -> Bool:
    if a.coords.__len__() < 2 or b.coords.__len__() < 2:
        return False
    for i in range(0, a.coords.__len__() - 1):
        var a1 = a.coords[i]
        var a2 = a.coords[i + 1]
        for j in range(0, b.coords.__len__() - 1):
            var b1 = b.coords[j]
            var b2 = b.coords[j + 1]
            if segments_intersect(a1, a2, b1, b2):
                return True
    return False


fn any_segment_intersection_coords(
    a_coords: List[Tuple[Float64, Float64]],
    b_coords: List[Tuple[Float64, Float64]],
) -> Bool:
    if a_coords.__len__() < 2 or b_coords.__len__() < 2:
        return False
    for i in range(0, a_coords.__len__() - 1):
        var a1 = a_coords[i]
        var a2 = a_coords[i + 1]
        for j in range(0, b_coords.__len__() - 1):
            var b1 = b_coords[j]
            var b2 = b_coords[j + 1]
            if segments_intersect(a1, a2, b1, b2):
                return True
    return False


fn signed_area_coords(coords: List[Tuple[Float64, Float64]]) -> Float64:
    if coords.__len__() < 2:
        return 0.0
    var s = 0.0
    for i in range(0, coords.__len__() - 1):
        var a = coords[i]
        var b = coords[i + 1]
        s += a[0] * b[1] - a[1] * b[0]
    return 0.5 * s


fn ring_is_ccw(r: LinearRing) -> Bool:
    return signed_area_coords(r.coords) > 0.0


fn segment_intersections(
    a1: Tuple[Float64, Float64],
    a2: Tuple[Float64, Float64],
    b1: Tuple[Float64, Float64],
    b2: Tuple[Float64, Float64],
    eps: Float64 = 1e-12,
) -> List[Tuple[Float64, Float64, Float64, Float64, Int32]]:
    var r_x = a2[0] - a1[0]
    var r_y = a2[1] - a1[1]
    var s_x = b2[0] - b1[0]
    var s_y = b2[1] - b1[1]
    var cross_rs = r_x * s_y - r_y * s_x

    fn absf(x: Float64) -> Float64:
        if x < 0.0:
            return -x
        return x

    var out = List[Tuple[Float64, Float64, Float64, Float64, Int32]]()
    if absf(cross_rs) < eps:
        var q_p_x = b1[0] - a1[0]
        var q_p_y = b1[1] - a1[1]
        var cross_qp_r = q_p_x * r_y - q_p_y * r_x
        if absf(cross_qp_r) >= eps:
            return out.copy()
        var r2 = r_x * r_x + r_y * r_y
        if r2 == 0.0:
            return out.copy()
        var t0 = (q_p_x * r_x + q_p_y * r_y) / r2
        var t1 = ((b2[0] - a1[0]) * r_x + (b2[1] - a1[1]) * r_y) / r2
        var tmin = t0 if t0 < t1 else t1
        var tmax = t1 if t1 > t0 else t0
        var lo = 0.0 if tmin < 0.0 else tmin
        var hi = 1.0 if tmax > 1.0 else tmax
        if hi < 0.0 or lo > 1.0:
            return out.copy()
        var add0 = lo >= 0.0 and lo <= 1.0
        var add1 = hi >= 0.0 and hi <= 1.0 and absf(hi - lo) > eps
        if add0:
            var ix = a1[0] + lo * r_x
            var iy = a1[1] + lo * r_y
            var tb = 0.0
            var kind2: Int32 = 2
            out.append((ix, iy, lo, tb, kind2))
        if add1:
            var ix2 = a1[0] + hi * r_x
            var iy2 = a1[1] + hi * r_y
            var tb2 = 0.0
            var kind3: Int32 = 2
            out.append((ix2, iy2, hi, tb2, kind3))
        return out.copy()
    else:
        var q_p_x = b1[0] - a1[0]
        var q_p_y = b1[1] - a1[1]
        var t = (q_p_x * s_y - q_p_y * s_x) / cross_rs
        var u = (q_p_x * r_y - q_p_y * r_x) / cross_rs
        if t >= -eps and t <= 1.0 + eps and u >= -eps and u <= 1.0 + eps:
            var tt = 0.0 if t < 0.0 else (1.0 if t > 1.0 else t)
            var uu = 0.0 if u < 0.0 else (1.0 if u > 1.0 else u)
            var ix = a1[0] + tt * r_x
            var iy = a1[1] + tt * r_y
            var kind: Int32 = 0
            if tt <= eps or tt >= 1.0 - eps or uu <= eps or uu >= 1.0 - eps:
                kind = 1
            out.append((ix, iy, tt, uu, kind))
        return out.copy()
