from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, Polygon, MultiLineString, MultiPolygon
from shapely.algorithms import point_in_polygon, any_segment_intersection


fn distance(ref a: Point, ref b: Geometry) -> Float64:
    return distance(Geometry(a.copy()), b)


fn distance(ref a: Geometry, ref b: Point) -> Float64:
    return distance(a, Geometry(b.copy()))


fn distance(ref a: LineString, ref b: Geometry) -> Float64:
    return distance(Geometry(a.copy()), b)


fn distance(ref a: Geometry, ref b: LineString) -> Float64:
    return distance(a, Geometry(b.copy()))


fn distance(ref a: Polygon, ref b: Geometry) -> Float64:
    return distance(Geometry(a.copy()), b)


fn distance(ref a: Geometry, ref b: Polygon) -> Float64:
    return distance(a, Geometry(b.copy()))


fn distance(ref a: Polygon, ref b: Point) -> Float64:
    return distance(b.copy(), a.copy())


fn distance(ref a: Point, ref b: Polygon) -> Float64:
    return distance(a.copy(), b.copy())


fn length(g: Geometry) -> Float64:
    if g.is_linestring():
        return length(g.as_linestring())
    if g.is_multilinestring():
        return length(g.as_multilinestring())
    return 0.0


fn area(g: Geometry) -> Float64:
    if g.is_polygon():
        return area(g.as_polygon())
    if g.is_multipolygon():
        return area(g.as_multipolygon())
    return 0.0


fn distance(a: Geometry, b: Geometry) -> Float64:
    if a.is_point() and b.is_point():
        return distance(a.as_point(), b.as_point())
    if a.is_point() and b.is_linestring():
        return distance(a.as_point(), b.as_linestring())
    if a.is_linestring() and b.is_point():
        return distance(a.as_linestring(), b.as_point())
    if a.is_point() and b.is_polygon():
        return distance(a.as_point(), b.as_polygon())
    if a.is_polygon() and b.is_point():
        return distance(a.as_polygon(), b.as_point())
    if a.is_linestring() and b.is_linestring():
        return distance(a.as_linestring(), b.as_linestring())
    if a.is_linestring() and b.is_polygon():
        return distance(a.as_linestring(), b.as_polygon())
    if a.is_polygon() and b.is_linestring():
        return distance(a.as_polygon(), b.as_linestring())
    if a.is_polygon() and b.is_polygon():
        return distance(a.as_polygon(), b.as_polygon())
    # MultiLineString folding
    if a.is_multilinestring():
        var mls = a.as_multilinestring()
        var best = 1.7976931348623157e308
        for ln in mls.lines:
            var d = distance(Geometry(ln.copy()), b)
            if d < best:
                best = d
        return best if best != 1.7976931348623157e308 else 0.0
    if b.is_multilinestring():
        var mls2 = b.as_multilinestring()
        var best2 = 1.7976931348623157e308
        for ln2 in mls2.lines:
            var d2 = distance(a, Geometry(ln2.copy()))
            if d2 < best2:
                best2 = d2
        return best2 if best2 != 1.7976931348623157e308 else 0.0
    # MultiPolygon folding
    if a.is_multipolygon():
        var mp = a.as_multipolygon()
        var best3 = 1.7976931348623157e308
        for p in mp.polys:
            var d3 = distance(Geometry(p.copy()), b)
            if d3 < best3:
                best3 = d3
        return best3 if best3 != 1.7976931348623157e308 else 0.0
    if b.is_multipolygon():
        var mp2 = b.as_multipolygon()
        var best4 = 1.7976931348623157e308
        for p2 in mp2.polys:
            var d4 = distance(a, Geometry(p2.copy()))
            if d4 < best4:
                best4 = d4
        return best4 if best4 != 1.7976931348623157e308 else 0.0
    return 0.0


fn length(line: LineString) -> Float64:
    if line.coords.__len__() <= 1:
        return 0.0
    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r
    var total = 0.0
    for i in range(0, line.coords.__len__() - 1):
        var a = line.coords[i]
        var b = line.coords[i + 1]
        var dx = b[0] - a[0]
        var dy = b[1] - a[1]
        total += sqrt_f64(dx * dx + dy * dy)
    return total


fn length(mls: MultiLineString) -> Float64:
    var s = 0.0
    for ln in mls.lines:
        s += length(ln)
    return s


fn area(poly: Polygon) -> Float64:
    # Shoelace over exterior minus holes
    ref ring = poly.shell
    if ring.coords.__len__() < 3:
        return 0.0
    var shell_sum = 0.0
    for i in range(0, ring.coords.__len__() - 1):
        ref a = ring.coords[i]
        ref b = ring.coords[i + 1]
        shell_sum += a[0] * b[1] - a[1] * b[0]
    var holes_sum = 0.0
    for h in poly.holes:
        if h.coords.__len__() >= 3:
            var hs = 0.0
            for i in range(0, h.coords.__len__() - 1):
                ref a = h.coords[i]
                ref b = h.coords[i + 1]
                hs += a[0] * b[1] - a[1] * b[0]
            holes_sum += 0.5 * abs(hs)
    var total = 0.5 * (shell_sum - holes_sum)
    if total < 0.0: return -total
    return total


fn area(mpoly: MultiPolygon) -> Float64:
    var s = 0.0
    for p in mpoly.polys:
        s += area(p)
    return s


fn distance(a: Point, b: Point) -> Float64:
    var dx = a.x - b.x
    var dy = a.y - b.y
    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r
    return sqrt_f64(dx * dx + dy * dy)


fn distance(a: Point, ls: LineString) -> Float64:
    if ls.coords.__len__() == 0:
        return 0.0
    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r
    var best = 1.7976931348623157e308
    for i in range(0, ls.coords.__len__() - 1):
        var a1 = ls.coords[i]
        var a2 = ls.coords[i + 1]
        var vx = a2[0] - a1[0]
        var vy = a2[1] - a1[1]
        var wx = a.x - a1[0]
        var wy = a.y - a1[1]
        var vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = (wx * vx + wy * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        var px = a1[0] + t * vx
        var py = a1[1] + t * vy
        var dx = a.x - px
        var dy = a.y - py
        var d = dx * dx + dy * dy
        if d < best: best = d
    return sqrt_f64(best)


fn distance(a: LineString, b: LineString) -> Float64:
    if a.coords.__len__() < 2 or b.coords.__len__() < 2:
        # fallback to endpoint distances
        if a.coords.__len__() == 0 or b.coords.__len__() == 0:
            return 0.0
        var pa = Point(a.coords[0][0], a.coords[0][1])
        var pb = Point(b.coords[0][0], b.coords[0][1])
        return distance(pa, pb)
    if any_segment_intersection(a, b):
        return 0.0
    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r
    fn pt_seg_d2(px: Float64, py: Float64, ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
        var vx = bx - ax
        var vy = by - ay
        var vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = ((px - ax) * vx + (py - ay) * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        var cx = ax + t * vx
        var cy = ay + t * vy
        var dx = px - cx
        var dy = py - cy
        return dx * dx + dy * dy
    var best = 1.7976931348623157e308
    for i in range(0, a.coords.__len__() - 1):
        var a1 = a.coords[i]
        var a2 = a.coords[i + 1]
        for j in range(0, b.coords.__len__() - 1):
            var b1 = b.coords[j]
            var b2 = b.coords[j + 1]
            var d2 = pt_seg_d2(a1[0], a1[1], b1[0], b1[1], b2[0], b2[1])
            if d2 < best: best = d2
            var d2b = pt_seg_d2(a2[0], a2[1], b1[0], b1[1], b2[0], b2[1])
            if d2b < best: best = d2b
            var d2c = pt_seg_d2(b1[0], b1[1], a1[0], a1[1], a2[0], a2[1])
            if d2c < best: best = d2c
            var d2d = pt_seg_d2(b2[0], b2[1], a1[0], a1[1], a2[0], a2[1])
            if d2d < best: best = d2d
    return sqrt_f64(best)


fn distance(ls: LineString, poly: Polygon) -> Float64:
    # If intersects or endpoint inside, distance is 0
    # Check shell/hole intersections by segments
    var shell_ls = LineString(poly.shell.coords)
    if any_segment_intersection(ls, shell_ls):
        return 0.0
    for h in poly.holes:
        var hls = LineString(h.coords)
        if any_segment_intersection(ls, hls):
            return 0.0
    if ls.coords.__len__() > 0:
        var p0 = Point(ls.coords[0][0], ls.coords[0][1])
        if point_in_polygon(p0, poly) != 0:
            return 0.0
    # otherwise compute min over segments to polygon edges
    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r
    fn pt_seg_d2(px: Float64, py: Float64, ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
        var vx = bx - ax
        var vy = by - ay
        var vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = ((px - ax) * vx + (py - ay) * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        var cx = ax + t * vx
        var cy = ay + t * vy
        var dx = px - cx
        var dy = py - cy
        return dx * dx + dy * dy
    var best = 1.7976931348623157e308
    # ls segments vs shell segments
    for i in range(0, ls.coords.__len__() - 1):
        var p1 = ls.coords[i]
        var p2 = ls.coords[i + 1]
        for j in range(0, poly.shell.coords.__len__() - 1):
            var a = poly.shell.coords[j]
            var b = poly.shell.coords[j + 1]
            var d2a = pt_seg_d2(p1[0], p1[1], a[0], a[1], b[0], b[1])
            if d2a < best: best = d2a
            var d2b = pt_seg_d2(p2[0], p2[1], a[0], a[1], b[0], b[1])
            if d2b < best: best = d2b
            var d2c = pt_seg_d2(a[0], a[1], p1[0], p1[1], p2[0], p2[1])
            if d2c < best: best = d2c
            var d2d = pt_seg_d2(b[0], b[1], p1[0], p1[1], p2[0], p2[1])
            if d2d < best: best = d2d
    # holes
    for h in poly.holes:
        for i in range(0, ls.coords.__len__() - 1):
            var p1 = ls.coords[i]
            var p2 = ls.coords[i + 1]
            for j in range(0, h.coords.__len__() - 1):
                var a = h.coords[j]
                var b = h.coords[j + 1]
                var d2a = pt_seg_d2(p1[0], p1[1], a[0], a[1], b[0], b[1])
                if d2a < best: best = d2a
                var d2b = pt_seg_d2(p2[0], p2[1], a[0], a[1], b[0], b[1])
                if d2b < best: best = d2b
                var d2c = pt_seg_d2(a[0], a[1], p1[0], p1[1], p2[0], p2[1])
                if d2c < best: best = d2c
                var d2d = pt_seg_d2(b[0], b[1], p1[0], p1[1], p2[0], p2[1])
                if d2d < best: best = d2d
    return sqrt_f64(best)


fn distance(poly: Polygon, ls: LineString) -> Float64:
    return distance(ls, poly)


fn distance(a: Polygon, b: Polygon) -> Float64:
    # zero if they intersect or one contains a vertex of the other
    var a_ls = LineString(a.shell.coords)
    var b_ls = LineString(b.shell.coords)
    if any_segment_intersection(a_ls, b_ls):
        return 0.0
    if a.shell.coords.__len__() > 0:
        var p = Point(a.shell.coords[0][0], a.shell.coords[0][1])
        if point_in_polygon(p, b) != 0: return 0.0
    if b.shell.coords.__len__() > 0:
        var q = Point(b.shell.coords[0][0], b.shell.coords[0][1])
        if point_in_polygon(q, a) != 0: return 0.0
    # otherwise min distance between shell segments
    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r
    fn pt_seg_d2(px: Float64, py: Float64, ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
        var vx = bx - ax
        var vy = by - ay
        var vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = ((px - ax) * vx + (py - ay) * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        var cx = ax + t * vx
        var cy = ay + t * vy
        var dx = px - cx
        var dy = py - cy
        return dx * dx + dy * dy
    var best = 1.7976931348623157e308
    for i in range(0, a.shell.coords.__len__() - 1):
        var a1 = a.shell.coords[i]
        var a2 = a.shell.coords[i + 1]
        for j in range(0, b.shell.coords.__len__() - 1):
            var b1 = b.shell.coords[j]
            var b2 = b.shell.coords[j + 1]
            var d2a = pt_seg_d2(a1[0], a1[1], b1[0], b1[1], b2[0], b2[1])
            if d2a < best: best = d2a
            var d2b = pt_seg_d2(a2[0], a2[1], b1[0], b1[1], b2[0], b2[1])
            if d2b < best: best = d2b
            var d2c = pt_seg_d2(b1[0], b1[1], a1[0], a1[1], a2[0], a2[1])
            if d2c < best: best = d2c
            var d2d = pt_seg_d2(b2[0], b2[1], a1[0], a1[1], a2[0], a2[1])
            if d2d < best: best = d2d
    return sqrt_f64(best)


fn distance(ls: LineString, p: Point) -> Float64:
    return distance(p, ls)


fn distance(p: Point, poly: Polygon) -> Float64:
    # inside or on boundary -> 0
    var rel = point_in_polygon(p, poly)
    if rel != 0:
        return 0.0

    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0:
            return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r

    var best = 1.7976931348623157e308

    # shell
    ref ring = poly.shell
    for i in range(0, ring.coords.__len__() - 1):
        var a1 = ring.coords[i]
        var a2 = ring.coords[i + 1]
        var vx = a2[0] - a1[0]
        var vy = a2[1] - a1[1]
        var wx = p.x - a1[0]
        var wy = p.y - a1[1]
        var vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = (wx * vx + wy * vy) / vlen2
        if t < 0.0:
            t = 0.0
        if t > 1.0:
            t = 1.0
        var px = a1[0] + t * vx
        var py = a1[1] + t * vy
        var dx = p.x - px
        var dy = p.y - py
        var d = dx * dx + dy * dy
        if d < best:
            best = d

    # holes
    for h in poly.holes:
        for i in range(0, h.coords.__len__() - 1):
            var a1 = h.coords[i]
            var a2 = h.coords[i + 1]
            var vx = a2[0] - a1[0]
            var vy = a2[1] - a1[1]
            var wx = p.x - a1[0]
            var wy = p.y - a1[1]
            var vlen2 = vx * vx + vy * vy
            var t = 0.0
            if vlen2 > 0.0:
                t = (wx * vx + wy * vy) / vlen2
            if t < 0.0:
                t = 0.0
            if t > 1.0:
                t = 1.0
            var px = a1[0] + t * vx
            var py = a1[1] + t * vy
            var dx = p.x - px
            var dy = p.y - py
            var d = dx * dx + dy * dy
            if d < best:
                best = d

    return sqrt_f64(best)
