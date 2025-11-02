from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, Polygon, MultiLineString, MultiPolygon
from shapely.algorithms import point_in_polygon, any_segment_intersection


fn length(g: Geometry) -> Float64:
    let t = g.__type_name__()
    if t == "LineString":
        let ls = unsafe_bitcast[LineString](g)
        return length(ls)
    if t == "MultiLineString":
        let mls = unsafe_bitcast[MultiLineString](g)
        return length(mls)
    return 0.0


fn area(g: Geometry) -> Float64:
    let t = g.__type_name__()
    if t == "Polygon":
        let p = unsafe_bitcast[Polygon](g)
        return area(p)
    if t == "MultiPolygon":
        let mp = unsafe_bitcast[MultiPolygon](g)
        return area(mp)
    return 0.0


fn distance(a: Geometry, b: Geometry) -> Float64:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Point" and tb == "Point":
        return distance(unsafe_bitcast[Point](a), unsafe_bitcast[Point](b))
    if ta == "Point" and tb == "LineString":
        return distance(unsafe_bitcast[Point](a), unsafe_bitcast[LineString](b))
    if ta == "LineString" and tb == "Point":
        return distance(unsafe_bitcast[LineString](a), unsafe_bitcast[Point](b))
    if ta == "Point" and tb == "Polygon":
        return distance(unsafe_bitcast[Point](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "Point":
        return distance(unsafe_bitcast[Polygon](a), unsafe_bitcast[Point](b))
    if ta == "LineString" and tb == "LineString":
        return distance(unsafe_bitcast[LineString](a), unsafe_bitcast[LineString](b))
    if ta == "LineString" and tb == "Polygon":
        return distance(unsafe_bitcast[LineString](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "LineString":
        return distance(unsafe_bitcast[Polygon](a), unsafe_bitcast[LineString](b))
    if ta == "Polygon" and tb == "Polygon":
        return distance(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    # MultiLineString folding
    if ta == "MultiLineString":
        let mls = unsafe_bitcast[MultiLineString](a)
        var best = 1.7976931348623157e308
        for ln in mls.lines:
            let d = distance(unsafe_bitcast[Geometry](ln), b)
            if d < best: best = d
        return best if best != 1.7976931348623157e308 else 0.0
    if tb == "MultiLineString":
        let mls2 = unsafe_bitcast[MultiLineString](b)
        var best2 = 1.7976931348623157e308
        for ln2 in mls2.lines:
            let d2 = distance(a, unsafe_bitcast[Geometry](ln2))
            if d2 < best2: best2 = d2
        return best2 if best2 != 1.7976931348623157e308 else 0.0
    # MultiPolygon folding
    if ta == "MultiPolygon":
        let mp = unsafe_bitcast[MultiPolygon](a)
        var best3 = 1.7976931348623157e308
        for p in mp.polys:
            let d3 = distance(unsafe_bitcast[Geometry](p), b)
            if d3 < best3: best3 = d3
        return best3 if best3 != 1.7976931348623157e308 else 0.0
    if tb == "MultiPolygon":
        let mp2 = unsafe_bitcast[MultiPolygon](b)
        var best4 = 1.7976931348623157e308
        for p2 in mp2.polys:
            let d4 = distance(a, unsafe_bitcast[Geometry](p2))
            if d4 < best4: best4 = d4
        return best4 if best4 != 1.7976931348623157e308 else 0.0
    return 0.0


fn length(line: LineString) -> Float64:
    if line.coords.size() <= 1:
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
    for i in range(0, line.coords.size() - 1):
        let a = line.coords[i]
        let b = line.coords[i + 1]
        let dx = b[0] - a[0]
        let dy = b[1] - a[1]
        total += sqrt_f64(dx * dx + dy * dy)
    return total


fn length(mls: MultiLineString) -> Float64:
    var s = 0.0
    for ln in mls.lines:
        s += length(ln)
    return s


fn area(poly: Polygon) -> Float64:
    # Shoelace over exterior minus holes
    let ring = poly.shell
    if ring.coords.size() < 3:
        return 0.0
    var shell_sum = 0.0
    for i in range(0, ring.coords.size() - 1):
        let a = ring.coords[i]
        let b = ring.coords[i + 1]
        shell_sum += a[0] * b[1] - a[1] * b[0]
    var holes_sum = 0.0
    for h in poly.holes:
        if h.coords.size() >= 3:
            var hs = 0.0
            for i in range(0, h.coords.size() - 1):
                let a = h.coords[i]
                let b = h.coords[i + 1]
                hs += a[0] * b[1] - a[1] * b[0]
            holes_sum += hs
    let total = 0.5 * (shell_sum - holes_sum)
    if total < 0.0: return -total
    return total


fn area(mpoly: MultiPolygon) -> Float64:
    var s = 0.0
    for p in mpoly.polys:
        s += area(p)
    return s


fn distance(a: Point, b: Point) -> Float64:
    let dx = a.x - b.x
    let dy = a.y - b.y
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
    if ls.coords.size() == 0:
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
    for i in range(0, ls.coords.size() - 1):
        let a1 = ls.coords[i]
        let a2 = ls.coords[i + 1]
        let vx = a2[0] - a1[0]
        let vy = a2[1] - a1[1]
        let wx = a.x - a1[0]
        let wy = a.y - a1[1]
        let vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = (wx * vx + wy * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        let px = a1[0] + t * vx
        let py = a1[1] + t * vy
        let dx = a.x - px
        let dy = a.y - py
        let d = dx * dx + dy * dy
        if d < best: best = d
    return sqrt_f64(best)


fn distance(a: LineString, b: LineString) -> Float64:
    if a.coords.size() < 2 or b.coords.size() < 2:
        # fallback to endpoint distances
        if a.coords.size() == 0 or b.coords.size() == 0:
            return 0.0
        let pa = Point(a.coords[0][0], a.coords[0][1])
        let pb = Point(b.coords[0][0], b.coords[0][1])
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
        let vx = bx - ax
        let vy = by - ay
        let vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = ((px - ax) * vx + (py - ay) * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        let cx = ax + t * vx
        let cy = ay + t * vy
        let dx = px - cx
        let dy = py - cy
        return dx * dx + dy * dy
    var best = 1.7976931348623157e308
    for i in range(0, a.coords.size() - 1):
        let a1 = a.coords[i]
        let a2 = a.coords[i + 1]
        for j in range(0, b.coords.size() - 1):
            let b1 = b.coords[j]
            let b2 = b.coords[j + 1]
            let d2 = pt_seg_d2(a1[0], a1[1], b1[0], b1[1], b2[0], b2[1])
            if d2 < best: best = d2
            let d2b = pt_seg_d2(a2[0], a2[1], b1[0], b1[1], b2[0], b2[1])
            if d2b < best: best = d2b
            let d2c = pt_seg_d2(b1[0], b1[1], a1[0], a1[1], a2[0], a2[1])
            if d2c < best: best = d2c
            let d2d = pt_seg_d2(b2[0], b2[1], a1[0], a1[1], a2[0], a2[1])
            if d2d < best: best = d2d
    return sqrt_f64(best)


fn distance(ls: LineString, poly: Polygon) -> Float64:
    # If intersects or endpoint inside, distance is 0
    # Check shell/hole intersections by segments
    let shell_ls = LineString(poly.shell.coords)
    if any_segment_intersection(ls, shell_ls):
        return 0.0
    for h in poly.holes:
        let hls = LineString(h.coords)
        if any_segment_intersection(ls, hls):
            return 0.0
    if ls.coords.size() > 0:
        let p0 = Point(ls.coords[0][0], ls.coords[0][1])
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
        let vx = bx - ax
        let vy = by - ay
        let vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = ((px - ax) * vx + (py - ay) * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        let cx = ax + t * vx
        let cy = ay + t * vy
        let dx = px - cx
        let dy = py - cy
        return dx * dx + dy * dy
    var best = 1.7976931348623157e308
    # ls segments vs shell segments
    for i in range(0, ls.coords.size() - 1):
        let p1 = ls.coords[i]
        let p2 = ls.coords[i + 1]
        for j in range(0, poly.shell.coords.size() - 1):
            let a = poly.shell.coords[j]
            let b = poly.shell.coords[j + 1]
            let d2a = pt_seg_d2(p1[0], p1[1], a[0], a[1], b[0], b[1])
            if d2a < best: best = d2a
            let d2b = pt_seg_d2(p2[0], p2[1], a[0], a[1], b[0], b[1])
            if d2b < best: best = d2b
            let d2c = pt_seg_d2(a[0], a[1], p1[0], p1[1], p2[0], p2[1])
            if d2c < best: best = d2c
            let d2d = pt_seg_d2(b[0], b[1], p1[0], p1[1], p2[0], p2[1])
            if d2d < best: best = d2d
    # holes
    for h in poly.holes:
        for i in range(0, ls.coords.size() - 1):
            let p1 = ls.coords[i]
            let p2 = ls.coords[i + 1]
            for j in range(0, h.coords.size() - 1):
                let a = h.coords[j]
                let b = h.coords[j + 1]
                let d2a = pt_seg_d2(p1[0], p1[1], a[0], a[1], b[0], b[1])
                if d2a < best: best = d2a
                let d2b = pt_seg_d2(p2[0], p2[1], a[0], a[1], b[0], b[1])
                if d2b < best: best = d2b
                let d2c = pt_seg_d2(a[0], a[1], p1[0], p1[1], p2[0], p2[1])
                if d2c < best: best = d2c
                let d2d = pt_seg_d2(b[0], b[1], p1[0], p1[1], p2[0], p2[1])
                if d2d < best: best = d2d
    return sqrt_f64(best)


fn distance(poly: Polygon, ls: LineString) -> Float64:
    return distance(ls, poly)


fn distance(a: Polygon, b: Polygon) -> Float64:
    # zero if they intersect or one contains a vertex of the other
    let a_ls = LineString(a.shell.coords)
    let b_ls = LineString(b.shell.coords)
    if any_segment_intersection(a_ls, b_ls):
        return 0.0
    if a.shell.coords.size() > 0:
        let p = Point(a.shell.coords[0][0], a.shell.coords[0][1])
        if point_in_polygon(p, b) != 0: return 0.0
    if b.shell.coords.size() > 0:
        let q = Point(b.shell.coords[0][0], b.shell.coords[0][1])
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
        let vx = bx - ax
        let vy = by - ay
        let vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = ((px - ax) * vx + (py - ay) * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        let cx = ax + t * vx
        let cy = ay + t * vy
        let dx = px - cx
        let dy = py - cy
        return dx * dx + dy * dy
    var best = 1.7976931348623157e308
    for i in range(0, a.shell.coords.size() - 1):
        let a1 = a.shell.coords[i]
        let a2 = a.shell.coords[i + 1]
        for j in range(0, b.shell.coords.size() - 1):
            let b1 = b.shell.coords[j]
            let b2 = b.shell.coords[j + 1]
            let d2a = pt_seg_d2(a1[0], a1[1], b1[0], b1[1], b2[0], b2[1])
            if d2a < best: best = d2a
            let d2b = pt_seg_d2(a2[0], a2[1], b1[0], b1[1], b2[0], b2[1])
            if d2b < best: best = d2b
            let d2c = pt_seg_d2(b1[0], b1[1], a1[0], a1[1], a2[0], a2[1])
            if d2c < best: best = d2c
            let d2d = pt_seg_d2(b2[0], b2[1], a1[0], a1[1], a2[0], a2[1])
            if d2d < best: best = d2d
    return sqrt_f64(best)


fn distance(ls: LineString, p: Point) -> Float64:
    return distance(p, ls)


fn distance(p: Point, poly: Polygon) -> Float64:
    # inside or on boundary -> 0
    let rel = point_in_polygon(p, poly)
    if rel != 0:
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
    # shell
    let ring = poly.shell
    for i in range(0, ring.coords.size() - 1):
        let a1 = ring.coords[i]
        let a2 = ring.coords[i + 1]
        let vx = a2[0] - a1[0]
        let vy = a2[1] - a1[1]
        let wx = p.x - a1[0]
        let wy = p.y - a1[1]
        let vlen2 = vx * vx + vy * vy
        var t = 0.0
        if vlen2 > 0.0:
            t = (wx * vx + wy * vy) / vlen2
        if t < 0.0: t = 0.0
        if t > 1.0: t = 1.0
        let px = a1[0] + t * vx
        let py = a1[1] + t * vy
        let dx = p.x - px
        let dy = p.y - py
        let d = dx * dx + dy * dy
        if d < best: best = d
    # holes
    for h in poly.holes:
        for i in range(0, h.coords.size() - 1):
            let a1 = h.coords[i]
            let a2 = h.coords[i + 1]
            let vx = a2[0] - a1[0]
            let vy = a2[1] - a1[1]
            let wx = p.x - a1[0]
            let wy = p.y - a1[1]
            let vlen2 = vx * vx + vy * vy
            var t = 0.0
            if vlen2 > 0.0:
                t = (wx * vx + wy * vy) / vlen2
            if t < 0.0: t = 0.0
            if t > 1.0: t = 1.0
            let px = a1[0] + t * vx
            let py = a1[1] + t * vy
            let dx = p.x - px
            let dy = p.y - py
            let d = dx * dx + dy * dy
            if d < best: best = d
    return sqrt_f64(best)
