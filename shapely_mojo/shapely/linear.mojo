from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, MultiLineString, GeometryCollection


fn line_interpolate_point(line: LineString, distance: Float64, normalized: Bool = False) -> Point:
    if line.coords.size() == 0:
        return Point(0.0, 0.0)
    if line.coords.size() == 1:
        return Point(line.coords[0][0], line.coords[0][1])

    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r

    # total length
    var total = 0.0
    for i in range(0, line.coords.size() - 1):
        var a = line.coords[i]
        var b = line.coords[i + 1]
        var dx = b[0] - a[0]
        var dy = b[1] - a[1]
        total += sqrt_f64(dx * dx + dy * dy)
    if total == 0.0:
        return Point(line.coords[0][0], line.coords[0][1])

    var d = distance
    if normalized:
        if d <= 0.0: return Point(line.coords[0][0], line.coords[0][1])
        if d >= 1.0: return Point(line.coords[line.coords.size() - 1][0], line.coords[line.coords.size() - 1][1])
        d = d * total
    else:
        if d < 0.0: d = total + d
        if d <= 0.0: return Point(line.coords[0][0], line.coords[0][1])
        if d >= total: return Point(line.coords[line.coords.size() - 1][0], line.coords[line.coords.size() - 1][1])

    # walk segments
    var acc = 0.0
    for i in range(0, line.coords.size() - 1):
        var a = line.coords[i]
        var b = line.coords[i + 1]
        var dx = b[0] - a[0]
        var dy = b[1] - a[1]
        var seg = sqrt_f64(dx * dx + dy * dy)
        if acc + seg >= d:
            var t = (d - acc) / seg
            return Point(a[0] + t * dx, a[1] + t * dy)
        acc += seg
    return Point(line.coords[line.coords.size() - 1][0], line.coords[line.coords.size() - 1][1])


fn line_locate_point(_line: LineString, _other: Point, normalized: Bool = False) -> Float64:
    if _line.coords.size() == 0:
        return 0.0
    if _line.coords.size() == 1:
        return 0.0

    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r

    # total length
    var total = 0.0
    for i in range(0, _line.coords.size() - 1):
        var a = _line.coords[i]
        var b = _line.coords[i + 1]
        var dx = b[0] - a[0]
        var dy = b[1] - a[1]
        total += sqrt_f64(dx * dx + dy * dy)
    if total == 0.0:
        return 0.0

    # walk segments
    var acc = 0.0
    for i in range(0, _line.coords.size() - 1):
        var a = _line.coords[i]
        var b = _line.coords[i + 1]
        var dx = b[0] - a[0]
        var dy = b[1] - a[1]
        var seg = sqrt_f64(dx * dx + dy * dy)
        var ap = Point(a[0], a[1])
        var bp = Point(b[0], b[1])
        var op = _other
        var ax = ap.x
        var ay = ap.y
        var bx = bp.x
        var by = bp.y
        var ox = op.x
        var oy = op.y
        var dot = ((ox - ax) * (bx - ax) + (oy - ay) * (by - ay))
        if dot <= 0.0:
            return acc
        var sq = (bx - ax) * (bx - ax) + (by - ay) * (by - ay)
        if dot >= sq:
            return acc + seg
        var proj = dot / sq
        var px = ax + proj * (bx - ax)
        var py = ay + proj * (by - ay)
        var pdx = ox - px
        var pdy = oy - py
        var pd = sqrt_f64(pdx * pdx + pdy * pdy)
        if pd <= 1e-9:
            return acc + proj * seg
        acc += seg
    return total


fn line_merge(line: LineString) -> LineString:
    # passthrough
    return line


fn line_merge(lines: MultiLineString) -> Geometry:
    # Greedy endpoint-to-endpoint merge of parts
    var parts = List[List[Tuple[Float64, Float64]]]()
    for ln in lines.lines:
        var seq = List[Tuple[Float64, Float64]]()
        for c in ln.coords: seq.append(c)
        parts.append(seq)

    fn first(ps: List[Tuple[Float64, Float64]]) -> Tuple[Float64, Float64]:
        return ps[0]
    fn last(ps: List[Tuple[Float64, Float64]]) -> Tuple[Float64, Float64]:
        return ps[ps.size() - 1]
    fn reverse_list(xs: List[Tuple[Float64, Float64]]) -> List[Tuple[Float64, Float64]]:
        var out = List[Tuple[Float64, Float64]]()
        var i = xs.size() - 1
        while True:
            out.append(xs[i])
            if i == 0: break
            i -= 1
        return out

    var changed = True
    while changed:
        changed = False
        var i = 0
        while i < parts.size():
            var j = i + 1
            while j < parts.size():
                var ai0 = first(parts[i])
                var ai1 = last(parts[i])
                var bj0 = first(parts[j])
                var bj1 = last(parts[j])
                if ai1[0] == bj0[0] and ai1[1] == bj0[1]:
                    # i ... + j ...
                    for k in range(1, parts[j].size()):
                        parts[i].append(parts[j][k])
                    parts.erase(j)
                    changed = True
                    continue
                if ai1[0] == bj1[0] and ai1[1] == bj1[1]:
                    var rev = reverse_list(parts[j])
                    for k in range(1, rev.size()):
                        parts[i].append(rev[k])
                    parts.erase(j)
                    changed = True
                    continue
                if ai0[0] == bj1[0] and ai0[1] == bj1[1]:
                    # j ... + i ...
                    var newp = parts[j]
                    for k in range(1, parts[i].size()):
                        newp.append(parts[i][k])
                    parts[i] = newp
                    parts.erase(j)
                    changed = True
                    continue
                if ai0[0] == bj0[0] and ai0[1] == bj0[1]:
                    var rev2 = reverse_list(parts[j])
                    var newp2 = rev2
                    for k in range(1, parts[i].size()):
                        newp2.append(parts[i][k])
                    parts[i] = newp2
                    parts.erase(j)
                    changed = True
                    continue
                j += 1
            i += 1

    if parts.size() == 0:
        return LineString([])
    if parts.size() == 1:
        return LineString(parts[0])
    var out = List[LineString]()
    for p in parts:
        out.append(LineString(p))
    return MultiLineString(out)


fn shared_paths(a: LineString, b: LineString) -> GeometryCollection:
    var forward = List[LineString]()
    var reverse = List[LineString]()
    for i in range(0, a.coords.size() - 1):
        var a1 = a.coords[i]
        var a2 = a.coords[i + 1]
        for j in range(0, b.coords.size() - 1):
            var b1 = b.coords[j]
            var b2 = b.coords[j + 1]
            if a1[0] == b1[0] and a1[1] == b1[1] and a2[0] == b2[0] and a2[1] == b2[1]:
                forward.append(LineString([a1, a2]))
            elif a1[0] == b2[0] and a1[1] == b2[1] and a2[0] == b1[0] and a2[1] == b1[1]:
                reverse.append(LineString([a1, a2]))
    return GeometryCollection([MultiLineString(forward), MultiLineString(reverse)])


fn shortest_line(_a: Geometry, _b: Geometry) -> LineString:
    # placeholder: origin-to-origin
    return LineString([(0.0, 0.0), (0.0, 0.0)])
