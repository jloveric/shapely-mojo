from shapely._geometry import Geometry
from shapely.geometry import GeometryCollection, LinearRing, Polygon, MultiPolygon, LineString, Point
from shapely.algorithms import orientation, point_in_polygon
from shapely.overlay import overlay_union, overlay_difference, overlay_intersection, overlay_xor


fn _empty_polygon() -> Polygon:
    return Polygon(LinearRing(List[Tuple[Float64, Float64]]()))


fn _is_empty_polygon(p: Polygon) -> Bool:
    return p.shell.coords.__len__() < 3


fn union(a: Geometry, b: Geometry) -> Geometry:
    if a.is_polygon() and b.is_polygon():
        return union(a.as_polygon(), b.as_polygon())
    if a.is_multipolygon() and b.is_polygon():
        return union(a.as_multipolygon(), b.as_polygon())
    if a.is_polygon() and b.is_multipolygon():
        return union(a.as_polygon(), b.as_multipolygon())
    if a.is_multipolygon() and b.is_multipolygon():
        return union(a.as_multipolygon(), b.as_multipolygon())
    return Geometry(_empty_polygon())

fn union(a: Geometry, b: Polygon) -> Geometry:
    if a.is_polygon():
        return union(a.as_polygon(), b)
    if a.is_multipolygon():
        return union(a.as_multipolygon(), b)
    return Geometry(_empty_polygon())

fn union(a: Polygon, b: Geometry) -> Geometry:
    if b.is_polygon():
        return union(a, b.as_polygon())
    if b.is_multipolygon():
        return union(a, b.as_multipolygon())
    return Geometry(_empty_polygon())


fn union(a: Polygon, b: Polygon) -> Geometry:
    return overlay_union(a, b)


fn union(mp: MultiPolygon, p: Polygon) -> Geometry:
    if mp.polys.__len__() == 0:
        return Geometry(p.copy())
    var g: Geometry = overlay_union(mp.polys[0], p)
    var i = 1
    while i < mp.polys.__len__():
        g = union(g, mp.polys[i])
        i += 1
    return g.copy()


fn union(p: Polygon, mp: MultiPolygon) -> Geometry:
    return union(mp, p)


fn union(a: MultiPolygon, b: MultiPolygon) -> Geometry:
    if a.polys.__len__() == 0:
        if b.polys.__len__() == 0:
            return Geometry(MultiPolygon([]))
        var g0: Geometry = Geometry(b.polys[0].copy())
        var jj = 1
        while jj < b.polys.__len__():
            g0 = union(g0, b.polys[jj])
            jj += 1
        return g0.copy()
    var g: Geometry = Geometry(a.polys[0].copy())
    var i = 1
    while i < a.polys.__len__():
        g = union(g, a.polys[i])
        i += 1
    var j = 0
    while j < b.polys.__len__():
        g = union(g, b.polys[j])
        j += 1
    return g.copy()


fn _cross(ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
    return ax * by - ay * bx


fn _intersect_point(s1: Tuple[Float64, Float64], s2: Tuple[Float64, Float64], c1: Tuple[Float64, Float64], c2: Tuple[Float64, Float64]) -> (Tuple[Float64, Float64], Bool):
    var r_x = s2[0] - s1[0]
    var r_y = s2[1] - s1[1]
    var s_x = c2[0] - c1[0]
    var s_y = c2[1] - c1[1]
    var denom = _cross(r_x, r_y, s_x, s_y)
    if denom == 0.0:
        return ((0.0, 0.0), False)
    var t = _cross(c1[0] - s1[0], c1[1] - s1[1], s_x, s_y) / denom
    return ((s1[0] + t * r_x, s1[1] + t * r_y), True)


fn _is_inside(p: Tuple[Float64, Float64], a: Tuple[Float64, Float64], b: Tuple[Float64, Float64]) -> Bool:
    # left-of test: inside if p is to left of a->b
    return orientation(a[0], a[1], b[0], b[1], p[0], p[1]) >= 0


fn _suth_hodg(subject: List[Tuple[Float64, Float64]], clip: List[Tuple[Float64, Float64]]) -> List[Tuple[Float64, Float64]]:
    var output = subject
    if output.__len__() == 0:
        return output
    var i = 0
    while i < clip.__len__() - 1:
        var A = clip[i]
        var B = clip[i + 1]
        var input = output
        output = List[Tuple[Float64, Float64]]()
        if input.__len__() == 0:
            break
        var S = input[input.__len__() - 1]
        for E in input:
            if _is_inside(E, A, B):
                if not _is_inside(S, A, B):
                    var rv = _intersect_point(S, E, A, B)
                    if rv[1]: output.append(rv[0])
                output.append(E)
            else:
                if _is_inside(S, A, B):
                    var rv2 = _intersect_point(S, E, A, B)
                    if rv2[1]: output.append(rv2[0])
            S = E
        i += 1
    return output


fn intersection(a: Polygon, b: Polygon) -> Geometry:
    return overlay_intersection(a, b)


fn intersection(a: Geometry, b: Geometry) -> Geometry:
    if a.is_polygon() and b.is_polygon():
        return intersection(a.as_polygon(), b.as_polygon())
    if a.is_multipolygon() and b.is_polygon():
        return intersection(a.as_multipolygon().polys[0], b) if a.as_multipolygon().polys.__len__() > 0 else Geometry(_empty_polygon())
    if a.is_polygon() and b.is_multipolygon():
        return intersection(a, b.as_multipolygon().polys[0]) if b.as_multipolygon().polys.__len__() > 0 else Geometry(_empty_polygon())
    return Geometry(_empty_polygon())

fn intersection(a: Geometry, b: Polygon) -> Geometry:
    if a.is_polygon():
        return intersection(a.as_polygon(), b)
    return Geometry(_empty_polygon())

fn intersection(a: Polygon, b: Geometry) -> Geometry:
    if b.is_polygon():
        return intersection(a, b.as_polygon())
    return Geometry(_empty_polygon())


fn difference(a: Polygon, b: Polygon) -> Geometry:
    return overlay_difference(a, b)


fn difference(a: Geometry, b: Geometry) -> Geometry:
    if a.is_polygon() and b.is_polygon():
        return difference(a.as_polygon(), b.as_polygon())
    if a.is_multipolygon() and b.is_polygon():
        return difference(a.as_multipolygon(), b.as_polygon())
    if a.is_polygon() and b.is_multipolygon():
        return difference(a.as_polygon(), b.as_multipolygon())
    if a.is_multipolygon() and b.is_multipolygon():
        return difference(a.as_multipolygon(), b.as_multipolygon())
    return Geometry(_empty_polygon())

fn difference(a: Geometry, b: Polygon) -> Geometry:
    if a.is_polygon():
        return difference(a.as_polygon(), b)
    if a.is_multipolygon():
        return difference(a.as_multipolygon(), b)
    return Geometry(_empty_polygon())

fn difference(a: Polygon, b: Geometry) -> Geometry:
    if b.is_polygon():
        return difference(a, b.as_polygon())
    if b.is_multipolygon():
        return difference(a, b.as_multipolygon())
    return Geometry(_empty_polygon())


fn symmetric_difference(a: Geometry, b: Geometry) -> Geometry:
    if a.is_polygon() and b.is_polygon():
        return symmetric_difference(a.as_polygon(), b.as_polygon())
    if a.is_multipolygon() and b.is_polygon():
        return symmetric_difference(a.as_multipolygon(), b.as_polygon())
    if a.is_polygon() and b.is_multipolygon():
        return symmetric_difference(a.as_polygon(), b.as_multipolygon())
    if a.is_multipolygon() and b.is_multipolygon():
        return symmetric_difference(a.as_multipolygon(), b.as_multipolygon())
    return Geometry(_empty_polygon())

fn symmetric_difference(a: Geometry, b: Polygon) -> Geometry:
    if a.is_polygon():
        return symmetric_difference(a.as_polygon(), b)
    if a.is_multipolygon():
        return symmetric_difference(a.as_multipolygon(), b)
    return Geometry(_empty_polygon())

fn symmetric_difference(a: Polygon, b: Geometry) -> Geometry:
    if b.is_polygon():
        return symmetric_difference(a, b.as_polygon())
    if b.is_multipolygon():
        return symmetric_difference(a, b.as_multipolygon())
    return Geometry(_empty_polygon())


fn unary_union(geoms: List[Geometry]) -> Geometry:
    if geoms.__len__() == 0:
        return Geometry(_empty_polygon())
    var acc: Geometry = geoms[0].copy()
    var i = 1
    while i < geoms.__len__():
        acc = union(acc, geoms[i])
        i += 1
    return acc.copy()


fn symmetric_difference(a: Polygon, b: Polygon) -> Geometry:
    return overlay_xor(a, b)


fn difference(a: MultiPolygon, p: Polygon) -> Geometry:
    var parts = List[Polygon]()
    for q in a.polys:
        var dg = difference(q, p)
        if dg.is_multipolygon():
            var mp = dg.as_multipolygon()
            for x in mp.polys:
                if not _is_empty_polygon(x):
                    parts.append(x.copy())
    if parts.__len__() == 0:
        return Geometry(_empty_polygon())
    if parts.__len__() == 1:
        return Geometry(parts[0].copy())
    return Geometry(MultiPolygon(parts))


fn difference(p: Polygon, b: MultiPolygon) -> Geometry:
    var acc = p.copy()
    for q in b.polys:
        var dg = difference(acc, q)
        if not dg.is_multipolygon():
            return Geometry(_empty_polygon())
        var mp = dg.as_multipolygon()
        if mp.polys.__len__() == 0:
            return Geometry(_empty_polygon())
        acc = mp.polys[0].copy()
    return Geometry(acc.copy())


fn difference(a: MultiPolygon, b: MultiPolygon) -> Geometry:
    var acc = List[Polygon]()
    for p in a.polys:
        for q in b.polys:
            var pg2 = difference(p, q)
            if pg2.is_multipolygon():
                var mp2 = pg2.as_multipolygon()
                for x in mp2.polys:
                    if not _is_empty_polygon(x):
                        acc.append(x.copy())
    if acc.__len__() == 0:
        return Geometry(_empty_polygon())
    if acc.__len__() == 1:
        return Geometry(acc[0].copy())
    return Geometry(MultiPolygon(acc))


fn symmetric_difference(a: MultiPolygon, p: Polygon) -> Geometry:
    # fold XOR across all polygons: ((p1 XOR p) XOR p2) ...
    if a.polys.__len__() == 0:
        return Geometry(p.copy())
    var acc: Geometry = symmetric_difference(a.polys[0], p)
    var i = 1
    while i < a.polys.__len__():
        acc = symmetric_difference(acc, a.polys[i])
        i += 1
    return acc.copy()


fn symmetric_difference(p: Polygon, a: MultiPolygon) -> Geometry:
    return symmetric_difference(a, p)


fn symmetric_difference(a: MultiPolygon, b: MultiPolygon) -> Geometry:
    if a.polys.__len__() == 0:
        if b.polys.__len__() == 0:
            return Geometry(_empty_polygon())
        var acc: Geometry = Geometry(b.polys[0].copy())
        var j = 1
        while j < b.polys.__len__():
            acc = symmetric_difference(acc, b.polys[j])
            j += 1
        return acc.copy()
    var acc2: Geometry = Geometry(a.polys[0].copy())
    var i = 1
    while i < a.polys.__len__():
        acc2 = symmetric_difference(acc2, a.polys[i])
        i += 1
    var k = 0
    while k < b.polys.__len__():
        acc2 = symmetric_difference(acc2, b.polys[k])
        k += 1
    return acc2.copy()
