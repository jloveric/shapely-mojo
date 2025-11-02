from shapely._geometry import Geometry
from shapely.geometry import GeometryCollection, LinearRing, Polygon, MultiPolygon, LineString, Point
from shapely.algorithms import orientation, point_in_polygon
from shapely.overlay import overlay_union, overlay_difference, overlay_intersection, overlay_xor


fn _empty_polygon() -> Polygon:
    return Polygon(LinearRing(List[Tuple[Float64, Float64]]()))


fn _is_empty_polygon(p: Polygon) -> Bool:
    return p.shell.coords.size() < 3


fn union(a: Geometry, b: Geometry) -> Geometry:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Polygon":
        return overlay_union(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    if ta == "MultiPolygon" and (tb == "Polygon" or tb == "MultiPolygon"):
        var acc: Geometry = unsafe_bitcast[Geometry](a)
        if tb == "Polygon":
            let p = unsafe_bitcast[Polygon](b)
            # fold: union acc with polygon p
            var res = GeometryCollection([]) as Geometry  # placeholder init
            # start with first polygon
            let mp = unsafe_bitcast[MultiPolygon](acc)
            if mp.polys.size() == 0:
                return p
            var g: Geometry = overlay_union(mp.polys[0], p)
            var i = 1
            while i < mp.polys.size():
                g = union(g, mp.polys[i])
                i += 1
            return g
        else:
            let mpb = unsafe_bitcast[MultiPolygon](b)
            var g2: Geometry = unsafe_bitcast[Geometry](a)
            var j = 0
            while j < mpb.polys.size():
                g2 = union(g2, mpb.polys[j])
                j += 1
            return g2
    if tb == "MultiPolygon" and ta == "Polygon":
        # symmetric to case above
        var g3: Geometry = unsafe_bitcast[Geometry](b)
        return union(g3, a)
    if ta == "MultiPolygon" and tb == "MultiPolygon":
        var g4: Geometry = unsafe_bitcast[Geometry](a)
        let mpb2 = unsafe_bitcast[MultiPolygon](b)
        var k = 0
        while k < mpb2.polys.size():
            g4 = union(g4, mpb2.polys[k])
            k += 1
        return g4
    return GeometryCollection([a, b])


fn union(a: Polygon, b: Polygon) -> Geometry:
    return overlay_union(a, b)


fn union(mp: MultiPolygon, p: Polygon) -> Geometry:
    if mp.polys.size() == 0: return p
    var g: Geometry = overlay_union(mp.polys[0], p)
    var i = 1
    while i < mp.polys.size():
        g = union(g, mp.polys[i])
        i += 1
    return g


fn union(p: Polygon, mp: MultiPolygon) -> Geometry:
    return union(mp, p)


fn union(a: MultiPolygon, b: MultiPolygon) -> Geometry:
    if a.polys.size() == 0:
        if b.polys.size() == 0: return MultiPolygon([])
        # return b as Geometry
        var g0: Geometry = b.polys[0]
        var jj = 1
        while jj < b.polys.size():
            g0 = union(g0, b.polys[jj])
            jj += 1
        return g0
    var g: Geometry = a.polys[0]
    var i = 1
    while i < a.polys.size():
        g = union(g, a.polys[i])
        i += 1
    var j = 0
    while j < b.polys.size():
        g = union(g, b.polys[j])
        j += 1
    return g


fn _cross(ax: Float64, ay: Float64, bx: Float64, by: Float64) -> Float64:
    return ax * by - ay * bx


fn _intersect_point(s1: Tuple[Float64, Float64], s2: Tuple[Float64, Float64], c1: Tuple[Float64, Float64], c2: Tuple[Float64, Float64]) -> (Tuple[Float64, Float64], Bool):
    let r_x = s2[0] - s1[0]
    let r_y = s2[1] - s1[1]
    let s_x = c2[0] - c1[0]
    let s_y = c2[1] - c1[1]
    let denom = _cross(r_x, r_y, s_x, s_y)
    if denom == 0.0:
        return ((0.0, 0.0), False)
    let t = _cross(c1[0] - s1[0], c1[1] - s1[1], s_x, s_y) / denom
    return ((s1[0] + t * r_x, s1[1] + t * r_y), True)


fn _is_inside(p: Tuple[Float64, Float64], a: Tuple[Float64, Float64], b: Tuple[Float64, Float64]) -> Bool:
    # left-of test: inside if p is to left of a->b
    return orientation(a[0], a[1], b[0], b[1], p[0], p[1]) >= 0


fn _suth_hodg(subject: List[Tuple[Float64, Float64]], clip: List[Tuple[Float64, Float64]]) -> List[Tuple[Float64, Float64]]:
    var output = subject
    if output.size() == 0: return output
    var i = 0
    while i < clip.size() - 1:
        let A = clip[i]
        let B = clip[i + 1]
        var input = output
        output = List[Tuple[Float64, Float64]]()
        if input.size() == 0: break
        var S = input[input.size() - 1]
        for E in input:
            if _is_inside(E, A, B):
                if not _is_inside(S, A, B):
                    let (I, ok) = _intersect_point(S, E, A, B)
                    if ok: output.push_back(I)
                output.push_back(E)
            else:
                if _is_inside(S, A, B):
                    let (I2, ok2) = _intersect_point(S, E, A, B)
                    if ok2: output.push_back(I2)
            S = E
        i += 1
    return output


fn intersection(a: Polygon, b: Polygon) -> Geometry:
    return overlay_intersection(a, b)


fn intersection(a: Geometry, b: Geometry) -> Geometry:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Polygon":
        return overlay_intersection(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    if ta == "MultiPolygon" and tb == "Polygon":
        return intersection(unsafe_bitcast[MultiPolygon](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "MultiPolygon":
        return intersection(unsafe_bitcast[Polygon](a), unsafe_bitcast[MultiPolygon](b))
    if ta == "MultiPolygon" and tb == "MultiPolygon":
        return intersection(unsafe_bitcast[MultiPolygon](a), unsafe_bitcast[MultiPolygon](b))
    return GeometryCollection([a, b])


fn difference(a: Polygon, b: Polygon) -> Geometry:
    return overlay_difference(a, b)


fn difference(a: Geometry, b: Geometry) -> Geometry:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Polygon":
        return difference(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    if ta == "MultiPolygon" and tb == "Polygon":
        return difference(unsafe_bitcast[MultiPolygon](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "MultiPolygon":
        return difference(unsafe_bitcast[Polygon](a), unsafe_bitcast[MultiPolygon](b))
    if ta == "MultiPolygon" and tb == "MultiPolygon":
        return difference(unsafe_bitcast[MultiPolygon](a), unsafe_bitcast[MultiPolygon](b))
    return a


fn symmetric_difference(a: Geometry, b: Geometry) -> Geometry:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Polygon":
        return symmetric_difference(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    if ta == "MultiPolygon" and tb == "Polygon":
        return symmetric_difference(unsafe_bitcast[MultiPolygon](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "MultiPolygon":
        return symmetric_difference(unsafe_bitcast[Polygon](a), unsafe_bitcast[MultiPolygon](b))
    if ta == "MultiPolygon" and tb == "MultiPolygon":
        return symmetric_difference(unsafe_bitcast[MultiPolygon](a), unsafe_bitcast[MultiPolygon](b))
    return GeometryCollection([a, b])


fn unary_union(geoms: List[Geometry]) -> Geometry:
    var polys = List[Polygon]()
    var others = List[Geometry]()
    for g in geoms:
        let t = g.__type_name__()
        if t == "Polygon":
            polys.push_back(unsafe_bitcast[Polygon](g))
        elif t == "MultiPolygon":
            for p in unsafe_bitcast[MultiPolygon](g).polys:
                polys.push_back(p)
        else:
            others.push_back(g)
    var result: Geometry = GeometryCollection([])
    if polys.size() > 0:
        var acc: Geometry = polys[0]
        var i = 1
        while i < polys.size():
            acc = union(acc, polys[i])
            i += 1
        result = acc
    if others.size() == 0:
        return result
    if result.__type_name__() == "Polygon" or result.__type_name__() == "MultiPolygon":
        var items = List[Geometry]()
        if result.__type_name__() == "Polygon":
            items.push_back(result)
        else:
            let mpr = unsafe_bitcast[MultiPolygon](result)
            items.push_back(mpr)
        for og in others: items.push_back(og)
        return GeometryCollection(items)
    return GeometryCollection(others)


fn symmetric_difference(a: Polygon, b: Polygon) -> Geometry:
    return overlay_xor(a, b)


fn intersection(a: MultiPolygon, p: Polygon) -> Geometry:
    var parts = List[Polygon]()
    for q in a.polys:
        let g = overlay_intersection(q, p)
        if g.__type_name__() == "Polygon":
            let pg = unsafe_bitcast[Polygon](g)
            if not _is_empty_polygon(pg): parts.push_back(pg)
        elif g.__type_name__() == "MultiPolygon":
            for x in unsafe_bitcast[MultiPolygon](g).polys:
                if not _is_empty_polygon(x): parts.push_back(x)
    if parts.size() == 0: return _empty_polygon()
    if parts.size() == 1: return parts[0]
    return MultiPolygon(parts)


fn intersection(p: Polygon, b: MultiPolygon) -> Geometry:
    return intersection(b, p)


fn intersection(a: MultiPolygon, b: MultiPolygon) -> Geometry:
    var parts = List[Polygon]()
    for q in a.polys:
        for r in b.polys:
            let g = overlay_intersection(q, r)
            if g.__type_name__() == "Polygon":
                let pg = unsafe_bitcast[Polygon](g)
                if not _is_empty_polygon(pg): parts.push_back(pg)
            elif g.__type_name__() == "MultiPolygon":
                for x in unsafe_bitcast[MultiPolygon](g).polys:
                    if not _is_empty_polygon(x): parts.push_back(x)
    if parts.size() == 0: return _empty_polygon()
    if parts.size() == 1: return parts[0]
    return MultiPolygon(parts)


fn symmetric_difference(a: MultiPolygon, p: Polygon) -> Geometry:
    # fold XOR across all polygons: ((p1 XOR p) XOR p2) ...
    if a.polys.size() == 0:
        return p
    var acc: Geometry = symmetric_difference(a.polys[0], p)
    var i = 1
    while i < a.polys.size():
        acc = symmetric_difference(acc, a.polys[i])
        i += 1
    return acc


fn symmetric_difference(p: Polygon, a: MultiPolygon) -> Geometry:
    return symmetric_difference(a, p)


fn symmetric_difference(a: MultiPolygon, b: MultiPolygon) -> Geometry:
    if a.polys.size() == 0:
        if b.polys.size() == 0: return _empty_polygon()
        var acc: Geometry = b.polys[0]
        var j = 1
        while j < b.polys.size():
            acc = symmetric_difference(acc, b.polys[j])
            j += 1
        return acc
    var acc2: Geometry = a.polys[0]
    var i = 1
    while i < a.polys.size():
        acc2 = symmetric_difference(acc2, a.polys[i])
        i += 1
    var k = 0
    while k < b.polys.size():
        acc2 = symmetric_difference(acc2, b.polys[k])
        k += 1
    return acc2


fn difference(a: MultiPolygon, p: Polygon) -> Geometry:
    var parts = List[Polygon]()
    for q in a.polys:
        let dg = difference(q, p)
        if dg.__type_name__() == "Polygon":
            let d = unsafe_bitcast[Polygon](dg)
            if not _is_empty_polygon(d): parts.push_back(d)
        elif dg.__type_name__() == "MultiPolygon":
            for x in unsafe_bitcast[MultiPolygon](dg).polys:
                if not _is_empty_polygon(x): parts.push_back(x)
    if parts.size() == 0: return _empty_polygon()
    if parts.size() == 1: return parts[0]
    return MultiPolygon(parts)


fn difference(p: Polygon, b: MultiPolygon) -> Geometry:
    var acc = p
    for q in b.polys:
        let dg = difference(acc, q)
        if dg.__type_name__() == "Polygon":
            acc = unsafe_bitcast[Polygon](dg)
            if _is_empty_polygon(acc): return _empty_polygon()
        elif dg.__type_name__() == "MultiPolygon":
            # collapsing MultiPolygon by returning first polygon; future: retain full result
            let mp = unsafe_bitcast[MultiPolygon](dg)
            if mp.polys.size() == 0: return _empty_polygon()
            acc = mp.polys[0]
    return acc


fn difference(a: MultiPolygon, b: MultiPolygon) -> Geometry:
    var acc = List[Polygon]()
    for p in a.polys:
        var pg = p as Geometry
        for q in b.polys:
            pg = difference(unsafe_bitcast[Polygon](pg), q)
        if pg.__type_name__() == "Polygon":
            let pp = unsafe_bitcast[Polygon](pg)
            if not _is_empty_polygon(pp): acc.push_back(pp)
        elif pg.__type_name__() == "MultiPolygon":
            for x in unsafe_bitcast[MultiPolygon](pg).polys:
                if not _is_empty_polygon(x): acc.push_back(x)
    if acc.size() == 0: return _empty_polygon()
    if acc.size() == 1: return acc[0]
    return MultiPolygon(acc)
