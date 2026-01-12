from shapely._geometry import Geometry
from shapely.geometry import LineString, MultiLineString, LinearRing, Polygon, MultiPolygon, GeometryCollection
from shapely.ops import polygonize_full


fn _close_ring_coords(coords: List[Tuple[Float64, Float64]]) -> List[Tuple[Float64, Float64]]:
    if coords.__len__() == 0:
        return coords.copy()
    var out = coords.copy()
    var first = out[0]
    var last = out[out.__len__() - 1]
    if first[0] != last[0] or first[1] != last[1]:
        out.append(first)
    return out.copy()


fn _make_valid_polygon(p: Polygon) -> Polygon:
    var shell_coords = _close_ring_coords(p.shell.coords)
    var shell = LinearRing(shell_coords)

    var holes = List[LinearRing]()
    for h in p.holes:
        holes.append(LinearRing(_close_ring_coords(h.coords)))

    return Polygon(shell, holes)


fn _polygonize_boundary(p: Polygon) -> Geometry:
    # Build linework from polygon boundary and polygonize it to resolve self-intersections.
    var segs = List[LineString]()

    fn add_ring(coords_in: List[Tuple[Float64, Float64]], mut segs: List[LineString]):
        var coords = _close_ring_coords(coords_in)
        if coords.__len__() < 4:
            return
        var i = 0
        while i < coords.__len__() - 1:
            var pts = List[Tuple[Float64, Float64]]()
            pts.append(coords[i])
            pts.append(coords[i + 1])
            segs.append(LineString(pts))
            i += 1

    add_ring(p.shell.coords, segs)
    for h in p.holes:
        add_ring(h.coords, segs)

    if segs.__len__() == 0:
        return Geometry(GeometryCollection([]))

    var mls = MultiLineString(segs)
    var res = polygonize_full(Geometry(mls.copy()))
    ref polys = res[0]

    var out_polys = List[Polygon]()
    for g in polys.geoms:
        if g.is_polygon():
            out_polys.append(g.as_polygon())
        elif g.is_multipolygon():
            for pp in g.as_multipolygon().polys:
                out_polys.append(pp.copy())

    fn _try_split_self_touch(shell: LinearRing) -> List[LinearRing]:
        var coords = _close_ring_coords(shell.coords)
        if coords.__len__() < 7:
            return List[LinearRing]()
        # Find a repeated vertex (excluding the closing vertex).
        var seen_pt = List[Tuple[Float64, Float64, Int]]()
        var i = 0
        while i < coords.__len__() - 1:
            var x = coords[i][0]
            var y = coords[i][1]
            var j = 0
            var found = -1
            while j < seen_pt.__len__():
                if seen_pt[j][0] == x and seen_pt[j][1] == y:
                    found = seen_pt[j][2]
                    break
                j += 1
            if found != -1 and i - found >= 3:
                # Split into two rings at (found -> i)
                var r1 = List[Tuple[Float64, Float64]]()
                var a = found
                while a <= i:
                    r1.append(coords[a])
                    a += 1
                if r1.__len__() > 0:
                    var f1 = r1[0]
                    var l1 = r1[r1.__len__() - 1]
                    if f1[0] != l1[0] or f1[1] != l1[1]:
                        r1.append(f1)

                var r2 = List[Tuple[Float64, Float64]]()
                a = i
                while a < coords.__len__() - 1:
                    r2.append(coords[a])
                    a += 1
                a = 0
                while a <= found:
                    r2.append(coords[a])
                    a += 1
                if r2.__len__() > 0:
                    var f2 = r2[0]
                    var l2 = r2[r2.__len__() - 1]
                    if f2[0] != l2[0] or f2[1] != l2[1]:
                        r2.append(f2)

                var out = List[LinearRing]()
                if r1.__len__() >= 4:
                    out.append(LinearRing(r1))
                if r2.__len__() >= 4:
                    out.append(LinearRing(r2))
                return out.copy()
            seen_pt.append((x, y, i))
            i += 1
        return List[LinearRing]()

    if out_polys.__len__() == 0:
        return Geometry(GeometryCollection([]))
    if out_polys.__len__() == 1:
        # If polygonization produced a single self-touching ring (e.g. bowtie), split it.
        var p0 = out_polys[0].copy()
        var split = _try_split_self_touch(p0.shell)
        if split.__len__() == 2:
            return Geometry(MultiPolygon([Polygon(split[0].copy()), Polygon(split[1].copy())]))
        return Geometry(p0.copy())
    return Geometry(MultiPolygon(out_polys))


fn make_valid(geom: Geometry) -> Geometry:
    # Best-effort Shapely-like make_valid for polygons: close rings and
    # polygonize boundary linework to resolve self-intersections.
    if geom.is_polygon():
        var p = _make_valid_polygon(geom.as_polygon())
        return _polygonize_boundary(p)
    if geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        var out = List[Polygon]()
        for p in mp.polys:
            var fixed = _make_valid_polygon(p.copy())
            var g = _polygonize_boundary(fixed)
            if g.is_polygon():
                out.append(g.as_polygon())
            elif g.is_multipolygon():
                for pp in g.as_multipolygon().polys:
                    out.append(pp.copy())
        if out.__len__() == 0:
            return Geometry(GeometryCollection([]))
        return Geometry(MultiPolygon(out))
    return geom.copy()
