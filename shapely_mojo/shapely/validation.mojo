from shapely._geometry import Geometry
from shapely.geometry import LinearRing, Polygon, MultiPolygon


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


fn make_valid(geom: Geometry) -> Geometry:
    # Minimal placeholder: ensure polygon rings are explicitly closed.
    # This is enough for many downstream computations that assume closure.
    if geom.is_polygon():
        return Geometry(_make_valid_polygon(geom.as_polygon()))
    if geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        var polys = List[Polygon]()
        for p in mp.polys:
            polys.append(_make_valid_polygon(p.copy()))
        return Geometry(MultiPolygon(polys))
    return geom.copy()
