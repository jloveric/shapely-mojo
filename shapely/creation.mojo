from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection


fn points(coords: Tuple[Float64, Float64]) -> Point:
    return Point(coords[0], coords[1])


fn linestrings(coords: List[Tuple[Float64, Float64]]) -> LineString:
    return LineString(coords)


fn linearrings(coords: List[Tuple[Float64, Float64]]) -> LinearRing:
    # Ensure closed ring (repeat start if needed)
    if coords.__len__() > 0:
        var first = coords[0]
        var last = coords[coords.__len__() - 1]
        if first[0] != last[0] or first[1] != last[1]:
            var closed = List[Tuple[Float64, Float64]]()
            for c in coords: closed.append(c)
            closed.append(first)
            return LinearRing(closed)
    return LinearRing(coords)


fn polygons(shell_coords: List[Tuple[Float64, Float64]], holes: List[List[Tuple[Float64, Float64]]] = List[List[Tuple[Float64, Float64]]]()) -> Polygon:
    var shell = linearrings(shell_coords)
    var ring_holes = List[LinearRing]()
    for h in holes:
        ring_holes.append(linearrings(h))
    return Polygon(shell, ring_holes)


fn multipoints(points_in: List[Point]) -> MultiPoint:
    return MultiPoint(points_in)


fn multilinestrings(lines: List[LineString]) -> MultiLineString:
    return MultiLineString(lines)


fn multipolygons(polys: List[Polygon]) -> MultiPolygon:
    return MultiPolygon(polys)


fn geometrycollections(geoms: List[Geometry]) -> GeometryCollection:
    return GeometryCollection(geoms)


fn box(xmin: Float64, ymin: Float64, xmax: Float64, ymax: Float64, ccw: Bool = True) -> Polygon:
    if ccw:
        return polygons([(xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin), (xmax, ymin)])
    else:
        return polygons([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)])


fn prepare(_geometry) -> None:
    return


fn destroy_prepared(_geometry) -> None:
    return


fn empty_point_array(n: Int32) -> List[Point]:
    return List[Point]()
