from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection


def points(coords: Tuple[Float64, Float64]) -> Point:
    return Point(coords[0], coords[1])


def linestrings(coords: List[Tuple[Float64, Float64]]) -> LineString:
    return LineString(coords)


def linearrings(coords: List[Tuple[Float64, Float64]]) -> LinearRing:
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


def polygons(shell_coords: List[Tuple[Float64, Float64]], holes: List[List[Tuple[Float64, Float64]]] = List[List[Tuple[Float64, Float64]]]()) -> Polygon:
    var shell = linearrings(shell_coords)
    var ring_holes = List[LinearRing]()
    for h in holes:
        ring_holes.append(linearrings(h))
    return Polygon(shell, ring_holes)


def multipoints(points_in: List[Point]) -> MultiPoint:
    return MultiPoint(points_in)


def multilinestrings(lines: List[LineString]) -> MultiLineString:
    return MultiLineString(lines)


def multipolygons(polys: List[Polygon]) -> MultiPolygon:
    return MultiPolygon(polys)


def geometrycollections(geoms: List[Geometry]) -> GeometryCollection:
    return GeometryCollection(geoms)


def box(xmin: Float64, ymin: Float64, xmax: Float64, ymax: Float64, ccw: Bool = True) -> Polygon:
    if ccw:
        return polygons([(xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin), (xmax, ymin)])
    else:
        return polygons([(xmin, ymin), (xmin, ymax), (xmax, ymax), (xmax, ymin), (xmin, ymin)])


def prepare(_geometry) -> None:
    return


def destroy_prepared(_geometry) -> None:
    return


def empty_point_array(n: Int32) -> List[Point]:
    return List[Point]()
