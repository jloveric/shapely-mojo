from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection


def transform[T](geometry: T, _transformation, include_z: Bool = False, *, interleaved: Bool = True) -> T:
    return geometry


def count_coordinates(_geometry) -> Int32:
    return 0


def get_coordinates(_geometry, include_z: Bool = False, return_index: Bool = False, *, include_m: Bool = False):
    return []


def set_coordinates[T](geometry: T, _coordinates) -> T:
    return geometry
