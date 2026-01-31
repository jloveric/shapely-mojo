from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection


fn transform[T](geometry: T, _transformation, include_z: Bool = False, *, interleaved: Bool = True) -> T:
    return geometry


fn count_coordinates(_geometry) -> Int32:
    return 0


fn get_coordinates(_geometry, include_z: Bool = False, return_index: Bool = False, *, include_m: Bool = False):
    return []


fn set_coordinates[T](geometry: T, _coordinates) -> T:
    return geometry
