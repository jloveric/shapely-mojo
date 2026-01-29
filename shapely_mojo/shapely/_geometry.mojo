from utils.variant import Variant
from shapely.geometry import (
    Point,
    LineString,
    Polygon,
    GeometryCollection,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
)


alias GeometryPayload = Variant[
    Point,
    LineString,
    Polygon,
    GeometryCollection,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
]


struct Geometry(Copyable, Movable):
    var payload: GeometryPayload

    fn __init__(out self, var value: Point):
        self.payload = GeometryPayload(value.copy())

    fn __init__(out self, var value: LineString):
        self.payload = GeometryPayload(value.copy())

    fn __init__(out self, var value: Polygon):
        self.payload = GeometryPayload(value.copy())

    fn __init__(out self, var value: GeometryCollection):
        self.payload = GeometryPayload(value.copy())

    fn __init__(out self, var value: MultiPoint):
        self.payload = GeometryPayload(value.copy())

    fn __init__(out self, var value: MultiLineString):
        self.payload = GeometryPayload(value.copy())

    fn __init__(out self, var value: MultiPolygon):
        self.payload = GeometryPayload(value.copy())

    fn is_point(self) -> Bool:
        return self.payload.isa[Point]()

    fn is_linestring(self) -> Bool:
        return self.payload.isa[LineString]()

    fn is_polygon(self) -> Bool:
        return self.payload.isa[Polygon]()

    fn is_multipoint(self) -> Bool:
        return self.payload.isa[MultiPoint]()

    fn is_geometrycollection(self) -> Bool:
        return self.payload.isa[GeometryCollection]()

    fn is_multilinestring(self) -> Bool:
        return self.payload.isa[MultiLineString]()

    fn is_multipolygon(self) -> Bool:
        return self.payload.isa[MultiPolygon]()

    fn as_point(self) -> Point:
        return self.payload[Point].copy()

    fn as_linestring(self) -> LineString:
        return self.payload[LineString].copy()

    fn as_polygon(self) -> Polygon:
        return self.payload[Polygon].copy()

    fn as_geometrycollection(self) -> GeometryCollection:
        return self.payload[GeometryCollection].copy()

    fn as_multipoint(self) -> MultiPoint:
        return self.payload[MultiPoint].copy()

    fn as_multilinestring(self) -> MultiLineString:
        return self.payload[MultiLineString].copy()

    fn as_multipolygon(self) -> MultiPolygon:
        return self.payload[MultiPolygon].copy()

    fn is_empty(self) -> Bool:
        if self.payload.isa[Point]():
            return self.payload[Point].is_empty()
        if self.payload.isa[LineString]():
            return self.payload[LineString].is_empty()
        if self.payload.isa[Polygon]():
            return self.payload[Polygon].is_empty()
        if self.payload.isa[GeometryCollection]():
            return self.payload[GeometryCollection].is_empty()
        if self.payload.isa[MultiPoint]():
            return self.payload[MultiPoint].is_empty()
        if self.payload.isa[MultiLineString]():
            return self.payload[MultiLineString].is_empty()
        if self.payload.isa[MultiPolygon]():
            return self.payload[MultiPolygon].is_empty()
        return False

    fn to_wkt(self) -> String:
        if self.payload.isa[Point]():
            return self.payload[Point].to_wkt()
        if self.payload.isa[LineString]():
            return self.payload[LineString].to_wkt()
        if self.payload.isa[Polygon]():
            return self.payload[Polygon].to_wkt()
        if self.payload.isa[GeometryCollection]():
            return self.payload[GeometryCollection].to_wkt()
        if self.payload.isa[MultiPoint]():
            return self.payload[MultiPoint].to_wkt()
        if self.payload.isa[MultiLineString]():
            return self.payload[MultiLineString].to_wkt()
        if self.payload.isa[MultiPolygon]():
            return self.payload[MultiPolygon].to_wkt()
        return "GEOMETRYCOLLECTION EMPTY"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        if self.payload.isa[Point]():
            return self.payload[Point].bounds()
        if self.payload.isa[LineString]():
            return self.payload[LineString].bounds()
        if self.payload.isa[Polygon]():
            return self.payload[Polygon].bounds()
        if self.payload.isa[GeometryCollection]():
            return self.payload[GeometryCollection].bounds()
        if self.payload.isa[MultiPoint]():
            return self.payload[MultiPoint].bounds()
        if self.payload.isa[MultiLineString]():
            return self.payload[MultiLineString].bounds()
        if self.payload.isa[MultiPolygon]():
            return self.payload[MultiPolygon].bounds()
        return (0.0, 0.0, 0.0, 0.0)

    fn area(self) -> Float64:
        if self.payload.isa[Polygon]():
            return self.payload[Polygon].area()
        if self.payload.isa[MultiPolygon]():
            return self.payload[MultiPolygon].area()
        if self.payload.isa[GeometryCollection]():
            var gc = self.payload[GeometryCollection].copy()
            var s = 0.0
            for g in gc.geoms:
                s += g.area()
            return s
        return 0.0

    fn length(self) -> Float64:
        if self.payload.isa[LineString]():
            return self.payload[LineString].length()
        if self.payload.isa[MultiLineString]():
            return self.payload[MultiLineString].length()
        if self.payload.isa[Polygon]():
            return self.payload[Polygon].length()
        if self.payload.isa[MultiPolygon]():
            return self.payload[MultiPolygon].length()
        if self.payload.isa[GeometryCollection]():
            var gc = self.payload[GeometryCollection].copy()
            var s = 0.0
            for g in gc.geoms:
                s += g.length()
            return s
        return 0.0


struct GEOSException:
    fn __init__(out self):
        return


fn geos_version() -> Tuple[Int32, Int32, Int32]:
    return (0, 0, 0)


fn geos_version_string() -> String:
    return "0.0.0"


fn geos_capi_version() -> Tuple[Int32, Int32, Int32]:
    return (0, 0, 0)


fn geos_capi_version_string() -> String:
    return "0.0.0"
