from std.utils.variant import Variant
from shapely.geometry import (
    Point,
    LineString,
    Polygon,
    GeometryCollection,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
)


comptime GeometryCollectionMember = Variant[
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
]

comptime GeometryPayload = Variant[
    Point,
    LineString,
    Polygon,
    GeometryCollection,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
]


def geometry_from_member(member: GeometryCollectionMember) -> Geometry:
    if member.isa[Point]():
        return Geometry(member[Point].copy())
    if member.isa[LineString]():
        return Geometry(member[LineString].copy())
    if member.isa[Polygon]():
        return Geometry(member[Polygon].copy())
    if member.isa[MultiPoint]():
        return Geometry(member[MultiPoint].copy())
    if member.isa[MultiLineString]():
        return Geometry(member[MultiLineString].copy())
    if member.isa[MultiPolygon]():
        return Geometry(member[MultiPolygon].copy())
    return Geometry(Point(0.0, 0.0))


def _geometry_to_member(g: Geometry) -> GeometryCollectionMember:
    if g.is_point():
        return GeometryCollectionMember(g.as_point().copy())
    if g.is_linestring():
        return GeometryCollectionMember(g.as_linestring().copy())
    if g.is_polygon():
        return GeometryCollectionMember(g.as_polygon().copy())
    if g.is_multipoint():
        return GeometryCollectionMember(g.as_multipoint().copy())
    if g.is_multilinestring():
        return GeometryCollectionMember(g.as_multilinestring().copy())
    if g.is_multipolygon():
        return GeometryCollectionMember(g.as_multipolygon().copy())
    return GeometryCollectionMember(Point(0.0, 0.0))


struct Geometry(Copyable, Movable, Deinitable):
    var payload: GeometryPayload

    def __init__(out self, var value: Point):
        self.payload = GeometryPayload(value.copy())

    def __init__(out self, var value: LineString):
        self.payload = GeometryPayload(value.copy())

    def __init__(out self, var value: Polygon):
        self.payload = GeometryPayload(value.copy())

    def __init__(out self, var value: GeometryCollection):
        self.payload = GeometryPayload(value.copy())

    def __init__(out self, var value: MultiPoint):
        self.payload = GeometryPayload(value.copy())

    def __init__(out self, var value: MultiLineString):
        self.payload = GeometryPayload(value.copy())

    def __init__(out self, var value: MultiPolygon):
        self.payload = GeometryPayload(value.copy())

    def is_point(self) -> Bool:
        return self.payload.isa[Point]()

    def is_linestring(self) -> Bool:
        return self.payload.isa[LineString]()

    def is_polygon(self) -> Bool:
        return self.payload.isa[Polygon]()

    def is_multipoint(self) -> Bool:
        return self.payload.isa[MultiPoint]()

    def is_geometrycollection(self) -> Bool:
        return self.payload.isa[GeometryCollection]()

    def is_multilinestring(self) -> Bool:
        return self.payload.isa[MultiLineString]()

    def is_multipolygon(self) -> Bool:
        return self.payload.isa[MultiPolygon]()

    def as_point(self) -> Point:
        return self.payload[Point].copy()

    def as_linestring(self) -> LineString:
        return self.payload[LineString].copy()

    def as_polygon(self) -> Polygon:
        return self.payload[Polygon].copy()

    def as_geometrycollection(self) -> GeometryCollection:
        return self.payload[GeometryCollection].copy()

    def as_multipoint(self) -> MultiPoint:
        return self.payload[MultiPoint].copy()

    def as_multilinestring(self) -> MultiLineString:
        return self.payload[MultiLineString].copy()

    def as_multipolygon(self) -> MultiPolygon:
        return self.payload[MultiPolygon].copy()

    def is_empty(self) -> Bool:
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

    def to_wkt(self) -> String:
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

    def bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
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

    def area(self) -> Float64:
        if self.payload.isa[Polygon]():
            return self.payload[Polygon].area()
        if self.payload.isa[MultiPolygon]():
            return self.payload[MultiPolygon].area()
        if self.payload.isa[GeometryCollection]():
            var gc = self.payload[GeometryCollection].copy()
            var s = 0.0
            for p in gc.geoms:
                s += geometry_from_member(p).area()
            return s
        return 0.0

    def length(self) -> Float64:
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
            for p in gc.geoms:
                s += geometry_from_member(p).length()
            return s
        return 0.0


struct GEOSException:
    def __init__(out self):
        return


def geos_version() -> Tuple[Int32, Int32, Int32]:
    return (0, 0, 0)


def geos_version_string() -> String:
    return "0.0.0"


def geos_capi_version() -> Tuple[Int32, Int32, Int32]:
    return (0, 0, 0)


def geos_capi_version_string() -> String:
    return "0.0.0"
