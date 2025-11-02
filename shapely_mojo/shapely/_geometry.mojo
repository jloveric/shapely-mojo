struct _GeomCodes:
    var POINT: Int32
    var LINESTRING: Int32
    var LINEARRING: Int32
    var POLYGON: Int32
    var MULTIPOINT: Int32
    var MULTILINESTRING: Int32
    var MULTIPOLYGON: Int32
    var GEOMETRYCOLLECTION: Int32
    var MISSING: Int32

    fn __init__(self):
        self.POINT = 0
        self.LINESTRING = 1
        self.LINEARRING = 2
        self.POLYGON = 3
        self.MULTIPOINT = 4
        self.MULTILINESTRING = 5
        self.MULTIPOLYGON = 6
        self.GEOMETRYCOLLECTION = 7
        self.MISSING = 255

var GeometryType = _GeomCodes()


struct Geometry:
    var geom_type: Int32

    fn __init__(self, geom_type: Int32 = GeometryType.MISSING):
        self.geom_type = geom_type

    fn is_empty(self) -> Bool:
        return False

    fn to_wkt(self) -> String:
        return "GEOMETRYCOLLECTION EMPTY"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        return (0.0, 0.0, 0.0, 0.0)


struct GEOSException:
    fn __init__(self):
        return


fn geos_version() -> (Int32, Int32, Int32):
    return (0, 0, 0)


fn geos_version_string() -> String:
    return "0.0.0"


fn geos_capi_version() -> (Int32, Int32, Int32):
    return (0, 0, 0)


fn geos_capi_version_string() -> String:
    return "0.0.0"
