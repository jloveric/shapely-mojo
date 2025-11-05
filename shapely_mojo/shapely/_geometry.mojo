struct Geometry(Copyable, Movable):
    fn __init__(out self):
        return

    fn is_empty(self) -> Bool:
        return False

    fn to_wkt(self) -> String:
        return "GEOMETRYCOLLECTION EMPTY"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        return (0.0, 0.0, 0.0, 0.0)


struct GEOSException:
    fn __init__(out self):
        return


fn geos_version() -> (Int32, Int32, Int32):
    return (0, 0, 0)


fn geos_version_string() -> String:
    return "0.0.0"


fn geos_capi_version() -> (Int32, Int32, Int32):
    return (0, 0, 0)


fn geos_capi_version_string() -> String:
    return "0.0.0"
