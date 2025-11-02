from shapely._geometry import Geometry


fn affine_transform[T: Geometry](geom: T, _matrix: Tuple[Float64, ...]) -> T:
    return geom


fn rotate[T: Geometry](geom: T, _angle: Float64, origin: String = "center", use_radians: Bool = False) -> T:
    return geom


fn scale[T: Geometry](geom: T, xfact: Float64 = 1.0, yfact: Float64 = 1.0, zfact: Float64 = 1.0, origin: String = "center") -> T:
    return geom


fn skew[T: Geometry](geom: T, xs: Float64 = 0.0, ys: Float64 = 0.0, origin: String = "center", use_radians: Bool = False) -> T:
    return geom


fn translate[T: Geometry](geom: T, xoff: Float64 = 0.0, yoff: Float64 = 0.0, zoff: Float64 = 0.0) -> T:
    return geom
