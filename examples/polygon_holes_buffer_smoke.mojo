from shapely._geometry import Geometry
from shapely.geometry import LinearRing, Polygon
from shapely.constructive import buffer, circle, JOIN_BEVEL, JOIN_MITRE, JOIN_ROUND


fn _poly_with_holes() -> Geometry:
    var shell = List[Tuple[Float64, Float64]]()
    shell.append((0.0, 0.0))
    shell.append((10.0, 0.0))
    shell.append((10.0, 7.0))
    shell.append((0.0, 7.0))
    shell.append((0.0, 0.0))

    var holes = List[LinearRing]()
    holes.append(circle(3.0, 2.0, 0.9, 16).as_polygon().shell.copy())
    holes.append(circle(7.0, 2.2, 1.1, 16).as_polygon().shell.copy())
    holes.append(circle(5.0, 5.0, 0.7, 16).as_polygon().shell.copy())

    return Geometry(Polygon(LinearRing(shell), holes))


fn main():
    var g = _poly_with_holes()
    var d = 1.2

    var b0 = buffer(g.copy(), d, 16, 1, Int32(JOIN_ROUND), 5.0)
    print("round:", b0.to_wkt()[0:20])

    var b1 = buffer(g.copy(), d, 16, 1, Int32(JOIN_BEVEL), 5.0)
    print("bevel:", b1.to_wkt()[0:20])

    var b2 = buffer(g.copy(), d, 16, 1, Int32(JOIN_MITRE), 5.0)
    print("mitre:", b2.to_wkt()[0:20])
