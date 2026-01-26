from shapely._geometry import Geometry
from shapely.geometry import LinearRing, Polygon
from shapely.constructive import buffer, JOIN_BEVEL, JOIN_MITRE
from shapely.validation import make_valid


fn _base_poly() -> Geometry:
    var shell = List[Tuple[Float64, Float64]]()
    shell.append((0.0, 0.0))
    shell.append((10.0, 0.0))
    shell.append((10.0, 7.0))
    shell.append((0.0, 7.0))
    shell.append((0.0, 0.0))
    return Geometry(Polygon(LinearRing(shell)))


fn main():
    var g = _base_poly()
    var d = 1.2

    var b1 = buffer(g.copy(), d, 16, 1, Int32(JOIN_BEVEL), 5.0)
    print("bevel type:", b1.to_wkt()[0:12])
    var mv1 = make_valid(b1.copy())
    print("bevel valid type:", mv1.to_wkt()[0:12])

    var b2 = buffer(g.copy(), d, 16, 1, Int32(JOIN_MITRE), 5.0)
    print("mitre type:", b2.to_wkt()[0:12])
    var mv2 = make_valid(b2.copy())
    print("mitre valid type:", mv2.to_wkt()[0:12])
