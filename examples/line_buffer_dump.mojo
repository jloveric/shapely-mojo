from shapely._geometry import Geometry
from shapely.geometry import LineString, MultiLineString
from shapely.constructive import buffer


fn _coord_list_json(coords: List[Tuple[Float64, Float64]]) -> String:
    var s = "["
    var i = 0
    while i < coords.__len__():
        if i != 0:
            s = s + ","
        s = s + "[" + coords[i][0].__str__() + "," + coords[i][1].__str__() + "]"
        i += 1
    return s + "]"


fn _lines_json(lines: List[LineString]) -> String:
    var s = "["
    var i = 0
    while i < lines.__len__():
        if i != 0:
            s = s + ","
        s = s + _coord_list_json(lines[i].coords)
        i += 1
    return s + "]"


fn _polygons_json(geom: Geometry) -> String:
    # polygons serialized as a list of shells
    var s = "["
    if geom.is_polygon():
        var poly = geom.as_polygon()
        s = s + _coord_list_json(poly.shell.coords)
    elif geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        var i = 0
        while i < mp.polys.__len__():
            if i != 0:
                s = s + ","
            s = s + _coord_list_json(mp.polys[i].shell.coords)
            i += 1
    return s + "]"


fn main():
    var l1 = LineString([(0.0, 0.0), (2.0, 0.5), (4.0, 0.0)])
    var l2 = LineString([(0.0, 2.0), (1.0, 3.0), (2.0, 2.5), (4.0, 3.0)])

    var mls = MultiLineString([l1.copy(), l2.copy()])
    var buf = buffer(Geometry(mls.copy()), 0.25)

    var s = "{\"lines\":" + _lines_json(mls.lines) + ",\"buffer\":" + _polygons_json(buf) + "}"
    print(s)
