from shapely._geometry import Geometry, GeometryType


struct Point(Geometry):
    var x: Float64
    var y: Float64
    var has_z: Bool
    var z: Float64

    fn __init__(self, x: Float64, y: Float64):
        self.x = x
        self.y = y
        self.has_z = False
        self.z = 0.0
        super().__init__(GeometryType.POINT)

    fn with_z(self, z: Float64) -> Point:
        var p = Point(self.x, self.y)
        p.has_z = True
        p.z = z
        return p

    fn is_empty(self) -> Bool:
        return False

    fn to_wkt(self) -> String:
        if self.has_z:
            return "POINT Z (" + self.x.__str__() + " " + self.y.__str__() + " " + self.z.__str__() + ")"
        return "POINT (" + self.x.__str__() + " " + self.y.__str__() + ")"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        return (self.x, self.y, self.x, self.y)


struct LineString(Geometry):
    var coords: List[Tuple[Float64, Float64]]

    fn __init__(self, coords: List[Tuple[Float64, Float64]]):
        self.coords = coords
        super().__init__(GeometryType.LINESTRING)

    fn is_empty(self) -> Bool:
        return self.coords.size() == 0

    fn to_wkt(self) -> String:
        var s = "LINESTRING ("
        var first = True
        for c in self.coords:
            if not first: s = s + ", "
            first = False
            s = s + c[0].__str__() + " " + c[1].__str__()
        return s + ")"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        if self.coords.size() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var minx = self.coords[0][0]
        var miny = self.coords[0][1]
        var maxx = minx
        var maxy = miny
        for c in self.coords:
            if c[0] < minx: minx = c[0]
            if c[0] > maxx: maxx = c[0]
            if c[1] < miny: miny = c[1]
            if c[1] > maxy: maxy = c[1]
        return (minx, miny, maxx, maxy)


struct LinearRing(Geometry):
    var coords: List[Tuple[Float64, Float64]]

    fn __init__(self, coords: List[Tuple[Float64, Float64]]):
        self.coords = coords
        super().__init__(GeometryType.LINEARRING)

    fn is_empty(self) -> Bool:
        return self.coords.is_empty()

    fn to_wkt(self) -> String:
        var s = "LINEARRING ("
        var first = True
        for c in self.coords:
            if not first: s = s + ", "
            first = False
            s = s + c[0].__str__() + " " + c[1].__str__()
        return s + ")"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        if self.coords.size() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var minx = self.coords[0][0]
        var miny = self.coords[0][1]
        var maxx = minx
        var maxy = miny
        for c in self.coords:
            if c[0] < minx: minx = c[0]
            if c[0] > maxx: maxx = c[0]
            if c[1] < miny: miny = c[1]
            if c[1] > maxy: maxy = c[1]
        return (minx, miny, maxx, maxy)


struct Polygon(Geometry):
    var shell: LinearRing
    var holes: List[LinearRing]

    fn __init__(self, shell: LinearRing, holes: List[LinearRing] = List[LinearRing]()):
        self.shell = shell
        self.holes = holes
        super().__init__(GeometryType.POLYGON)

    fn is_empty(self) -> Bool:
        return self.shell.coords.size() == 0

    fn to_wkt(self) -> String:
        var s = "POLYGON ("
        # shell
        s = s + "("
        var first = True
        for c in self.shell.coords:
            if not first: s = s + ", "
            first = False
            s = s + c[0].__str__() + " " + c[1].__str__()
        s = s + ")"
        # holes
        for h in self.holes:
            s = s + ", ("
            var first_h = True
            for c in h.coords:
                if not first_h: s = s + ", "
                first_h = False
                s = s + c[0].__str__() + " " + c[1].__str__()
            s = s + ")"
        return s + ")"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        if self.shell.coords.size() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var minx = self.shell.coords[0][0]
        var miny = self.shell.coords[0][1]
        var maxx = minx
        var maxy = miny
        for c in self.shell.coords:
            if c[0] < minx: minx = c[0]
            if c[0] > maxx: maxx = c[0]
            if c[1] < miny: miny = c[1]
            if c[1] > maxy: maxy = c[1]
        return (minx, miny, maxx, maxy)


struct GeometryCollection(Geometry):
    var geoms: List[Geometry]

    fn __init__(self, geoms: List[Geometry]):
        self.geoms = geoms
        super().__init__(GeometryType.GEOMETRYCOLLECTION)

    fn is_empty(self) -> Bool:
        return self.geoms.size() == 0

    fn to_wkt(self) -> String:
        var s = "GEOMETRYCOLLECTION ("
        var first = True
        for g in self.geoms:
            if not first: s = s + ", "
            first = False
            s = s + g.to_wkt()
        return s + ")"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        if self.geoms.size() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var inited = False
        var minx = 0.0
        var miny = 0.0
        var maxx = 0.0
        var maxy = 0.0
        for g in self.geoms:
            var b = g.bounds()
            if not inited:
                minx = b[0]; miny = b[1]; maxx = b[2]; maxy = b[3]
                inited = True
            else:
                if b[0] < minx: minx = b[0]
                if b[2] > maxx: maxx = b[2]
                if b[1] < miny: miny = b[1]
                if b[3] > maxy: maxy = b[3]
        return (minx, miny, maxx, maxy)


struct MultiPoint(Geometry):
    var points: List[Point]
    fn __init__(self, points: List[Point]):
        self.points = points
        super().__init__(GeometryType.MULTIPOINT)

    fn to_wkt(self) -> String:
        var s = "MULTIPOINT ("
        var first = True
        for p in self.points:
            if not first: s = s + ", "
            first = False
            s = s + p.x.__str__() + " " + p.y.__str__()
        return s + ")"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        if self.points.size() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var minx = self.points[0].x
        var miny = self.points[0].y
        var maxx = minx
        var maxy = miny
        for p in self.points:
            if p.x < minx: minx = p.x
            if p.x > maxx: maxx = p.x
            if p.y < miny: miny = p.y
            if p.y > maxy: maxy = p.y
        return (minx, miny, maxx, maxy)


struct MultiLineString(Geometry):
    var lines: List[LineString]
    fn __init__(self, lines: List[LineString]):
        self.lines = lines
        super().__init__(GeometryType.MULTILINESTRING)

    fn to_wkt(self) -> String:
        var s = "MULTILINESTRING ("
        var first = True
        for ln in self.lines:
            if not first: s = s + ", "
            first = False
            s = s + "("
            var firstc = True
            for c in ln.coords:
                if not firstc: s = s + ", "
                firstc = False
                s = s + c[0].__str__() + " " + c[1].__str__()
            s = s + ")"
        return s + ")"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        if self.lines.size() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var inited = False
        var minx = 0.0
        var miny = 0.0
        var maxx = 0.0
        var maxy = 0.0
        for ln in self.lines:
            var b = ln.bounds()
            if not inited:
                minx = b[0]; miny = b[1]; maxx = b[2]; maxy = b[3]
                inited = True
            else:
                if b[0] < minx: minx = b[0]
                if b[2] > maxx: maxx = b[2]
                if b[1] < miny: miny = b[1]
                if b[3] > maxy: maxy = b[3]
        return (minx, miny, maxx, maxy)


struct MultiPolygon(Geometry):
    var polys: List[Polygon]
    fn __init__(self, polys: List[Polygon]):
        self.polys = polys
        super().__init__(GeometryType.MULTIPOLYGON)

    fn to_wkt(self) -> String:
        var s = "MULTIPOLYGON ("
        var firstp = True
        for poly in self.polys:
            if not firstp: s = s + ", "
            firstp = False
            s = s + "(("
            var firstc = True
            for c in poly.shell.coords:
                if not firstc: s = s + ", "
                firstc = False
                s = s + c[0].__str__() + " " + c[1].__str__()
            s = s + ")"
            for h in poly.holes:
                s = s + ", ("
                var firsth = True
                for c in h.coords:
                    if not firsth: s = s + ", "
                    firsth = False
                    s = s + c[0].__str__() + " " + c[1].__str__()
                s = s + ")"
            s = s + ")"
        return s + ")"

    fn bounds(self) -> (Float64, Float64, Float64, Float64):
        if self.polys.size() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var inited = False
        var minx = 0.0
        var miny = 0.0
        var maxx = 0.0
        var maxy = 0.0
        for poly in self.polys:
            var b = poly.bounds()
            if not inited:
                minx = b[0]; miny = b[1]; maxx = b[2]; maxy = b[3]
                inited = True
            else:
                if b[0] < minx: minx = b[0]
                if b[2] > maxx: maxx = b[2]
                if b[1] < miny: miny = b[1]
                if b[3] > maxy: maxy = b[3]
        return (minx, miny, maxx, maxy)
