from shapely._geometry import Geometry


fn sqrt_f64(x: Float64) -> Float64:
    if x <= 0.0:
        return 0.0
    var r = x
    var i = 0
    while i < 12:
        r = 0.5 * (r + x / r)
        i += 1
    return r


struct Point(Copyable, Movable):
    var x: Float64
    var y: Float64
    var has_z: Bool
    var z: Float64

    fn __init__(out self, x: Float64, y: Float64):
        self.x = x
        self.y = y
        self.has_z = False
        self.z = 0.0

    fn with_z(self, z: Float64) -> Point:
        var p = Point(self.x, self.y)
        p.has_z = True
        p.z = z
        return p

    fn is_empty(self) -> Bool:
        return False

    fn to_wkt(self) -> String:
        if self.has_z:
            return (
                "POINT Z ("
                + self.x.__str__()
                + " "
                + self.y.__str__()
                + " "
                + self.z.__str__()
                + ")"
            )
        return "POINT (" + self.x.__str__() + " " + self.y.__str__() + ")"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        return (self.x, self.y, self.x, self.y)


struct LineString(Copyable, Movable):
    var coords: List[Tuple[Float64, Float64]]

    fn __init__(out self, coords: List[Tuple[Float64, Float64]]):
        self.coords = coords.copy()

    fn is_empty(self) -> Bool:
        return self.coords.__len__() == 0

    fn to_wkt(self) -> String:
        var s = "LINESTRING ("
        var first = True
        for c in self.coords:
            if not first:
                s = s + ", "
            first = False
            s = s + c[0].__str__() + " " + c[1].__str__()
        return s + ")"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        if self.coords.__len__() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var minx = self.coords[0][0]
        var miny = self.coords[0][1]
        var maxx = minx
        var maxy = miny
        for c in self.coords:
            if c[0] < minx:
                minx = c[0]
            if c[0] > maxx:
                maxx = c[0]
            if c[1] < miny:
                miny = c[1]
            if c[1] > maxy:
                maxy = c[1]
        return (minx, miny, maxx, maxy)

    fn length(self) -> Float64:
        if self.coords.__len__() <= 1:
            return 0.0
        var total = 0.0
        for i in range(0, self.coords.__len__() - 1):
            var a = self.coords[i]
            var b = self.coords[i + 1]
            var dx = b[0] - a[0]
            var dy = b[1] - a[1]
            total += sqrt_f64(dx * dx + dy * dy)
        return total


struct LinearRing(Copyable, Movable):
    var coords: List[Tuple[Float64, Float64]]

    fn __init__(out self, coords: List[Tuple[Float64, Float64]]):
        self.coords = coords.copy()

    fn is_empty(self) -> Bool:
        return self.coords.__len__() == 0

    fn to_wkt(self) -> String:
        var s = "LINEARRING ("
        var first = True
        for c in self.coords:
            if not first:
                s = s + ", "
            first = False
            s = s + c[0].__str__() + " " + c[1].__str__()
        return s + ")"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        if self.coords.__len__() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var minx = self.coords[0][0]
        var miny = self.coords[0][1]
        var maxx = minx
        var maxy = miny
        for c in self.coords:
            if c[0] < minx:
                minx = c[0]
            if c[0] > maxx:
                maxx = c[0]
            if c[1] < miny:
                miny = c[1]
            if c[1] > maxy:
                maxy = c[1]
        return (minx, miny, maxx, maxy)

    fn length(self) -> Float64:
        if self.coords.__len__() <= 1:
            return 0.0
        var total = 0.0
        for i in range(0, self.coords.__len__() - 1):
            var a = self.coords[i]
            var b = self.coords[i + 1]
            var dx = b[0] - a[0]
            var dy = b[1] - a[1]
            total += sqrt_f64(dx * dx + dy * dy)
        return total


struct Polygon(Copyable, Movable):
    var shell: LinearRing
    var holes: List[LinearRing]

    fn __init__(
        out self,
        shell: LinearRing,
        holes: List[LinearRing] = List[LinearRing](),
    ):
        self.shell = shell.copy()
        self.holes = holes.copy()

    fn is_empty(self) -> Bool:
        return self.shell.coords.__len__() == 0

    fn to_wkt(self) -> String:
        var s = "POLYGON ("
        # shell
        s = s + "("
        var first = True
        for c in self.shell.coords:
            if not first:
                s = s + ", "
            first = False
            s = s + c[0].__str__() + " " + c[1].__str__()
        s = s + ")"
        # holes
        for h in self.holes:
            s = s + ", ("
            var first_h = True
            for c in h.coords:
                if not first_h:
                    s = s + ", "
                first_h = False
                s = s + c[0].__str__() + " " + c[1].__str__()
            s = s + ")"
        return s + ")"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        if self.shell.coords.__len__() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var minx = self.shell.coords[0][0]
        var miny = self.shell.coords[0][1]
        var maxx = minx
        var maxy = miny
        for c in self.shell.coords:
            if c[0] < minx:
                minx = c[0]
            if c[0] > maxx:
                maxx = c[0]
            if c[1] < miny:
                miny = c[1]
            if c[1] > maxy:
                maxy = c[1]
        return (minx, miny, maxx, maxy)

    fn area(self) -> Float64:
        if self.shell.coords.__len__() < 3:
            return 0.0
        var shell_sum = 0.0
        for i in range(0, self.shell.coords.__len__() - 1):
            var a = self.shell.coords[i]
            var b = self.shell.coords[i + 1]
            shell_sum += a[0] * b[1] - a[1] * b[0]
        var holes_sum = 0.0
        for h in self.holes:
            if h.coords.__len__() >= 3:
                var hs = 0.0
                for i in range(0, h.coords.__len__() - 1):
                    var a = h.coords[i]
                    var b = h.coords[i + 1]
                    hs += a[0] * b[1] - a[1] * b[0]
                holes_sum += 0.5 * abs(hs)
        var total = 0.5 * (shell_sum - holes_sum)
        if total < 0.0:
            return -total
        return total

    fn length(self) -> Float64:
        var s = self.shell.length()
        for h in self.holes:
            s += h.length()
        return s


struct GeometryCollection(Copyable, Movable):
    var geoms: List[Geometry]

    fn __init__(out self, geoms: List[Geometry]):
        self.geoms = geoms.copy()

    fn is_empty(self) -> Bool:
        return self.geoms.__len__() == 0

    fn to_wkt(self) -> String:
        var s = "GEOMETRYCOLLECTION ("
        var first = True
        for g in self.geoms:
            if not first:
                s = s + ", "
            first = False
            s = s + g.to_wkt()
        return s + ")"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        if self.geoms.__len__() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var inited = False
        var minx = 0.0
        var miny = 0.0
        var maxx = 0.0
        var maxy = 0.0
        for g in self.geoms:
            var b = g.bounds()
            if not inited:
                minx = b[0]
                miny = b[1]
                maxx = b[2]
                maxy = b[3]
                inited = True
            else:
                if b[0] < minx:
                    minx = b[0]
                if b[2] > maxx:
                    maxx = b[2]
                if b[1] < miny:
                    miny = b[1]
                if b[3] > maxy:
                    maxy = b[3]
        return (minx, miny, maxx, maxy)


struct MultiPoint(Copyable, Movable):
    var points: List[Point]

    fn __init__(out self, points: List[Point]):
        self.points = points.copy()

    fn is_empty(self) -> Bool:
        return self.points.__len__() == 0

    fn to_wkt(self) -> String:
        var s = "MULTIPOINT ("
        var first = True
        for p in self.points:
            if not first:
                s = s + ", "
            first = False
            s = s + p.x.__str__() + " " + p.y.__str__()
        return s + ")"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        if self.points.__len__() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var minx = self.points[0].x
        var miny = self.points[0].y
        var maxx = minx
        var maxy = miny
        for p in self.points:
            if p.x < minx:
                minx = p.x
            if p.x > maxx:
                maxx = p.x
            if p.y < miny:
                miny = p.y
            if p.y > maxy:
                maxy = p.y
        return (minx, miny, maxx, maxy)


struct MultiLineString(Copyable, Movable):
    var lines: List[LineString]

    fn __init__(out self, lines: List[LineString]):
        self.lines = lines.copy()

    fn is_empty(self) -> Bool:
        return self.lines.__len__() == 0

    fn to_wkt(self) -> String:
        var s = "MULTILINESTRING ("
        var first = True
        for ln in self.lines:
            if not first:
                s = s + ", "
            first = False
            s = s + "("
            var firstc = True
            for c in ln.coords:
                if not firstc:
                    s = s + ", "
                firstc = False
                s = s + c[0].__str__() + " " + c[1].__str__()
            s = s + ")"
        return s + ")"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        if self.lines.__len__() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var inited = False
        var minx = 0.0
        var miny = 0.0
        var maxx = 0.0
        var maxy = 0.0
        for ln in self.lines:
            var b = ln.bounds()
            if not inited:
                minx = b[0]
                miny = b[1]
                maxx = b[2]
                maxy = b[3]
                inited = True
            else:
                if b[0] < minx:
                    minx = b[0]
                if b[2] > maxx:
                    maxx = b[2]
                if b[1] < miny:
                    miny = b[1]
                if b[3] > maxy:
                    maxy = b[3]
        return (minx, miny, maxx, maxy)


    fn length(self) -> Float64:
        var s = 0.0
        for ln in self.lines:
            s += ln.length()
        return s


struct MultiPolygon(Copyable, Movable):
    var polys: List[Polygon]

    fn __init__(out self, polys: List[Polygon]):
        self.polys = polys.copy()

    fn is_empty(self) -> Bool:
        return self.polys.__len__() == 0

    fn to_wkt(self) -> String:
        var s = "MULTIPOLYGON ("
        var firstp = True
        for poly in self.polys:
            if not firstp:
                s = s + ", "
            firstp = False
            s = s + "(("
            var firstc = True
            for c in poly.shell.coords:
                if not firstc:
                    s = s + ", "
                firstc = False
                s = s + c[0].__str__() + " " + c[1].__str__()
            s = s + ")"
            for h in poly.holes:
                s = s + ", ("
                var firsth = True
                for c in h.coords:
                    if not firsth:
                        s = s + ", "
                    firsth = False
                    s = s + c[0].__str__() + " " + c[1].__str__()
                s = s + ")"
            s = s + ")"
        return s + ")"

    fn bounds(self) -> Tuple[Float64, Float64, Float64, Float64]:
        if self.polys.__len__() == 0:
            return (0.0, 0.0, 0.0, 0.0)
        var inited = False
        var minx = 0.0
        var miny = 0.0
        var maxx = 0.0
        var maxy = 0.0
        for poly in self.polys:
            var b = poly.bounds()
            if not inited:
                minx = b[0]
                miny = b[1]
                maxx = b[2]
                maxy = b[3]
                inited = True
            else:
                if b[0] < minx:
                    minx = b[0]
                if b[2] > maxx:
                    maxx = b[2]
                if b[1] < miny:
                    miny = b[1]
                if b[3] > maxy:
                    maxy = b[3]
        return (minx, miny, maxx, maxy)


    fn area(self) -> Float64:
        var s = 0.0
        for p in self.polys:
            s += p.area()
        return s


    fn length(self) -> Float64:
        var s = 0.0
        for p in self.polys:
            s += p.length()
        return s
