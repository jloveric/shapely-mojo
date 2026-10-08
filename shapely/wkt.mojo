from shapely._geometry import Geometry
from shapely.geometry import Point, LinearRing, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon


def to_wkt(geom: Geometry) -> String:
    return geom.to_wkt()


def _wkt_strip(s: String) -> String:
    return String(s.strip())


def _wkt_slice(s: String, start: Int, end: Int) -> String:
    return String(s[byte=start:end])


def _starts_with(s: String, prefix: String) -> Bool:
    var n = prefix.byte_length()
    if n > s.byte_length():
        return False
    return _wkt_slice(s, 0, n) == prefix


def _ends_with(s: String, suffix: String) -> Bool:
    var n = suffix.byte_length()
    if n > s.byte_length():
        return False
    return _wkt_slice(s, s.byte_length() - n, s.byte_length()) == suffix


def _parse_float(slice: StringSlice) raises -> Float64:
    return Float64(String(slice))


def _parse_point(body: String) raises -> Point:
    # expects "x y" or "x y z"; only use x y
    var parts = body.split(" ")
    if len(parts) < 2:
        return Point(0.0, 0.0)
    return Point(_parse_float(parts[0]), _parse_float(parts[1]))


def _parse_ring(body: String) raises -> LinearRing:
    # expects "x y, x y, ..."
    var coords = List[Tuple[Float64, Float64]]()
    for seg in body.split(","):
        var trimmed = seg.strip()
        var xy = trimmed.split(" ")
        if len(xy) >= 2:
            coords.append((_parse_float(xy[0]), _parse_float(xy[1])))
    return LinearRing(coords)


def from_wkt(wkt: String) raises -> Geometry:
    var s = _wkt_strip(wkt)
    if _starts_with(s, "POINT"):
        var l = s.find("(")
        var r = s.rfind(")")
        if l >= 0 and r > l:
            var body = _wkt_strip(_wkt_slice(s, l + 1, r))
            return Geometry(_parse_point(body))
    if _starts_with(s, "LINEARRING"):
        var l = s.find("(")
        var r = s.rfind(")")
        if l >= 0 and r > l:
            var body = _wkt_strip(_wkt_slice(s, l + 1, r))
            var ring = _parse_ring(body)
            return Geometry(LineString(ring.coords.copy()))
    if _starts_with(s, "LINESTRING"):
        var l = s.find("(")
        var r = s.rfind(")")
        if l >= 0 and r > l:
            var body = _wkt_strip(_wkt_slice(s, l + 1, r))
            var coords = List[Tuple[Float64, Float64]]()
            for seg in body.split(","):
                var xy = seg.strip().split(" ")
                if len(xy) >= 2:
                    coords.append((_parse_float(xy[0]), _parse_float(xy[1])))
            return Geometry(LineString(coords))
    if _starts_with(s, "MULTILINESTRING"):
        var l = s.find("(")
        var r = s.rfind(")")
        if l >= 0 and r > l:
            var inner = _wkt_strip(_wkt_slice(s, l + 1, r))
            var grouped = inner.replace("), (", ")|("
                ).replace("),(", ")|("
                ).replace("( (", "(("
                ).replace(" )", ")")
            var lines = List[LineString]()
            for g in grouped.split("|"):
                var gg = _wkt_strip(String(g))
                var ring_body = gg
                if _starts_with(gg, "(") and _ends_with(gg, ")"):
                    ring_body = _wkt_slice(gg, 1, gg.byte_length() - 1)
                var coords = List[Tuple[Float64, Float64]]()
                for seg in ring_body.split(","):
                    var xy = seg.strip().split(" ")
                    if len(xy) >= 2:
                        coords.append((_parse_float(xy[0]), _parse_float(xy[1])))
                lines.append(LineString(coords))
            return Geometry(MultiLineString(lines))
    if _starts_with(s, "MULTIPOINT"):
        var l = s.find("(")
        var r = s.rfind(")")
        if l >= 0 and r > l:
            var body = _wkt_strip(_wkt_slice(s, l + 1, r))
            var pts = List[Point]()
            if body.find("(") >= 0:
                for chunk in body.split(")"):
                    var inner = chunk.replace("(", "").replace(",", " ").strip()
                    if String(inner).byte_length() == 0: continue
                    var xy = inner.split(" ")
                    if len(xy) >= 2:
                        pts.append(Point(_parse_float(xy[0]), _parse_float(xy[1])))
            else:
                for seg in body.split(","):
                    var xy = seg.strip().split(" ")
                    if len(xy) >= 2:
                        pts.append(Point(_parse_float(xy[0]), _parse_float(xy[1])))
            return Geometry(MultiPoint(pts))
    if _starts_with(s, "POLYGON"):
        var l = s.find("(")
        var r = s.rfind(")")
        if l >= 0 and r > l:
            var inner = _wkt_strip(_wkt_slice(s, l + 1, r))
            # split rings by '), (' patterns
            var grouped = inner.replace("), (", ")|("
                ).replace("),(", ")|("
                ).strip()
            var rings = List[LinearRing]()
            for g in grouped.split("|"):
                var gg = _wkt_strip(String(g))
                var ring_body = gg
                if _starts_with(gg, "(") and _ends_with(gg, ")"):
                    ring_body = _wkt_slice(gg, 1, gg.byte_length() - 1)
                rings.append(_parse_ring(ring_body))
            if len(rings) == 0:
                return Geometry(Polygon(LinearRing(List[Tuple[Float64, Float64]]())))
            var shell = rings[0].copy()
            var holes = List[LinearRing]()
            var i = 1
            while i < len(rings):
                holes.append(rings[i].copy())
                i += 1
            return Geometry(Polygon(shell, holes))
    if _starts_with(s, "MULTIPOLYGON"):
        var l = s.find("(")
        var r = s.rfind(")")
        if l >= 0 and r > l:
            var inner = _wkt_strip(_wkt_slice(s, l + 1, r))
            # split polygons by ')), ((' boundaries
            var grouped = inner.replace(")), ((", "))|(("
                ).replace(")),((", "))|(("
                ).strip()
            var polys = List[Polygon]()
            for pg in grouped.split("|"):
                var pgs = _wkt_strip(String(pg))
                var body = pgs
                if _starts_with(pgs, "((") and _ends_with(pgs, "))"):
                    body = _wkt_slice(pgs, 1, pgs.byte_length() - 1)
                # Now body contains '(ring),(ring),...'
                var rings_grouped = body.replace("), (", ")|("
                    ).replace("),(", ")|("
                    ).strip()
                var rings = List[LinearRing]()
                for g in rings_grouped.split("|"):
                    var gg = _wkt_strip(String(g))
                    var ring_body = gg
                    if _starts_with(gg, "(") and _ends_with(gg, ")"):
                        ring_body = _wkt_slice(gg, 1, gg.byte_length() - 1)
                    rings.append(_parse_ring(ring_body))
                if len(rings) > 0:
                    var shell = rings[0].copy()
                    var holes = List[LinearRing]()
                    var i = 1
                    while i < len(rings):
                        holes.append(rings[i].copy())
                        i += 1
                    polys.append(Polygon(shell, holes))
            return Geometry(MultiPolygon(polys))
    # default fallback
    return Geometry(Polygon(LinearRing(List[Tuple[Float64, Float64]]())))
