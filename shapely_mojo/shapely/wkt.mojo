from shapely._geometry import Geometry
from shapely.geometry import Point, LinearRing, LineString, Polygon, MultiPoint, MultiLineString, MultiPolygon


fn to_wkt(geom: Geometry) -> String:
    return geom.to_wkt()


fn _parse_point(body: String) -> Point:
    # expects "x y" or "x y z"; only use x y
    let parts = body.split(" ")
    if parts.size() < 2:
        return Point(0.0, 0.0)
    return Point(parts[0].to_float64(), parts[1].to_float64())


fn _parse_ring(body: String) -> LinearRing:
    # expects "x y, x y, ..."
    var coords = List[Tuple[Float64, Float64]]()
    for seg in body.split(","):
        let trimmed = seg.strip()
        let xy = trimmed.split(" ")
        if xy.size() >= 2:
            coords.append((xy[0].to_float64(), xy[1].to_float64()))
    return LinearRing(coords)


fn from_wkt(wkt: String) -> Geometry:
    let s = wkt.strip()
    if s.starts_with("POINT"):
        let l = s.find("(")
        let r = s.rfind(")")
        if l >= 0 and r > l:
            let body = s.slice(l + 1, r).strip()
            return _parse_point(body)
    if s.starts_with("LINEARRING"):
        let l = s.find("(")
        let r = s.rfind(")")
        if l >= 0 and r > l:
            let body = s.slice(l + 1, r).strip()
            return _parse_ring(body)
    if s.starts_with("LINESTRING"):
        let l = s.find("(")
        let r = s.rfind(")")
        if l >= 0 and r > l:
            let body = s.slice(l + 1, r).strip()
            var coords = List[Tuple[Float64, Float64]]()
            for seg in body.split(","):
                let xy = seg.strip().split(" ")
                if xy.size() >= 2:
                    coords.append((xy[0].to_float64(), xy[1].to_float64()))
            return LineString(coords)
    if s.starts_with("MULTILINESTRING"):
        let l = s.find("(")
        let r = s.rfind(")")
        if l >= 0 and r > l:
            let inner = s.slice(l + 1, r).strip()
            let grouped = inner.replace("), (", ")|("
                ).replace("),(", ")|("
                ).replace("( (", "(("
                ).replace(" )", ")")
            var lines = List[LineString]()
            for g in grouped.split("|"):
                let gg = g.strip()
                var ring_body = gg
                if gg.starts_with("(") and gg.ends_with(")"):
                    ring_body = gg.slice(1, gg.size() - 1)
                var coords = List[Tuple[Float64, Float64]]()
                for seg in ring_body.split(","):
                    let xy = seg.strip().split(" ")
                    if xy.size() >= 2:
                        coords.append((xy[0].to_float64(), xy[1].to_float64()))
                lines.append(LineString(coords))
            return MultiLineString(lines)
    if s.starts_with("MULTIPOINT"):
        let l = s.find("(")
        let r = s.rfind(")")
        if l >= 0 and r > l:
            let body = s.slice(l + 1, r).strip()
            var pts = List[Point]()
            if body.find("(") >= 0:
                for chunk in body.split(")"):
                    let inner = chunk.replace("(", "").replace(",", " ").strip()
                    if inner.size() == 0: continue
                    let xy = inner.split(" ")
                    if xy.size() >= 2:
                        pts.append(Point(xy[0].to_float64(), xy[1].to_float64()))
            else:
                for seg in body.split(","):
                    let xy = seg.strip().split(" ")
                    if xy.size() >= 2:
                        pts.append(Point(xy[0].to_float64(), xy[1].to_float64()))
            return MultiPoint(pts)
    if s.starts_with("POLYGON"):
        let l = s.find("(")
        let r = s.rfind(")")
        if l >= 0 and r > l:
            let inner = s.slice(l + 1, r).strip()
            # split rings by '), (' patterns
            let grouped = inner.replace("), (", ")|("
                ).replace("),(", ")|("
                ).strip()
            var rings = List[LinearRing]()
            for g in grouped.split("|"):
                let gg = g.strip()
                var ring_body = gg
                if gg.starts_with("(") and gg.ends_with(")"):
                    ring_body = gg.slice(1, gg.size() - 1)
                rings.append(_parse_ring(ring_body))
            if rings.size() == 0:
                return Polygon(LinearRing(List[Tuple[Float64, Float64]]()))
            let shell = rings[0]
            var holes = List[LinearRing]()
            var i = 1
            while i < rings.size():
                holes.append(rings[i])
                i += 1
            return Polygon(shell, holes)
    if s.starts_with("MULTIPOLYGON"):
        let l = s.find("(")
        let r = s.rfind(")")
        if l >= 0 and r > l:
            let inner = s.slice(l + 1, r).strip()
            # split polygons by ')), ((' boundaries
            let grouped = inner.replace(")), ((", "))|(("
                ).replace(")),((", "))|(("
                ).strip()
            var polys = List[Polygon]()
            for pg in grouped.split("|"):
                let pgs = pg.strip()
                var body = pgs
                if pgs.starts_with("((") and pgs.ends_with("))"):
                    body = pgs.slice(1, pgs.size() - 1)  # remove one paren from each side -> '(...) , (...)'
                # Now body contains '(ring),(ring),...'
                let rings_grouped = body.replace("), (", ")|("
                    ).replace("),(", ")|("
                    ).strip()
                var rings = List[LinearRing]()
                for g in rings_grouped.split("|"):
                    let gg = g.strip()
                    var ring_body = gg
                    if gg.starts_with("(") and gg.ends_with(")"):
                        ring_body = gg.slice(1, gg.size() - 1)
                    rings.append(_parse_ring(ring_body))
                if rings.size() > 0:
                    let shell = rings[0]
                    var holes = List[LinearRing]()
                    var i = 1
                    while i < rings.size():
                        holes.append(rings[i])
                        i += 1
                    polys.append(Polygon(shell, holes))
            return MultiPolygon(polys)
    # default fallback
    return Polygon(LinearRing(List[Tuple[Float64, Float64]]()))
