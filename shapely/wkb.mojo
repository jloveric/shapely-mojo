from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPoint, MultiLineString, MultiPolygon, GeometryCollection


fn _write_u8(out: List[UInt8], v: UInt8) -> List[UInt8]:
    var o = out
    o.append(v)
    return o


fn _write_u32_le(out: List[UInt8], v: UInt32) -> List[UInt8]:
    var o = out
    o.append((v & 0xFF) as UInt8)
    o.append(((v >> 8) & 0xFF) as UInt8)
    o.append(((v >> 16) & 0xFF) as UInt8)
    o.append(((v >> 24) & 0xFF) as UInt8)
    return o


fn _write_u64_le(out: List[UInt8], v: UInt64) -> List[UInt8]:
    var o = out
    o.append((v & 0xFF) as UInt8)
    o.append(((v >> 8) & 0xFF) as UInt8)
    o.append(((v >> 16) & 0xFF) as UInt8)
    o.append(((v >> 24) & 0xFF) as UInt8)
    o.append(((v >> 32) & 0xFF) as UInt8)
    o.append(((v >> 40) & 0xFF) as UInt8)
    o.append(((v >> 48) & 0xFF) as UInt8)
    o.append(((v >> 56) & 0xFF) as UInt8)
    return o


fn _f64_to_le_bytes(out: List[UInt8], x: Float64) -> List[UInt8]:
    let bits = unsafe_bitcast[UInt64](x)
    return _write_u64_le(out, bits)


fn _read_u8(buf: List[UInt8], pos: Int) -> (UInt8, Int):
    let v = buf[pos]
    return (v, pos + 1)


fn _read_u32(buf: List[UInt8], pos: Int, little: Bool) -> (UInt32, Int):
    var b0 = buf[pos] as UInt32
    var b1 = buf[pos + 1] as UInt32
    var b2 = buf[pos + 2] as UInt32
    var b3 = buf[pos + 3] as UInt32
    let np = pos + 4
    if little:
        return (b0 | (b1 << 8) | (b2 << 16) | (b3 << 24), np)
    else:
        return (b3 | (b2 << 8) | (b1 << 16) | (b0 << 24), np)


fn _read_u64(buf: List[UInt8], pos: Int, little: Bool) -> (UInt64, Int):
    var v: UInt64 = 0
    if little:
        var i = 0
        while i < 8:
            v = v | ((buf[pos + i] as UInt64) << (i * 8))
            i += 1
    else:
        var i2 = 0
        while i2 < 8:
            v = (v << 8) | (buf[pos + i2] as UInt64)
            i2 += 1
    return (v, pos + 8)


fn _read_f64(buf: List[UInt8], pos: Int, little: Bool) -> (Float64, Int):
    let (bits, np) = _read_u64(buf, pos, little)
    return (unsafe_bitcast[Float64](bits), np)


fn _to_wkb_point(p: Point) -> List[UInt8]:
    var out = List[UInt8]()
    out = _write_u8(out, 1)  # little endian
    out = _write_u32_le(out, 1)  # Point type
    out = _f64_to_le_bytes(out, p.x)
    out = _f64_to_le_bytes(out, p.y)
    return out


fn _to_wkb_linestring(ls: LineString) -> List[UInt8]:
    var out = List[UInt8]()
    out = _write_u8(out, 1)
    out = _write_u32_le(out, 2)
    out = _write_u32_le(out, ls.coords.size() as UInt32)
    for c in ls.coords:
        out = _f64_to_le_bytes(out, c[0])
        out = _f64_to_le_bytes(out, c[1])
    return out


fn _to_wkb_ring(out: List[UInt8], ring: LinearRing) -> List[UInt8]:
    var o = out
    o = _write_u32_le(o, ring.coords.size() as UInt32)
    for c in ring.coords:
        o = _f64_to_le_bytes(o, c[0])
        o = _f64_to_le_bytes(o, c[1])
    return o


fn _to_wkb_polygon(poly: Polygon) -> List[UInt8]:
    var out = List[UInt8]()
    out = _write_u8(out, 1)
    out = _write_u32_le(out, 3)
    let nrings = (1 + poly.holes.size()) as UInt32
    out = _write_u32_le(out, nrings)
    out = _to_wkb_ring(out, poly.shell)
    for h in poly.holes:
        out = _to_wkb_ring(out, h)
    return out


fn _to_wkb_multipoint(mp: MultiPoint) -> List[UInt8]:
    var out = List[UInt8]()
    out = _write_u8(out, 1)
    out = _write_u32_le(out, 4)
    out = _write_u32_le(out, mp.points.size() as UInt32)
    for p in mp.points:
        let sub = _to_wkb_point(p)
        for b in sub: out.append(b)
    return out


fn _to_wkb_multilinestring(mls: MultiLineString) -> List[UInt8]:
    var out = List[UInt8]()
    out = _write_u8(out, 1)
    out = _write_u32_le(out, 5)
    out = _write_u32_le(out, mls.lines.size() as UInt32)
    for ln in mls.lines:
        let sub = _to_wkb_linestring(ln)
        for b in sub: out.append(b)
    return out


fn _to_wkb_multipolygon(mp: MultiPolygon) -> List[UInt8]:
    var out = List[UInt8]()
    out = _write_u8(out, 1)
    out = _write_u32_le(out, 6)
    out = _write_u32_le(out, mp.polys.size() as UInt32)
    for p in mp.polys:
        let sub = _to_wkb_polygon(p)
        for b in sub: out.append(b)
    return out


fn to_wkb(g: Geometry) -> List[UInt8]:
    let t = g.__type_name__()
    if t == "Point":
        return _to_wkb_point(unsafe_bitcast[Point](g))
    if t == "LineString":
        return _to_wkb_linestring(unsafe_bitcast[LineString](g))
    if t == "Polygon":
        return _to_wkb_polygon(unsafe_bitcast[Polygon](g))
    if t == "MultiPoint":
        return _to_wkb_multipoint(unsafe_bitcast[MultiPoint](g))
    if t == "MultiLineString":
        return _to_wkb_multilinestring(unsafe_bitcast[MultiLineString](g))
    if t == "MultiPolygon":
        return _to_wkb_multipolygon(unsafe_bitcast[MultiPolygon](g))
    # GeometryCollection: write collection of nested
    if t == "GeometryCollection":
        let gc = unsafe_bitcast[GeometryCollection](g)
        var out = List[UInt8]()
        _write_u8(out, 1)
        _write_u32_le(out, 7)
        _write_u32_le(out, gc.geoms.size() as UInt32)
        for gg in gc.geoms:
            let sub = to_wkb(gg)
            for b in sub: out.append(b)
        return out
    return List[UInt8]()


fn _from_wkb_point(buf: List[UInt8], pos: Int, little: Bool) -> (Point, Int):
    let (x, p1) = _read_f64(buf, pos, little)
    let (y, p2) = _read_f64(buf, p1, little)
    return (Point(x, y), p2)


fn _from_wkb_linestring(buf: List[UInt8], pos: Int, little: Bool) -> (LineString, Int):
    let (nv, p0) = _read_u32(buf, pos, little)
    let n = nv as Int
    var coords = List[Tuple[Float64, Float64]]()
    var i = 0
    var p = p0
    while i < n:
        let (x, p1) = _read_f64(buf, p, little)
        let (y, p2) = _read_f64(buf, p1, little)
        coords.append((x, y))
        i += 1
        p = p2
    return (LineString(coords), p)


fn _from_wkb_ring(buf: List[UInt8], pos: Int, little: Bool) -> (LinearRing, Int):
    let (nv, p0) = _read_u32(buf, pos, little)
    let n = nv as Int
    var coords = List[Tuple[Float64, Float64]]()
    var i = 0
    var p = p0
    while i < n:
        let (x, p1) = _read_f64(buf, p, little)
        let (y, p2) = _read_f64(buf, p1, little)
        coords.append((x, y))
        i += 1
        p = p2
    return (LinearRing(coords), p)


fn _from_wkb_polygon(buf: List[UInt8], pos: Int, little: Bool) -> (Polygon, Int):
    let (nrv, p0) = _read_u32(buf, pos, little)
    let nrings = nrv as Int
    if nrings == 0:
        return (Polygon(LinearRing(List[Tuple[Float64, Float64]]())), p0)
    let (shell, p1) = _from_wkb_ring(buf, p0, little)
    var holes = List[LinearRing]()
    var i = 1
    var p = p1
    while i < nrings:
        let (hr, p2) = _from_wkb_ring(buf, p, little)
        holes.append(hr)
        i += 1
        p = p2
    return (Polygon(shell, holes), p)


fn _from_wkb_multipoint(buf: List[UInt8], pos: Int, little: Bool) -> (MultiPoint, Int):
    let (nv, p0) = _read_u32(buf, pos, little)
    let n = nv as Int
    var pts = List[Point]()
    var i = 0
    var p = p0
    while i < n:
        let (g, p2) = _read_geom(buf, p)
        if g.__type_name__() == "Point":
            pts.append(unsafe_bitcast[Point](g))
        i += 1
        p = p2
    return (MultiPoint(pts), p)


fn _from_wkb_multilinestring(buf: List[UInt8], pos: Int, little: Bool) -> (MultiLineString, Int):
    let (nv, p0) = _read_u32(buf, pos, little)
    let n = nv as Int
    var lines = List[LineString]()
    var i = 0
    var p = p0
    while i < n:
        let (g, p2) = _read_geom(buf, p)
        if g.__type_name__() == "LineString":
            lines.append(unsafe_bitcast[LineString](g))
        i += 1
        p = p2
    return (MultiLineString(lines), p)


fn _from_wkb_multipolygon(buf: List[UInt8], pos: Int, little: Bool) -> (MultiPolygon, Int):
    let (nv, p0) = _read_u32(buf, pos, little)
    let n = nv as Int
    var polys = List[Polygon]()
    var i = 0
    var p = p0
    while i < n:
        let (g, p2) = _read_geom(buf, p)
        if g.__type_name__() == "Polygon":
            polys.append(unsafe_bitcast[Polygon](g))
        i += 1
        p = p2
    return (MultiPolygon(polys), p)


fn _read_geom(buf: List[UInt8], pos: Int) -> (Geometry, Int):
    if pos + 5 > buf.size():
        return (Polygon(LinearRing(List[Tuple[Float64, Float64]]())), pos)
    let (e, p0) = _read_u8(buf, pos)
    let little = (e == 1)
    let (t, p1) = _read_u32(buf, p0, little)
    if t == 1:
        return _from_wkb_point(buf, p1, little)
    if t == 2:
        return _from_wkb_linestring(buf, p1, little)
    if t == 3:
        return _from_wkb_polygon(buf, p1, little)
    if t == 4:
        return _from_wkb_multipoint(buf, p1, little)
    if t == 5:
        return _from_wkb_multilinestring(buf, p1, little)
    if t == 6:
        return _from_wkb_multipolygon(buf, p1, little)
    if t == 7:
        let (nv, p2) = _read_u32(buf, p1, little)
        let n = nv as Int
        var geoms = List[Geometry]()
        var i = 0
        var p = p2
        while i < n:
            let (g, pn) = _read_geom(buf, p)
            geoms.append(g)
            i += 1
            p = pn
        return (GeometryCollection(geoms), p)
    return (Polygon(LinearRing(List[Tuple[Float64, Float64]]())), p1)


fn from_wkb(buf: List[UInt8]) -> Geometry:
    var pos = 0
    let (g, _np) = _read_geom(buf, pos)
    return g
