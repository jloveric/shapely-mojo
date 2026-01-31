from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, MultiLineString, MultiPoint, Polygon, GeometryCollection, LinearRing
from shapely.linear import line_merge as _line_merge, shared_paths as _shared_paths
from shapely.set_operations import unary_union as _unary_union, intersection as _poly_intersection
from shapely.algorithms import on_segment, segment_intersections, point_in_ring, signed_area_coords


struct PolygonizeSeg(Copyable, Movable):
    var ax: Float64
    var ay: Float64
    var bx: Float64
    var by: Float64

    fn __init__(out self, ax: Float64, ay: Float64, bx: Float64, by: Float64):
        self.ax = ax
        self.ay = ay
        self.bx = bx
        self.by = by


struct PolygonizeDEdge(Copyable, Movable):
    var src: Int32
    var dst: Int32
    var dx: Float64
    var dy: Float64
    var alive: Bool
    var used: Bool

    fn __init__(out self, src: Int32, dst: Int32, dx: Float64, dy: Float64, alive: Bool, used: Bool):
        self.src = src
        self.dst = dst
        self.dx = dx
        self.dy = dy
        self.alive = alive
        self.used = used


fn _absf(x: Float64) -> Float64:
    if x < 0.0:
        return -x
    return x


fn _polygonize_get_vid(x: Float64, y: Float64, mut verts: List[Tuple[Float64, Float64]], eps: Float64 = 1e-12) -> Int32:
    var k = 0
    while k < verts.__len__():
        var v = verts[k]
        if _absf(v[0] - x) <= eps and _absf(v[1] - y) <= eps:
            return Int32(k)
        k += 1
    verts.append((x, y))
    return Int32(verts.__len__() - 1)


fn _polygonize_ensure_adj(mut adj: List[List[Int32]], vid: Int32):
    while adj.__len__() <= Int(vid):
        adj.append(List[Int32]())


fn _polygonize_add_unique_edge(
    a: Int32,
    b: Int32,
    mut seen: List[Tuple[Int32, Int32]],
    mut out_lines: List[LineString],
    ref verts: List[Tuple[Float64, Float64]],
):
    var u0 = a
    var u1 = b
    if u1 < u0:
        var tmp = u0
        u0 = u1
        u1 = tmp
    var found = False
    var i = 0
    while i < seen.__len__():
        if seen[i][0] == u0 and seen[i][1] == u1:
            found = True
            break
        i += 1
    if not found:
        seen.append((u0, u1))
        out_lines.append(LineString([(verts[u0][0], verts[u0][1]), (verts[u1][0], verts[u1][1])]))


fn _polygonize_next_edge(ref adj: List[List[Int32]], ref edges: List[PolygonizeDEdge], at_vertex: Int32, bx: Float64, by: Float64) -> Int32:
    if Int(at_vertex) >= adj.__len__():
        return -1
    ref cand = adj[at_vertex]
    # Prefer not to immediately backtrack along the reverse edge we just came from.
    # The reverse edge is characterized by being colinear with (bx, by) and having
    # positive dot-product.
    var best_idx: Int32 = -1
    var best_cross = 1.0e308
    var best_dot = -1.0e308
    var best_idx_fallback: Int32 = -1
    var best_cross_fallback = 1.0e308
    var best_dot_fallback = -1.0e308
    var i = 0
    while i < cand.__len__():
        var ei = cand[i]
        if not edges[ei].alive or edges[ei].used:
            i += 1
            continue
        var w_x = edges[ei].dx
        var w_y = edges[ei].dy
        var cr = bx * w_y - by * w_x
        var dt = bx * w_x + by * w_y
        # track fallback (including back-edge)
        if cr < 0.0:
            # prefer smallest clockwise turn (cr closest to 0 from below)
            if best_idx_fallback == -1 or best_cross_fallback >= 0.0 or cr > best_cross_fallback or (cr == best_cross_fallback and dt > best_dot_fallback):
                best_idx_fallback = ei
                best_cross_fallback = cr
                best_dot_fallback = dt
        elif best_idx_fallback == -1 and best_cross_fallback >= 0.0:
            if cr < best_cross_fallback or (cr == best_cross_fallback and dt > best_dot_fallback):
                best_idx_fallback = ei
                best_cross_fallback = cr
                best_dot_fallback = dt

        # skip immediate backtracking if possible
        if cr == 0.0 and dt > 0.0:
            i += 1
            continue

        if cr < 0.0:
            if best_idx == -1 or best_cross >= 0.0 or cr > best_cross or (cr == best_cross and dt > best_dot):
                best_idx = ei
                best_cross = cr
                best_dot = dt
        elif best_idx == -1 and best_cross >= 0.0:
            if cr < best_cross or (cr == best_cross and dt > best_dot):
                best_idx = ei
                best_cross = cr
                best_dot = dt
        i += 1
    if best_idx != -1:
        return best_idx
    return best_idx_fallback


fn linemerge(lines) -> Geometry:
    return _line_merge(lines)


fn polygonize(lines) -> GeometryCollection:
    # Collect all LineStrings from input
    var lns = List[LineString]()
    if lines.is_linestring():
        lns.append(lines.as_linestring())
    elif lines.is_multilinestring():
        var mls = lines.as_multilinestring()
        for ln in mls.lines:
            lns.append(ln.copy())
    elif lines.is_geometrycollection():
        var gc = lines.as_geometrycollection()
        for g in gc.geoms:
            if g.is_linestring():
                lns.append(g.as_linestring())
            elif g.is_multilinestring():
                var mls2 = g.as_multilinestring()
                for ln in mls2.lines:
                    lns.append(ln.copy())

    var segs = List[PolygonizeSeg]()
    for ln in lns:
        if ln.coords.__len__() < 2: continue
        var i = 0
        while i < ln.coords.__len__() - 1:
            var a = ln.coords[i]
            var b = ln.coords[i + 1]
            segs.append(PolygonizeSeg(a[0], a[1], b[0], b[1]))
            i += 1

    # split segments at intersections
    var hit_seg = List[Int32]()
    var hit_t = List[Float64]()
    var i = 0
    while i < segs.__len__():
        var j = i + 1
        while j < segs.__len__():
            var pts = segment_intersections((segs[i].ax, segs[i].ay), (segs[i].bx, segs[i].by), (segs[j].ax, segs[j].ay), (segs[j].bx, segs[j].by))
            for P in pts:
                var ti = P[2]
                var tj = P[3]
                if ti > 0.0 and ti < 1.0:
                    hit_seg.append(Int32(i))
                    hit_t.append(ti)
                if tj > 0.0 and tj < 1.0:
                    hit_seg.append(Int32(j))
                    hit_t.append(tj)
            j += 1
        i += 1

    var verts = List[Tuple[Float64, Float64]]()
    var edges = List[PolygonizeDEdge]()
    var adj = List[List[Int32]]()

    # helper: sort ts and emit pieces
    var sidx = 0
    while sidx < segs.__len__():
        var ts = List[Float64]()
        ts.append(0.0)
        ts.append(1.0)
        var hh = 0
        while hh < hit_seg.__len__():
            if hit_seg[hh] == Int32(sidx):
                ts.append(hit_t[hh])
            hh += 1
        # insertion sort
        var a = 1
        while a < ts.__len__():
            var v = ts[a]
            var b = a - 1
            while b >= 0 and ts[b] > v:
                ts[b + 1] = ts[b]
                b -= 1
            ts[b + 1] = v
            a += 1
        var m = 0
        while m < ts.__len__() - 1:
            var t0 = ts[m]
            var t1 = ts[m + 1]
            if t1 - t0 > 1e-12:
                var ax = segs[sidx].ax + (segs[sidx].bx - segs[sidx].ax) * t0
                var ay = segs[sidx].ay + (segs[sidx].by - segs[sidx].ay) * t0
                var bx = segs[sidx].ax + (segs[sidx].bx - segs[sidx].ax) * t1
                var by = segs[sidx].ay + (segs[sidx].by - segs[sidx].ay) * t1
                var s_id = _polygonize_get_vid(ax, ay, verts)
                var d_id = _polygonize_get_vid(bx, by, verts)
                _polygonize_ensure_adj(adj, s_id)
                _polygonize_ensure_adj(adj, d_id)
                var dx = bx - ax
                var dy = by - ay
                var ei = Int32(edges.__len__())
                edges.append(PolygonizeDEdge(s_id, d_id, dx, dy, True, False))
                adj[s_id].append(ei)
                # add reverse too to allow traversal both ways
                var eri = Int32(edges.__len__())
                edges.append(PolygonizeDEdge(d_id, s_id, -dx, -dy, True, False))
                adj[d_id].append(eri)
            m += 1
        sidx += 1

    # prune dangles (degree < 2)
    var deg = List[Int32]()
    var v = 0
    while v < adj.__len__():
        deg.append(Int32(adj[v].__len__()))
        v += 1
    var changed = True
    while changed:
        changed = False
        var vi = 0
        while vi < adj.__len__():
            if deg[vi] > 0 and deg[vi] < 2:
                # remove all edges from this vertex
                var p = 0
                while p < adj[vi].__len__():
                    var eidx = adj[vi][p]
                    if edges[eidx].alive:
                        edges[eidx].alive = False
                        deg[vi] -= 1
                        var to = edges[eidx].dst
                        # Also remove the reverse directed edge to keep degree bookkeeping consistent.
                        var rp = 0
                        while rp < adj[to].__len__():
                            var ridx = adj[to][rp]
                            if edges[ridx].alive and edges[ridx].dst == Int32(vi):
                                edges[ridx].alive = False
                                if deg[to] > 0:
                                    deg[to] -= 1
                                break
                            rp += 1
                        changed = True
                    p += 1
                adj[vi] = List[Int32]()
            vi += 1

    var rings = List[List[Tuple[Float64, Float64]]]()
    var e = 0
    while e < edges.__len__():
        if not edges[e].alive or edges[e].used:
            e += 1
            continue
        var ring = List[Tuple[Float64, Float64]]()
        var start_e = Int32(e)
        var cur_e = start_e
        var start_v = edges[cur_e].src
        var bx = -edges[cur_e].dx
        var by = -edges[cur_e].dy
        var closed = False
        while True:
            edges[cur_e].used = True
            var vsrc = edges[cur_e].src
            ring.append((verts[vsrc][0], verts[vsrc][1]))
            var vdst = edges[cur_e].dst
            if vdst == start_v and ring.__len__() >= 4:
                ring.append((verts[start_v][0], verts[start_v][1]))
                closed = True
                break
            var ne = _polygonize_next_edge(adj, edges, vdst, bx, by)
            if ne == -1:
                break
            bx = -edges[ne].dx
            by = -edges[ne].dy
            cur_e = ne
        if closed and ring.__len__() >= 4:
            rings.append(ring)
        e += 1

    # assemble polygons from rings
    var shells = List[LinearRing]()
    var holes = List[LinearRing]()
    for r in rings:
        var area = signed_area_coords(r)
        if area > 0.0:
            shells.append(LinearRing(r))
        else:
            holes.append(LinearRing(r))

    var polys = List[Polygon]()
    var used_hole = List[Bool]()
    for _ in holes: used_hole.append(False)
    var si = 0
    while si < shells.__len__():
        var sh = shells[si]
        var sh_holes = List[LinearRing]()
        var hi = 0
        while hi < holes.__len__():
            if not used_hole[hi]:
                var hr = holes[hi]
                if hr.coords.__len__() > 0:
                    var pt = hr.coords[0]
                    var inside = point_in_ring(Point(pt[0], pt[1]), sh) != 0
                    if inside:
                        sh_holes.append(hr)
                        used_hole[hi] = True
            hi += 1
        polys.append(Polygon(sh, sh_holes))
        si += 1
    return GeometryCollection(polys)


fn polygonize_full(lines: Geometry) -> Tuple[GeometryCollection, MultiLineString, MultiLineString, MultiLineString]:
    # Collect all LineStrings from input
    var lns = List[LineString]()
    if lines.is_linestring():
        lns.append(lines.as_linestring())
    elif lines.is_multilinestring():
        var mls_in = lines.as_multilinestring()
        for ln in mls_in.lines:
            lns.append(ln.copy())
    elif lines.is_geometrycollection():
        var gc = lines.as_geometrycollection()
        for g in gc.geoms:
            if g.is_linestring():
                lns.append(g.as_linestring())
            elif g.is_multilinestring():
                var mls2 = g.as_multilinestring()
                for ln in mls2.lines:
                    lns.append(ln.copy())

    var segs = List[PolygonizeSeg]()
    for ln in lns:
        if ln.coords.__len__() < 2: continue
        var i = 0
        while i < ln.coords.__len__() - 1:
            var a = ln.coords[i]
            var b = ln.coords[i + 1]
            segs.append(PolygonizeSeg(a[0], a[1], b[0], b[1]))
            i += 1

    # split segments at intersections
    var hit_seg = List[Int32]()
    var hit_t = List[Float64]()
    var si = 0
    while si < segs.__len__():
        var sj = si + 1
        while sj < segs.__len__():
            var pts = segment_intersections((segs[si].ax, segs[si].ay), (segs[si].bx, segs[si].by), (segs[sj].ax, segs[sj].ay), (segs[sj].bx, segs[sj].by))
            for P in pts:
                var ti = P[2]
                var tj = P[3]
                if ti > 0.0 and ti < 1.0:
                    hit_seg.append(Int32(si))
                    hit_t.append(ti)
                if tj > 0.0 and tj < 1.0:
                    hit_seg.append(Int32(sj))
                    hit_t.append(tj)
            sj += 1
        si += 1

    # build atomic directed edges

    var verts = List[Tuple[Float64, Float64]]()
    var edges = List[PolygonizeDEdge]()
    var adj = List[List[Int32]]()

    # helper: sort ts and emit pieces
    var sidx = 0
    while sidx < segs.__len__():
        var ts = List[Float64]()
        ts.append(0.0)
        ts.append(1.0)
        var hh = 0
        while hh < hit_seg.__len__():
            if hit_seg[hh] == Int32(sidx):
                ts.append(hit_t[hh])
            hh += 1
        # insertion sort
        var a = 1
        while a < ts.__len__():
            var v = ts[a]
            var b = a - 1
            while b >= 0 and ts[b] > v:
                ts[b + 1] = ts[b]
                b -= 1
            ts[b + 1] = v
            a += 1
        var m = 0
        while m < ts.__len__() - 1:
            var t0 = ts[m]
            var t1 = ts[m + 1]
            if t1 - t0 > 1e-12:
                var ax = segs[sidx].ax + (segs[sidx].bx - segs[sidx].ax) * t0
                var ay = segs[sidx].ay + (segs[sidx].by - segs[sidx].ay) * t0
                var bx = segs[sidx].ax + (segs[sidx].bx - segs[sidx].ax) * t1
                var by = segs[sidx].ay + (segs[sidx].by - segs[sidx].ay) * t1
                var s_id = _polygonize_get_vid(ax, ay, verts)
                var d_id = _polygonize_get_vid(bx, by, verts)
                _polygonize_ensure_adj(adj, s_id)
                _polygonize_ensure_adj(adj, d_id)
                var dx = bx - ax
                var dy = by - ay
                var ei = Int32(edges.__len__())
                edges.append(PolygonizeDEdge(s_id, d_id, dx, dy, True, False))
                adj[s_id].append(ei)
                # add reverse too
                var eri = Int32(edges.__len__())
                edges.append(PolygonizeDEdge(d_id, s_id, -dx, -dy, True, False))
                adj[d_id].append(eri)
            m += 1
        sidx += 1

    # initial degrees for dangle detection
    var deg = List[Int32]()
    var v = 0
    while v < adj.__len__():
        deg.append(Int32(adj[v].__len__()))
        v += 1

    var dangles_seen = List[Tuple[Int32, Int32]]()
    var dangle_lines = List[LineString]()
    var cut_seen = List[Tuple[Int32, Int32]]()
    var cut_lines = List[LineString]()
    var ei = 0
    while ei < edges.__len__():
        ref e = edges[ei]
        if deg.__len__() > 0:
            if deg[e.src] == 1 or deg[e.dst] == 1:
                _polygonize_add_unique_edge(e.src, e.dst, dangles_seen, dangle_lines, verts)
        ei += 1

    # prune dangles (degree < 2)
    var changed = True
    while changed:
        changed = False
        var vi = 0
        while vi < adj.__len__():
            if deg[vi] > 0 and deg[vi] < 2:
                var p = 0
                while p < adj[vi].__len__():
                    var eidx = adj[vi][p]
                    if edges[eidx].alive:
                        edges[eidx].alive = False
                        deg[vi] -= 1
                        var to = edges[eidx].dst
                        # Also remove the reverse directed edge to keep degree bookkeeping consistent.
                        var rp = 0
                        while rp < adj[to].__len__():
                            var ridx = adj[to][rp]
                            if edges[ridx].alive and edges[ridx].dst == Int32(vi):
                                edges[ridx].alive = False
                                if deg[to] > 0:
                                    deg[to] -= 1
                                break
                            rp += 1
                        changed = True
                    p += 1
                adj[vi] = List[Int32]()
            vi += 1

    var rings = List[List[Tuple[Float64, Float64]]]()
    var e2 = 0
    while e2 < edges.__len__():
        if not edges[e2].alive or edges[e2].used:
            e2 += 1
            continue
        var ring = List[Tuple[Float64, Float64]]()
        var start_e = Int32(e2)
        var cur_e = start_e
        var start_v = edges[cur_e].src
        var bx = -edges[cur_e].dx
        var by = -edges[cur_e].dy
        var closed = False
        while True:
            edges[cur_e].used = True
            var vsrc = edges[cur_e].src
            ring.append((verts[vsrc][0], verts[vsrc][1]))
            var vdst = edges[cur_e].dst
            if vdst == start_v and ring.__len__() >= 4:
                ring.append((verts[start_v][0], verts[start_v][1]))
                closed = True
                break
            var ne = _polygonize_next_edge(adj, edges, vdst, bx, by)
            if ne == -1:
                break
            bx = -edges[ne].dx
            by = -edges[ne].dy
            cur_e = ne
        if closed and ring.__len__() >= 4:
            rings.append(ring.copy())
        e2 += 1

    # assemble polygons from rings (same as polygonize)
    var shells = List[LinearRing]()
    var holes = List[LinearRing]()
    for r in rings:
        var area = signed_area_coords(r)
        if area > 0.0:
            shells.append(LinearRing(r))
        else:
            holes.append(LinearRing(r))

    var polys = List[Polygon]()
    var used_hole = List[Bool]()
    for _ in holes: used_hole.append(False)
    var si2 = 0
    while si2 < shells.__len__():
        ref sh = shells[si2]
        var sh_holes = List[LinearRing]()
        var hi2 = 0
        while hi2 < holes.__len__():
            if not used_hole[hi2]:
                ref hr = holes[hi2]
                if hr.coords.__len__() > 0:
                    var pt = hr.coords[0]
                    var inside = point_in_ring(Point(pt[0], pt[1]), sh) == 1
                    if inside:
                        sh_holes.append(hr.copy())
                        used_hole[hi2] = True
            hi2 += 1
        polys.append(Polygon(sh.copy(), sh_holes))
        si2 += 1

    # classify remaining edges as cut-edges (alive but not used)
    var ce = 0
    while ce < edges.__len__():
        if edges[ce].alive and not edges[ce].used:
            _polygonize_add_unique_edge(edges[ce].src, edges[ce].dst, cut_seen, cut_lines, verts)
        ce += 1

    var poly_geoms = List[Geometry]()
    for p in polys:
        poly_geoms.append(Geometry(p.copy()))
    return (GeometryCollection(poly_geoms), MultiLineString(dangle_lines), MultiLineString(cut_lines), MultiLineString([]))


fn unary_union(geoms: List[Geometry]) -> Geometry:
    return _unary_union(geoms)


fn clamp01(t: Float64) -> Float64:
    if t < 0.0: return 0.0
    if t > 1.0: return 1.0
    return t


fn closest_on_seg(ax: Float64, ay: Float64, bx: Float64, by: Float64, px: Float64, py: Float64) -> Tuple[Float64, Float64, Float64]:
    var vx = bx - ax
    var vy = by - ay
    var vlen2 = vx * vx + vy * vy
    var t = 0.0
    if vlen2 > 0.0:
        t = ((px - ax) * vx + (py - ay) * vy) / vlen2
    t = clamp01(t)
    var cx = ax + t * vx
    var cy = ay + t * vy
    var dx = px - cx
    var dy = py - cy
    return (cx, cy, dx * dx + dy * dy)


fn nearest_points(a: Point, b: Point) -> Tuple[Point, Point]:
    return (a, b)


fn nearest_points(p: Point, ls: LineString) -> Tuple[Point, Point]:
    var best = 1.7976931348623157e308
    var bx = p.x
    var by = p.y
    for i in range(0, ls.coords.__len__() - 1):
        var a = ls.coords[i]
        var b = ls.coords[i + 1]
        var (cx, cy, d2) = closest_on_seg(a[0], a[1], b[0], b[1], p.x, p.y)
        if d2 < best:
            best = d2
            bx = cx
            by = cy
    return (Point(bx, by), p)


fn nearest_points(ls: LineString, p: Point) -> Tuple[Point, Point]:
    var (q, _p) = nearest_points(p, ls)
    return (q, p)


fn nearest_points(l1: LineString, l2: LineString) -> Tuple[Point, Point]:
    var best = 1.7976931348623157e308
    var p1 = Point(0.0, 0.0)
    var p2 = Point(0.0, 0.0)
    for i in range(0, l1.coords.__len__()):
        var px = l1.coords[i][0]
        var py = l1.coords[i][1]
        for j in range(0, l2.coords.__len__() - 1):
            var a = l2.coords[j]
            var b = l2.coords[j + 1]
            var (cx, cy, d2) = closest_on_seg(a[0], a[1], b[0], b[1], px, py)
            if d2 < best:
                best = d2
                p1 = Point(px, py)
                p2 = Point(cx, cy)
    for i in range(0, l2.coords.__len__()):
        var px = l2.coords[i][0]
        var py = l2.coords[i][1]
        for j in range(0, l1.coords.__len__() - 1):
            var a = l1.coords[j]
            var b = l1.coords[j + 1]
            var (cx, cy, d2) = closest_on_seg(a[0], a[1], b[0], b[1], px, py)
            if d2 < best:
                best = d2
                p1 = Point(px, py)
                p2 = Point(cx, cy)
    return (p1, p2)


fn clip_by_rect(poly: Polygon, xmin: Float64, ymin: Float64, xmax: Float64, ymax: Float64) -> Geometry:
    # Build rectangle polygon and intersect using our polygon clipper
    var rect = Polygon(LinearRing([
        (xmin, ymin), (xmax, ymin), (xmax, ymax), (xmin, ymax), (xmin, ymin)
    ]))
    return _poly_intersection(poly, rect)


fn clip_by_rect(ls: LineString, xmin: Float64, ymin: Float64, xmax: Float64, ymax: Float64) -> Geometry:
    # Liang–Barsky per segment; build one or more clipped polylines
    if ls.coords.__len__() < 2:
        return LineString(ls.coords)

    fn clip_seg(x0: Float64, y0: Float64, x1: Float64, y1: Float64,
                xmin: Float64, ymin: Float64, xmax: Float64, ymax: Float64,
               ) -> Tuple[Bool, Float64, Float64, Float64, Float64]:
        var t0 = 0.0
        var t1 = 1.0
        var dx = x1 - x0
        var dy = y1 - y0
        
        fn upd(p: Float64, q: Float64, mut t0: Float64, mut t1: Float64) -> Bool:
            if p == 0.0:
                if q < 0.0: return False
                return True
            var r = q / p
            if p < 0.0:
                if r > t1: return False
                if r > t0: t0 = r
            else:
                if r < t0: return False
                if r < t1: t1 = r
            return True
        
        if not upd(-dx, x0 - xmin, t0, t1): return (False, 0.0, 0.0, 0.0, 0.0)
        if not upd( dx, xmax - x0, t0, t1): return (False, 0.0, 0.0, 0.0, 0.0)
        if not upd(-dy, y0 - ymin, t0, t1): return (False, 0.0, 0.0, 0.0, 0.0)
        if not upd( dy, ymax - y0, t0, t1): return (False, 0.0, 0.0, 0.0, 0.0)

        var cx0 = x0 + t0 * dx
        var cy0 = y0 + t0 * dy
        var cx1 = x0 + t1 * dx
        var cy1 = y0 + t1 * dy
        return (True, cx0, cy0, cx1, cy1)

    var lines = List[List[Tuple[Float64, Float64]]]()
    var current = List[Tuple[Float64, Float64]]()
    for i in range(0, ls.coords.__len__() - 1):
        var a = ls.coords[i]
        var b = ls.coords[i + 1]
        var (ok, x0, y0, x1, y1) = clip_seg(a[0], a[1], b[0], b[1], xmin, ymin, xmax, ymax)
        if ok:
            # start a new part if current is empty or discontinuous from previous end
            if current.__len__() == 0:
                current.append((x0, y0))
                current.append((x1, y1))
            else:
                var last = current[current.__len__() - 1]
                if last[0] == x0 and last[1] == y0:
                    current.append((x1, y1))
                else:
                    lines.append(current)
                    current = List[Tuple[Float64, Float64]]()
                    current.append((x0, y0))
                    current.append((x1, y1))
        else:
            if current.__len__() > 0:
                lines.append(current)
                current = List[Tuple[Float64, Float64]]()
    if current.__len__() > 0:
        lines.append(current)

    if lines.__len__() == 0:
        return LineString([])
    if lines.__len__() == 1:
        return LineString(lines[0])
    var mls = List[LineString]()
    for part in lines:
        mls.append(LineString(part))
    return MultiLineString(mls)


fn orient(poly: Polygon, sign: Float64 = 1.0) -> Polygon:
    fn signed_area(coords: List[Tuple[Float64, Float64]]) -> Float64:
        if coords.__len__() < 2: return 0.0
        var s = 0.0
        for i in range(0, coords.__len__() - 1):
            ref a = coords[i]
            ref b = coords[i + 1]
            s += a[0] * b[1] - a[1] * b[0]
        return 0.5 * s

    fn ensure_orient(coords: List[Tuple[Float64, Float64]], want_positive: Bool) -> List[Tuple[Float64, Float64]]:
        var area = signed_area(coords)
        var is_positive = area > 0.0
        if want_positive == is_positive:
            return coords
        # reverse
        var out = List[Tuple[Float64, Float64]]()
        var i = coords.__len__() - 1
        while True:
            out.append(coords[i])
            if i == 0: break
            i -= 1
        return out

    var want_shell_pos = sign >= 0.0
    var new_shell = LinearRing(ensure_orient(poly.shell.coords, want_shell_pos))
    var new_holes = List[LinearRing]()
    for h in poly.holes:
        var hole_pos = not want_shell_pos
        new_holes.append(LinearRing(ensure_orient(h.coords, hole_pos)))
    return Polygon(new_shell, new_holes)


fn orient(geom: Geometry, _sign: Float64 = 1.0) -> Geometry:
    return geom


fn orient(mpoly: MultiPolygon, sign: Float64 = 1.0) -> MultiPolygon:
    var new_polys = List[Polygon]()
    for p in mpoly.polys:
        new_polys.append(orient(p, sign))
    return MultiPolygon(new_polys)


fn split(ls: LineString, pt: Point) -> GeometryCollection:
    # If point is not on the line, return original
    if ls.coords.size() < 2:
        return GeometryCollection([ls])
    var on_any = False
    var idx = 0
    for i in range(0, ls.coords.size() - 1):
        let a = ls.coords[i]
        let b = ls.coords[i + 1]
        if on_segment(a[0], a[1], b[0], b[1], pt.x, pt.y):
            on_any = True
            idx = i
            break
    if not on_any:
        return GeometryCollection([ls])
    # build two parts
    var left = List[Tuple[Float64, Float64]]()
    var right = List[Tuple[Float64, Float64]]()
    # left: from start to idx, then pt
    for i in range(0, idx + 1):
        left.append(ls.coords[i])
    if left.size() == 0 or left[left.size() - 1][0] != pt.x or left[left.size() - 1][1] != pt.y:
        left.append((pt.x, pt.y))
    # right: start from pt to end
    if ls.coords[idx + 1][0] != pt.x or ls.coords[idx + 1][1] != pt.y:
        right.append((pt.x, pt.y))
    for i in range(idx + 1, ls.coords.size()):
        right.append(ls.coords[i])
    var parts = List[Geometry]()
    if left.size() >= 2:
        parts.append(LineString(left))
    if right.size() >= 2:
        parts.append(LineString(right))
    if parts.size() == 0:
        return GeometryCollection([ls])
    return GeometryCollection(parts)


fn substring(geom: LineString, start_dist: Float64, end_dist: Float64, normalized: Bool = False) -> LineString:
    if geom.coords.size() == 0:
        return LineString([])
    if geom.coords.size() == 1:
        return LineString([geom.coords[0]])

    fn sqrt_f64(x: Float64) -> Float64:
        if x <= 0.0: return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r

    # total length
    var total = 0.0
    for i in range(0, geom.coords.size() - 1):
        let a = geom.coords[i]
        let b = geom.coords[i + 1]
        let dx = b[0] - a[0]
        let dy = b[1] - a[1]
        total += sqrt_f64(dx * dx + dy * dy)
    if total == 0.0:
        return LineString([geom.coords[0]])

    var s = start_dist
    var e = end_dist
    if normalized:
        s = s * total
        e = e * total
    # handle negatives as distances from end
    if s < 0.0: s = total + s
    if e < 0.0: e = total + e
    # clamp
    if s < 0.0: s = 0.0
    if e < 0.0: e = 0.0
    if s > total: s = total
    if e > total: e = total

    var reverse = False
    if s > e:
        let tmp = s
        s = e
        e = tmp
        reverse = True

    # collect points
    var acc = 0.0
    var out = List[Tuple[Float64, Float64]]()
    # find start point
    for i in range(0, geom.coords.size() - 1):
        let a = geom.coords[i]
        let b = geom.coords[i + 1]
        let dx = b[0] - a[0]
        let dy = b[1] - a[1]
        let seg = sqrt_f64(dx * dx + dy * dy)
        if acc + seg >= s:
            let t = if seg == 0.0 { 0.0 } else { (s - acc) / seg }
            let sx = a[0] + t * dx
            let sy = a[1] + t * dy
            out.append((sx, sy))
            # continue adding intermediate vertices until reaching end
            var here = acc + seg
            # add endpoints within (s, e)
            if here < e:
                out.append((b[0], b[1]))
            # add subsequent segments fully until overshoot
            var j = i + 1
            while j < geom.coords.size() - 1 and here < e:
                let na = geom.coords[j]
                let nb = geom.coords[j + 1]
                let ndx = nb[0] - na[0]
                let ndy = nb[1] - na[1]
                let nseg = sqrt_f64(ndx * ndx + ndy * ndy)
                if here + nseg <= e:
                    out.append((nb[0], nb[1]))
                    here += nseg
                    j += 1
                else:
                    let tt = if nseg == 0.0 { 0.0 } else { (e - here) / nseg }
                    let ex = na[0] + tt * ndx
                    let ey = na[1] + tt * ndy
                    out.append((ex, ey))
                    here = e
                    break
            break
        acc += seg

    if out.size() == 0:
        # start beyond total length after clamping -> return single endpoint at end
        out.append((geom.coords[geom.coords.size() - 1][0], geom.coords[geom.coords.size() - 1][1]))

    if reverse:
        # reverse coordinate order
        var rev = List[Tuple[Float64, Float64]]()
        var k = out.size() - 1
        while True:
            rev.append(out[k])
            if k == 0: break
            k -= 1
        out = rev

    return LineString(out)


fn shared_paths(a: LineString, b: LineString) -> GeometryCollection:
    return _shared_paths(a, b)


fn shortest_line(a: Point, b: Point) -> LineString:
    return LineString([(a.x, a.y), (b.x, b.y)])


fn shortest_line(p: Point, ls: LineString) -> LineString:
    let (q, _p) = nearest_points(p, ls)
    return LineString([(q.x, q.y), (p.x, p.y)])


fn shortest_line(ls: LineString, p: Point) -> LineString:
    let (q, _p) = nearest_points(ls, p)
    return LineString([(q.x, q.y), (p.x, p.y)])


fn shortest_line(l1: LineString, l2: LineString) -> LineString:
    let (a, b) = nearest_points(l1, l2)
    return LineString([(a.x, a.y), (b.x, b.y)])


fn triangulate(_geom: Geometry, _tolerance: Float64 = 0.0, _edges: Bool = False) -> GeometryCollection:
    return GeometryCollection([])


fn voronoi_diagram(_geom: Geometry, _envelope: Geometry = Geometry(), _tolerance: Float64 = 0.0, _edges: Bool = False) -> GeometryCollection:
    return GeometryCollection([])


fn validate(_geom: Geometry) -> Bool:
    return True
