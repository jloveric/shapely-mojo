from shapely._geometry import Geometry
from shapely.geometry import Point, LinearRing, Polygon, MultiPolygon
from shapely.algorithms import point_in_polygon, point_in_ring, segment_intersections, signed_area_coords


struct DEdge:
    var src: Int32
    var dst: Int32
    var dx: Float64
    var dy: Float64
    var include: Bool
    var used: Bool


fn abs_f64(x: Float64) -> Float64:
    if x < 0.0: return -x
    return x


fn sqrt_f64(x: Float64) -> Float64:
    if x <= 0.0: return 0.0
    var r = x
    var i = 0
    while i < 12:
        r = 0.5 * (r + x / r)
        i += 1
    return r


fn get_vid(x: Float64, y: Float64, inout verts: List[Tuple[Float64, Float64]], eps: Float64 = 1e-12) -> Int32:
    var i = 0
    while i < verts.size():
        let v = verts[i]
        if abs_f64(v[0] - x) <= eps and abs_f64(v[1] - y) <= eps:
            return i as Int32
        i += 1
    verts.push_back((x, y))
    return (verts.size() - 1) as Int32


struct Segment:
    var ax: Float64
    var ay: Float64
    var bx: Float64
    var by: Float64
    var owner: Int32
    var ts: List[Float64]

    fn __init__(inout self, ax: Float64, ay: Float64, bx: Float64, by: Float64, owner: Int32):
        self.ax = ax
        self.ay = ay
        self.bx = bx
        self.by = by
        self.owner = owner
        self.ts = [0.0, 1.0]


fn add_ring_segments(r: LinearRing, owner: Int32, inout segs: List[Segment]):
    if r.coords.size() < 2: return
    var i = 0
    while i < r.coords.size() - 1:
        let a = r.coords[i]
        let b = r.coords[i + 1]
        segs.push_back(Segment(a[0], a[1], b[0], b[1], owner))
        i += 1


fn make_point(x: Float64, y: Float64) -> Point:
    return Point(x, y)


fn add_intersections(inout segs: List[Segment]):
    var i = 0
    while i < segs.size():
        var j = i + 1
        while j < segs.size():
            let si = segs[i]
            let sj = segs[j]
            let pts = segment_intersections((si.ax, si.ay), (si.bx, si.by), (sj.ax, sj.ay), (sj.bx, sj.by))
            for P in pts:
                let t = P[2]
                let u = P[3]
                if t > 0.0 and t < 1.0:
                    segs[i].ts.push_back(t)
                if u > 0.0 and u < 1.0:
                    segs[j].ts.push_back(u)
            j += 1
        i += 1


fn build_edges(a: Polygon, b: Polygon, op: Int32) -> (List[Tuple[Float64, Float64]], List[DEdge], List[List[Int32]]):
    var segs = List[Segment]()
    add_ring_segments(a.shell, 0, segs)
    for h in a.holes: add_ring_segments(h, 0, segs)
    add_ring_segments(b.shell, 1, segs)
    for h2 in b.holes: add_ring_segments(h2, 1, segs)
    add_intersections(segs)

    var verts = List[Tuple[Float64, Float64]]()
    var edges = List[DEdge]()
    var adj = List[List[Int32]]()

    fn ensure_adj(inout adj: List[List[Int32]], vid: Int32):
        while adj.size() <= (vid as Int):
            adj.push_back(List[Int32]())

    let eps = 1e-9
    fn emit_edge(ax: Float64, ay: Float64, bx: Float64, by: Float64):
        let sx = ax
        let sy = ay
        let dx = bx - ax
        let dy = by - ay
        let len = sqrt_f64(dx * dx + dy * dy)
        if len == 0.0: return
        let mx = (ax + bx) * 0.5
        let my = (ay + by) * 0.5
        let nx = -dy / len
        let ny = dx / len
        let sxp = mx + nx * eps
        let syp = my + ny * eps
        let insideA = point_in_polygon(make_point(sxp, syp), a) != 0
        let insideB = point_in_polygon(make_point(sxp, syp), b) != 0
        var include = False
        if op == 0:
            include = insideA and insideB
        elif op == 1:
            include = insideA or insideB
        elif op == 2:
            include = insideA and not insideB
        else:
            include = (insideA and not insideB) or (insideB and not insideA)
        let s_id = get_vid(sx, sy, verts)
        let d_id = get_vid(bx, by, verts)
        ensure_adj(adj, s_id)
        ensure_adj(adj, d_id)
        let e_idx = edges.size() as Int32
        edges.push_back(DEdge(s_id, d_id, dx, dy, include, False))
        adj[s_id].push_back(e_idx)

    var i = 0
    while i < segs.size():
        var ts = segs[i].ts
        # simple insertion sort
        var k = 1
        while k < ts.size():
            var tval = ts[k]
            var m = k - 1
            while m >= 0 and ts[m] > tval:
                ts[m + 1] = ts[m]
                m -= 1
            ts[m + 1] = tval
            k += 1
        # dedupe nearly-equal parameters
        var ts2 = List[Float64]()
        if ts.size() > 0:
            ts2.push_back(ts[0])
            var q = 1
            while q < ts.size():
                if abs_f64(ts[q] - ts2[ts2.size() - 1]) > 1e-12:
                    ts2.push_back(ts[q])
                q += 1
        ts = ts2
        var j = 0
        while j < ts.size() - 1:
            let t0 = ts[j]
            let t1 = ts[j + 1]
            if t1 - t0 > 1e-12:
                let ax = segs[i].ax + (segs[i].bx - segs[i].ax) * t0
                let ay = segs[i].ay + (segs[i].by - segs[i].ay) * t0
                let bx = segs[i].ax + (segs[i].bx - segs[i].ax) * t1
                let by = segs[i].ay + (segs[i].by - segs[i].ay) * t1
                emit_edge(ax, ay, bx, by)
                emit_edge(bx, by, ax, ay)
            j += 1
        i += 1

    return (verts, edges, adj)


fn next_edge(adj: List[List[Int32]], edges: List[DEdge], at_vertex: Int32, bx: Float64, by: Float64) -> Int32:
    if (at_vertex as Int) >= adj.size():
        return -1
    var best_idx: Int32 = -1
    var best_cross = -1.0e308
    var best_dot = -1.0e308
    let cand = adj[at_vertex]
    var i = 0
    while i < cand.size():
        let ei = cand[i]
        if not edges[ei].include or edges[ei].used:
            i += 1
            continue
        let w_x = edges[ei].dx
        let w_y = edges[ei].dy
        let cr = bx * w_y - by * w_x
        let dt = bx * w_x + by * w_y
        if cr > 0.0:
            if best_idx == -1 or best_cross < 0.0 or dt > best_dot:
                best_idx = ei
                best_cross = cr
                best_dot = dt
        elif best_idx == -1 and best_cross < 0.0:
            if dt > best_dot:
                best_idx = ei
                best_cross = cr
                best_dot = dt
        i += 1
    return best_idx


fn build_rings(verts: List[Tuple[Float64, Float64]], inout edges: List[DEdge], adj: List[List[Int32]]) -> List[List[Tuple[Float64, Float64]]]:
    var rings = List[List[Tuple[Float64, Float64]]]()
    var e = 0
    while e < edges.size():
        if not edges[e].include or edges[e].used:
            e += 1
            continue
        var ring = List[Tuple[Float64, Float64]]()
        var start_e = e as Int32
        var cur_e = start_e
        let sx = verts[edges[cur_e].src][0]
        let sy = verts[edges[cur_e].src][1]
        var bx = -edges[cur_e].dx
        var by = -edges[cur_e].dy
        while True:
            edges[cur_e].used = True
            let vsrc = edges[cur_e].src
            ring.push_back((verts[vsrc][0], verts[vsrc][1]))
            let vdst = edges[cur_e].dst
            let ne = next_edge(adj, edges, vdst, bx, by)
            if ne == -1:
                break
            bx = -edges[ne].dx
            by = -edges[ne].dy
            if ne == start_e:
                ring.push_back((verts[edges[ne].src][0], verts[edges[ne].src][1]))
                break
            cur_e = ne
        if ring.size() >= 4:
            rings.push_back(ring)
        e += 1
    return rings


fn assemble_polygons(rings: List[List[Tuple[Float64, Float64]]]) -> Geometry:
    if rings.size() == 0:
        return Polygon(LinearRing(List[Tuple[Float64, Float64]]()))
    var shells = List[LinearRing]()
    var holes = List[LinearRing]()
    for r in rings:
        let area = signed_area_coords(r)
        if area > 0.0:
            shells.push_back(LinearRing(r))
        else:
            holes.push_back(LinearRing(r))
    if shells.size() == 0:
        return Polygon(LinearRing(List[Tuple[Float64, Float64]]()))
    # assign each hole to the smallest-area containing shell
    var polys = List[Polygon]()
    var used_hole = List[Bool]()
    for _ in holes: used_hole.push_back(False)
    # precompute shell areas (positive)
    var shell_areas = List[Float64]()
    for sh in shells:
        var ar = signed_area_coords(sh.coords)
        if ar < 0.0: ar = -ar
        shell_areas.push_back(ar)
    var i = 0
    while i < shells.size():
        let sh = shells[i]
        var sh_holes = List[LinearRing]()
        var j = 0
        while j < holes.size():
            if not used_hole[j]:
                let hr = holes[j]
                let pt = hr.coords[0]
                # find containing shell with minimal area
                var best_idx: Int32 = -1
                var best_area = 1.7976931348623157e308
                var si = 0
                while si < shells.size():
                    let inside = point_in_ring(Point(pt[0], pt[1]), shells[si]) != 0
                    if inside:
                        let sar = shell_areas[si]
                        if sar < best_area:
                            best_area = sar
                            best_idx = si as Int32
                    si += 1
                if best_idx == (i as Int32):
                    sh_holes.push_back(hr)
                    used_hole[j] = True
            j += 1
        polys.push_back(Polygon(sh, sh_holes))
        i += 1
    if polys.size() == 1:
        return polys[0]
    return MultiPolygon(polys)


fn overlay_intersection(a: Polygon, b: Polygon) -> Geometry:
    let (verts, edges, adj) = build_edges(a, b, 0)
    var e2 = edges
    let rings = build_rings(verts, e2, adj)
    return assemble_polygons(rings)


fn overlay_union(a: Polygon, b: Polygon) -> Geometry:
    let (verts, edges, adj) = build_edges(a, b, 1)
    var e2 = edges
    let rings = build_rings(verts, e2, adj)
    return assemble_polygons(rings)


fn overlay_difference(a: Polygon, b: Polygon) -> Geometry:
    let (verts, edges, adj) = build_edges(a, b, 2)
    var e2 = edges
    let rings = build_rings(verts, e2, adj)
    return assemble_polygons(rings)


fn overlay_xor(a: Polygon, b: Polygon) -> Geometry:
    let (verts, edges, adj) = build_edges(a, b, 3)
    var e2 = edges
    let rings = build_rings(verts, e2, adj)
    return assemble_polygons(rings)
