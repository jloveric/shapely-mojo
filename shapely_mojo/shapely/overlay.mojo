from shapely._geometry import Geometry
from shapely.geometry import Point, LinearRing, Polygon, MultiPolygon
from shapely.algorithms import (
    point_in_polygon,
    point_in_ring,
    segment_intersections,
    signed_area_coords,
)


struct DEdge(Copyable, Movable):
    var src: Int32
    var dst: Int32
    var dx: Float64
    var dy: Float64
    var include: Bool
    var used: Bool

    fn __init__(
        out self,
        src: Int32,
        dst: Int32,
        dx: Float64,
        dy: Float64,
        include: Bool,
        used: Bool,
    ):
        self.src = src
        self.dst = dst
        self.dx = dx
        self.dy = dy
        self.include = include
        self.used = used


fn abs_f64(x: Float64) -> Float64:
    if x < 0.0:
        return -x
    return x


fn sqrt_f64(x: Float64) -> Float64:
    if x <= 0.0:
        return 0.0
    var r = x
    var i = 0
    while i < 12:
        r = 0.5 * (r + x / r)
        i += 1
    return r


fn get_vid(
    x: Float64,
    y: Float64,
    mut verts: List[Tuple[Float64, Float64]],
    eps: Float64 = 1e-12,
) -> Int32:
    var i: Int = 0
    while i < verts.__len__():
        var v = verts[i]
        if abs_f64(v[0] - x) <= eps and abs_f64(v[1] - y) <= eps:
            return Int32(i)
        i += 1
    verts.append((x, y))
    return Int32(verts.__len__() - 1)


struct Segment(Copyable, Movable):
    var ax: Float64
    var ay: Float64
    var bx: Float64
    var by: Float64
    var owner: Int32
    var ts: List[Float64]

    fn __init__(
        out self,
        ax: Float64,
        ay: Float64,
        bx: Float64,
        by: Float64,
        owner: Int32,
    ):
        self.ax = ax
        self.ay = ay
        self.bx = bx
        self.by = by
        self.owner = owner
        self.ts = [0.0, 1.0]


fn add_ring_segments(r: LinearRing, owner: Int32, mut segs: List[Segment]):
    if r.coords.__len__() < 2:
        return
    var i = 0
    while i < r.coords.__len__() - 1:
        var a = r.coords[i]
        var b = r.coords[i + 1]
        segs.append(Segment(a[0], a[1], b[0], b[1], owner))
        i += 1


fn ensure_adj(mut adj: List[List[Int32]], vid: Int32):
    while adj.__len__() <= Int(vid):
        adj.append(List[Int32]())


fn emit_edge(
    ax: Float64,
    ay: Float64,
    bx: Float64,
    by: Float64,
    a: Polygon,
    b: Polygon,
    op: Int32,
    eps: Float64,
    mut verts: List[Tuple[Float64, Float64]],
    mut edges: List[DEdge],
    mut used: List[Bool],
    mut adj: List[List[Int32]],
):
    var sx = ax
    var sy = ay
    var dx = bx - ax
    var dy = by - ay
    var len = sqrt_f64(dx * dx + dy * dy)
    if len == 0.0:
        return
    var mx = (ax + bx) * 0.5
    var my = (ay + by) * 0.5
    var nx = -dy / len
    var ny = dx / len
    # Sample just to the left and right of the directed edge.
    # We include the directed edge only when the left side is inside the
    # operation result and the right side is outside, so that the traced
    # ring has the interior on the left.
    var lx = mx + nx * eps
    var ly = my + ny * eps
    var rx = mx - nx * eps
    var ry = my - ny * eps
    var insideA_L = point_in_polygon(make_point(lx, ly), a) != 0
    var insideB_L = point_in_polygon(make_point(lx, ly), b) != 0
    var insideA_R = point_in_polygon(make_point(rx, ry), a) != 0
    var insideB_R = point_in_polygon(make_point(rx, ry), b) != 0

    var include: Bool
    if op == 0:
        include = (insideA_L and insideB_L) and not (insideA_R and insideB_R)
    elif op == 1:
        include = (insideA_L or insideB_L) and not (insideA_R or insideB_R)
    elif op == 2:
        include = (insideA_L and not insideB_L) and not (insideA_R and not insideB_R)
    else:
        include = (
            ((insideA_L and not insideB_L) or (insideB_L and not insideA_L))
            and not ((insideA_R and not insideB_R) or (insideB_R and not insideA_R))
        )

    var s_id = get_vid(sx, sy, verts)
    var d_id = get_vid(bx, by, verts)
    ensure_adj(adj, s_id)
    ensure_adj(adj, d_id)
    var e_idx = Int32(edges.__len__())
    edges.append(
        DEdge(
            s_id,
            d_id,
            dx,
            dy,
            include,
            False,
        )
    )
    used.append(False)
    adj[s_id].append(e_idx)


fn make_point(x: Float64, y: Float64) -> Point:
    return Point(x, y)


fn compute_intersections(segs: List[Segment]) -> List[List[Float64]]:
    var adds = List[List[Float64]]()
    var idx = 0
    while idx < segs.__len__():
        adds.append(List[Float64]())
        idx += 1

    var i = 0
    while i < segs.__len__():
        var j = i + 1
        while j < segs.__len__():
            var pts = segment_intersections(
                (segs[i].ax, segs[i].ay),
                (segs[i].bx, segs[i].by),
                (segs[j].ax, segs[j].ay),
                (segs[j].bx, segs[j].by),
            )
            for P in pts:
                var t = P[2]
                var u = P[3]
                if t > 0.0 and t < 1.0:
                    adds[i].append(t)
                if u > 0.0 and u < 1.0:
                    adds[j].append(u)
            j += 1
        i += 1

    return adds.copy()


fn build_edges(
    a: Polygon, b: Polygon, op: Int32
) -> Tuple[
    List[Tuple[Float64, Float64]],
    List[DEdge],
    List[Bool],
    List[List[Int32]],
]:
    var segs = List[Segment]()
    add_ring_segments(a.shell, 0, segs)
    for h in a.holes:
        add_ring_segments(h, 0, segs)
    add_ring_segments(b.shell, 1, segs)
    for h2 in b.holes:
        add_ring_segments(h2, 1, segs)
    var adds = compute_intersections(segs)

    var verts = List[Tuple[Float64, Float64]]()
    var edges = List[DEdge]()
    var used = List[Bool]()
    var adj = List[List[Int32]]()

    var eps = 1e-9

    var i = 0
    while i < segs.__len__():
        var ts = List[Float64]()
        ts.append(0.0)
        ts.append(1.0)
        if i < adds.__len__():
            ref ap = adds[i]
            var ii = 0
            while ii < ap.__len__():
                ts.append(ap[ii])
                ii += 1
        # simple insertion sort
        var k = 1
        while k < ts.__len__():
            var tval = ts[k]
            var m = k - 1
            while m >= 0 and ts[m] > tval:
                ts[m + 1] = ts[m]
                m -= 1
            ts[m + 1] = tval
            k += 1
        # dedupe nearly-equal parameters
        var ts2 = List[Float64]()
        if ts.__len__() > 0:
            ts2.append(ts[0])
            var q = 1
            while q < ts.__len__():
                if abs_f64(ts[q] - ts2[ts2.__len__() - 1]) > 1e-12:
                    ts2.append(ts[q])
                q += 1
        ts = ts2.copy()
        var j = 0
        while j < ts.__len__() - 1:
            var t0 = ts[j]
            var t1 = ts[j + 1]
            if t1 - t0 > 1e-12:
                var ax = segs[i].ax + (segs[i].bx - segs[i].ax) * t0
                var ay = segs[i].ay + (segs[i].by - segs[i].ay) * t0
                var bx = segs[i].ax + (segs[i].bx - segs[i].ax) * t1
                var by = segs[i].ay + (segs[i].by - segs[i].ay) * t1
                emit_edge(
                    ax, ay, bx, by, a, b, op, eps, verts, edges, used, adj
                )
                emit_edge(
                    bx, by, ax, ay, a, b, op, eps, verts, edges, used, adj
                )
            j += 1
        i += 1

    return (verts.copy(), edges.copy(), used.copy(), adj.copy())


fn next_edge(
    adj: List[List[Int32]],
    edges: List[DEdge],
    used: List[Bool],
    at_vertex: Int32,
    bx: Float64,
    by: Float64,
) -> Int32:
    if Int(at_vertex) >= adj.__len__():
        return -1
    # Choose the outgoing edge that makes the smallest clockwise turn from
    # the incoming direction (bx, by). This traces faces with interior on the left.
    var best_idx: Int32 = -1
    var best_dot = -1.0e308
    var have_clockwise = False
    ref cand = adj[Int(at_vertex)]
    var i = 0
    while i < cand.__len__():
        var ei = cand[i]
        if not edges[Int(ei)].include or used[Int(ei)]:
            i += 1
            continue
        var w_x = edges[Int(ei)].dx
        var w_y = edges[Int(ei)].dy
        var cr = bx * w_y - by * w_x
        var dt = bx * w_x + by * w_y
        if cr < 0.0:
            # preferred: clockwise turn
            if not have_clockwise or dt > best_dot:
                best_idx = ei
                best_dot = dt
                have_clockwise = True
        elif not have_clockwise:
            # fallback: no clockwise options; choose the smallest counterclockwise
            # by maximizing dot as well.
            if best_idx == -1 or dt > best_dot:
                best_idx = ei
                best_dot = dt
        i += 1
    return best_idx


fn build_rings(
    verts: List[Tuple[Float64, Float64]],
    edges: List[DEdge],
    mut used: List[Bool],
    adj: List[List[Int32]],
) -> List[List[Tuple[Float64, Float64]]]:
    var rings = List[List[Tuple[Float64, Float64]]]()
    var e: Int = 0
    while e < edges.__len__():
        if not edges[e].include or used[e]:
            e += 1
            continue

        var ring = List[Tuple[Float64, Float64]]()
        var start_e = Int32(e)
        var cur_e = start_e
        var start_v = edges[Int(start_e)].src
        var bx = -edges[Int(cur_e)].dx
        var by = -edges[Int(cur_e)].dy
        var closed = False
        while True:
            var idx = Int(cur_e)
            used[idx] = True
            var vsrc = edges[idx].src
            ring.append((verts[vsrc][0], verts[vsrc][1]))
            var vdst = edges[idx].dst
            if vdst == start_v and ring.__len__() >= 3:
                ring.append((verts[start_v][0], verts[start_v][1]))
                closed = True
                break
            var ne = next_edge(adj, edges, used, vdst, bx, by)
            if ne == -1:
                break
            bx = -edges[Int(ne)].dx
            by = -edges[Int(ne)].dy
            cur_e = ne
        if closed and ring.__len__() >= 3:
            # Ensure explicit closure for downstream area computations.
            var first = ring[0]
            var last = ring[ring.__len__() - 1]
            if first[0] != last[0] or first[1] != last[1]:
                ring.append(first)
        if closed and ring.__len__() >= 4:
            rings.append(ring.copy())
        e += 1
    return rings.copy()


fn assemble_polygons(
    rings: List[List[Tuple[Float64, Float64]]]
) -> MultiPolygon:
    if rings.__len__() == 0:
        return MultiPolygon([])
    var shells = List[LinearRing]()
    var holes = List[LinearRing]()
    for r in rings:
        var area = signed_area_coords(r)
        if area > 0.0:
            shells.append(LinearRing(r))
        else:
            holes.append(LinearRing(r))
    if shells.__len__() == 0:
        return MultiPolygon([])
    # assign each hole to the smallest-area containing shell
    var polys = List[Polygon]()
    var used_hole = List[Bool]()
    for _ in holes:
        used_hole.append(False)
    # precompute shell areas (positive)
    var shell_areas = List[Float64]()
    for sh in shells:
        var ar = signed_area_coords(sh.coords)
        if ar < 0.0:
            ar = -ar
        shell_areas.append(ar)
    var i = 0
    while i < shells.__len__():
        ref sh = shells[i]
        var sh_holes = List[LinearRing]()
        var j = 0
        while j < holes.__len__():
            if not used_hole[j]:
                ref hr = holes[j]
                var pt = hr.coords[0]
                # find containing shell with minimal area
                var best_idx: Int32 = -1
                var best_area = 1.7976931348623157e308
                var si = 0
                while si < shells.__len__():
                    var inside = (
                        point_in_ring(Point(pt[0], pt[1]), shells[si]) != 0
                    )
                    if inside:
                        var sar = shell_areas[si]
                        if sar < best_area:
                            best_area = sar
                            best_idx = Int32(si)
                    si += 1
                if best_idx == Int32(i):
                    sh_holes.append(hr.copy())
                    used_hole[j] = True
            j += 1
        polys.append(Polygon(sh.copy(), sh_holes))
        i += 1
    if polys.__len__() == 1:
        return MultiPolygon([polys[0].copy()])
    return MultiPolygon(polys)


fn overlay_intersection(a: Polygon, b: Polygon) -> Geometry:
    var tmp0 = build_edges(a, b, 0)
    var verts = tmp0[0].copy()
    var edges = tmp0[1].copy()
    var used = tmp0[2].copy()
    var adj = tmp0[3].copy()
    var u2 = used.copy()
    var rings = build_rings(verts, edges, u2, adj)
    return Geometry(assemble_polygons(rings))


fn overlay_union(a: Polygon, b: Polygon) -> Geometry:
    var tmp1 = build_edges(a, b, 1)
    var verts = tmp1[0].copy()
    var edges = tmp1[1].copy()
    var used = tmp1[2].copy()
    var adj = tmp1[3].copy()
    var u2 = used.copy()
    var rings = build_rings(verts, edges, u2, adj)
    return Geometry(assemble_polygons(rings))


fn overlay_difference(a: Polygon, b: Polygon) -> Geometry:
    var tmp2 = build_edges(a, b, 2)
    var verts = tmp2[0].copy()
    var edges = tmp2[1].copy()
    var used = tmp2[2].copy()
    var adj = tmp2[3].copy()
    var u2 = used.copy()
    var rings = build_rings(verts, edges, u2, adj)
    return Geometry(assemble_polygons(rings))


fn overlay_xor(a: Polygon, b: Polygon) -> Geometry:
    var tmp3 = build_edges(a, b, 3)
    var verts = tmp3[0].copy()
    var edges = tmp3[1].copy()
    var used = tmp3[2].copy()
    var adj = tmp3[3].copy()
    var u2 = used.copy()
    var rings = build_rings(verts, edges, u2, adj)
    return Geometry(assemble_polygons(rings))
