from python import Python, PythonObject

from shapely._geometry import Geometry
from shapely.creation import box
from shapely.constructive import buffer
from shapely.set_operations import union, intersection, difference, symmetric_difference
from shapely.geometry import Point, LinearRing, Polygon
from shapely.strtree import STRtree


fn _mk_skew_poly(x0: Float64, y0: Float64, w: Float64, h: Float64) -> Polygon:
    var sx = 0.35 * w
    var pts = List[Tuple[Float64, Float64]]()
    pts.append((x0, y0))
    pts.append((x0 + w, y0 + 0.2 * h))
    pts.append((x0 + w + sx, y0 + h))
    pts.append((x0 + sx, y0 + 0.8 * h))
    pts.append((x0, y0))
    return Polygon(LinearRing(pts))


fn _now_ns() raises -> Int64:
    var time: PythonObject = Python.import_module("time")
    var t: PythonObject = time.perf_counter_ns()
    var ns = Python.py_long_as_ssize_t(t)
    return Int64(ns)


fn _print_result(name: String, iters: Int, elapsed_ns: Int64):
    var secs = Float64(elapsed_ns) / 1.0e9
    print("RESULT\t" + name + "\t" + secs.__str__() + "\t" + iters.__str__())


fn main() raises:
    # buffer
    var p = box(0.0, 0.0, 2.0, 2.0)
    var warm = 200
    var iters = 2000
    var i = 0
    while i < warm:
        var _ = buffer(p.copy(), 0.25, 8)
        i += 1
    var t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = buffer(p.copy(), 0.25, 8)
        i += 1
    var t1 = _now_ns()
    _print_result("buffer_box", iters, t1 - t0)

    var pg = Geometry(p.copy())
    i = 0
    while i < warm:
        var _ = buffer(pg.copy(), 0.25, 8)
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = buffer(pg.copy(), 0.25, 8)
        i += 1
    t1 = _now_ns()
    _print_result("buffer_box_geom", iters, t1 - t0)

    var coords = List[Tuple[Float64, Float64]]()
    coords.append((0.0, 0.0))
    coords.append((2.0, 0.0))
    coords.append((3.0, 1.0))
    coords.append((2.0, 2.0))
    coords.append((0.0, 2.0))
    coords.append((-1.0, 1.0))
    coords.append((0.0, 0.0))
    var conv = Polygon(LinearRing(coords))
    i = 0
    while i < warm:
        var _ = buffer(conv.copy(), 0.25, 8)
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = buffer(conv.copy(), 0.25, 8)
        i += 1
    t1 = _now_ns()
    _print_result("buffer_convex_poly", iters, t1 - t0)

    var coords2 = List[Tuple[Float64, Float64]]()
    coords2.append((0.0, 0.0))
    coords2.append((3.0, 0.0))
    coords2.append((3.0, 3.0))
    coords2.append((2.0, 3.0))
    coords2.append((2.0, 1.0))
    coords2.append((1.0, 1.0))
    coords2.append((1.0, 3.0))
    coords2.append((0.0, 3.0))
    coords2.append((0.0, 0.0))
    var conc = Polygon(LinearRing(coords2))
    i = 0
    while i < warm:
        var _ = buffer(conc.copy(), 0.25, 8)
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = buffer(conc.copy(), 0.25, 8)
        i += 1
    t1 = _now_ns()
    _print_result("buffer_concave_poly", iters, t1 - t0)

    var skew = _mk_skew_poly(0.0, 0.0, 3.0, 3.0)
    i = 0
    while i < warm:
        var _ = buffer(skew.copy(), 0.25, 8)
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = buffer(skew.copy(), 0.25, 8)
        i += 1
    t1 = _now_ns()
    _print_result("buffer_skew_poly", iters, t1 - t0)

    var shell = List[Tuple[Float64, Float64]]()
    shell.append((0.0, 0.0))
    shell.append((10.0, 0.0))
    shell.append((10.0, 10.0))
    shell.append((0.0, 10.0))
    shell.append((0.0, 0.0))
    var hole1 = List[Tuple[Float64, Float64]]()
    hole1.append((2.0, 2.0))
    hole1.append((4.0, 2.0))
    hole1.append((4.0, 4.0))
    hole1.append((2.0, 4.0))
    hole1.append((2.0, 2.0))
    var hole2 = List[Tuple[Float64, Float64]]()
    hole2.append((6.0, 2.0))
    hole2.append((8.0, 2.0))
    hole2.append((8.0, 4.0))
    hole2.append((6.0, 4.0))
    hole2.append((6.0, 2.0))
    var holes = List[LinearRing]()
    holes.append(LinearRing(hole1))
    holes.append(LinearRing(hole2))
    var poly_holes = Polygon(LinearRing(shell), holes)
    i = 0
    while i < warm:
        var _ = buffer(poly_holes.copy(), 0.25, 8)
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = buffer(poly_holes.copy(), 0.25, 8)
        i += 1
    t1 = _now_ns()
    _print_result("buffer_poly_holes", iters, t1 - t0)

    # boolean ops
    var a = Geometry(box(0.0, 0.0, 2.0, 2.0))
    var b = Geometry(box(1.0, 0.8, 3.0, 2.6))

    warm = 500
    iters = 5000

    i = 0
    while i < warm:
        var _ = union(a.copy(), b.copy())
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = union(a.copy(), b.copy())
        i += 1
    t1 = _now_ns()
    _print_result("union", iters, t1 - t0)

    i = 0
    while i < warm:
        var _ = intersection(a.copy(), b.copy())
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = intersection(a.copy(), b.copy())
        i += 1
    t1 = _now_ns()
    _print_result("intersection", iters, t1 - t0)

    i = 0
    while i < warm:
        var _ = difference(a.copy(), b.copy())
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = difference(a.copy(), b.copy())
        i += 1
    t1 = _now_ns()
    _print_result("difference", iters, t1 - t0)

    i = 0
    while i < warm:
        var _ = symmetric_difference(a.copy(), b.copy())
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = symmetric_difference(a.copy(), b.copy())
        i += 1
    t1 = _now_ns()
    _print_result("symmetric_difference", iters, t1 - t0)

    # STRtree: query + nearest
    var items = List[Geometry]()
    var x = 0
    while x < 40:
        var y = 0
        while y < 25:
            var xmin = Float64(x) * 1.5
            var ymin = Float64(y) * 1.5
            items.append(Geometry(box(xmin, ymin, xmin + 1.0, ymin + 1.0)))
            y += 1
        x += 1

    var tree = STRtree(items)
    var qpoly = Geometry(box(10.0, 10.0, 30.0, 25.0))
    var qpt = Point(31.2, 9.7)

    warm = 200
    iters = 2000

    i = 0
    while i < warm:
        var _ = tree.query_items(qpoly.copy(), "intersects")
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = tree.query_items(qpoly.copy(), "intersects")
        i += 1
    t1 = _now_ns()
    _print_result("strtree_query_intersects_idx", iters, t1 - t0)

    i = 0
    while i < warm:
        var _ = tree.query(qpoly.copy(), "intersects")
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = tree.query(qpoly.copy(), "intersects")
        i += 1
    t1 = _now_ns()
    _print_result("strtree_query_intersects_geom", iters, t1 - t0)

    i = 0
    while i < warm:
        var _ = tree.nearest_item(qpt.copy())
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = tree.nearest_item(qpt.copy())
        i += 1
    t1 = _now_ns()
    _print_result("strtree_nearest_point_idx", iters, t1 - t0)

    i = 0
    while i < warm:
        var _ = tree.nearest(qpt.copy())
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = tree.nearest(qpt.copy())
        i += 1
    t1 = _now_ns()
    _print_result("strtree_nearest_point_geom", iters, t1 - t0)

    var skew_items = List[Geometry]()
    x = 0
    while x < 40:
        var y = 0
        while y < 25:
            var xmin2 = Float64(x) * 1.5
            var ymin2 = Float64(y) * 1.5
            skew_items.append(Geometry(_mk_skew_poly(xmin2, ymin2, 1.0, 1.0)))
            y += 1
        x += 1

    var skew_tree = STRtree(skew_items)
    var skew_qpoly = Geometry(_mk_skew_poly(10.0, 10.0, 20.0, 15.0))

    i = 0
    while i < warm:
        var _ = skew_tree.query_items(skew_qpoly.copy(), "intersects")
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = skew_tree.query_items(skew_qpoly.copy(), "intersects")
        i += 1
    t1 = _now_ns()
    _print_result("strtree_query_intersects_skew_idx", iters, t1 - t0)

    i = 0
    while i < warm:
        var _ = skew_tree.query(skew_qpoly.copy(), "intersects")
        i += 1
    t0 = _now_ns()
    i = 0
    while i < iters:
        var _ = skew_tree.query(skew_qpoly.copy(), "intersects")
        i += 1
    t1 = _now_ns()
    _print_result("strtree_query_intersects_skew_geom", iters, t1 - t0)
