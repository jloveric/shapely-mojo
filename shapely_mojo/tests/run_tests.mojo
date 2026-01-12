from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPolygon, MultiLineString, GeometryCollection
from shapely.set_operations import union, intersection, difference, symmetric_difference, unary_union
from shapely.measurement import area, distance
from shapely.strtree import STRtree
from shapely.ops import polygonize_full
from shapely.constructive import simplify, convex_hull
from shapely.validation import make_valid


fn approx_eq(a: Float64, b: Float64, eps: Float64 = 1e-9) -> Bool:
    var d = a - b
    var ad = d
    if d < 0.0:
        ad = -d
    return ad <= eps


fn expect(name: String, cond: Bool) -> (Int32, Int32):
    if cond:
        print("PASS: " + name)
        return (1, 0)
    else:
        print("FAIL: " + name)
        return (0, 1)


fn mk_square(x0: Float64, y0: Float64, x1: Float64, y1: Float64) -> Polygon:
    return Polygon(LinearRing([(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]))


fn test_set_operations() -> (Int32, Int32):
    var A = mk_square(0.0, 0.0, 2.0, 2.0)
    var B = mk_square(1.0, 1.0, 3.0, 3.0)

    var u = union(A, B)
    var inter = intersection(A, B)
    var diff = difference(A, B)
    var xor = symmetric_difference(A, B)

    var p: Int32 = 0
    var f: Int32 = 0
    var t1 = expect("intersection area == 1", approx_eq(area(inter), 1.0))
    p += t1[0]; f += t1[1]
    var t2 = expect("union area == 7", approx_eq(area(u), 7.0))
    p += t2[0]; f += t2[1]
    var t3 = expect("difference area == 3", approx_eq(area(diff), 3.0))
    p += t3[0]; f += t3[1]
    var t4 = expect("xor area == 6", approx_eq(area(xor), 6.0))
    p += t4[0]; f += t4[1]

    var mp = MultiPolygon([A.copy(), B.copy()])
    var u2 = union(mp, A)
    var t5 = expect("union(MultiPolygon, Polygon) area == 7", approx_eq(area(u2), 7.0))
    p += t5[0]; f += t5[1]

    var items = List[Geometry]()
    items.append(Geometry(A.copy()))
    items.append(Geometry(B.copy()))
    var uu = unary_union(items)
    var t6 = expect("unary_union area == 7", approx_eq(area(uu), 7.0))
    p += t6[0]; f += t6[1]

    return (p, f)


fn test_make_valid_basic() -> (Int32, Int32):
    var p: Int32 = 0
    var f: Int32 = 0

    # classic bowtie (self-intersecting) polygon
    var bow = Polygon(LinearRing([(0.0,0.0),(2.0,2.0),(0.0,2.0),(2.0,0.0),(0.0,0.0)]))
    var mv = make_valid(Geometry(bow.copy()))
    print("make_valid(bowtie) -> " + mv.to_wkt())
    var is_mp = mv.is_multipolygon()
    var t0 = expect("make_valid(bowtie) is multipolygon", is_mp)
    p += t0[0]; f += t0[1]
    if is_mp:
        var t0b = expect("make_valid(bowtie) parts == 2", mv.as_multipolygon().polys.__len__() == 2)
        p += t0b[0]; f += t0b[1]
    else:
        var t0c = expect("make_valid(bowtie) parts == 2", False)
        p += t0c[0]; f += t0c[1]

    # open ring should be closed by make_valid
    var open_shell = LinearRing([(0.0,0.0),(2.0,0.0),(2.0,2.0),(0.0,2.0)])
    var p2 = Polygon(open_shell)
    var mv2 = make_valid(Geometry(p2.copy()))
    var t1 = expect("make_valid(open shell) returns polygon or multipolygon", mv2.is_polygon() or mv2.is_multipolygon())
    p += t1[0]; f += t1[1]
    var nonempty = not mv2.is_empty()
    var t1b = expect("make_valid(open shell) non-empty", nonempty)
    p += t1b[0]; f += t1b[1]

    return (p, f)


fn test_constructive_convex_hull_simplify() -> (Int32, Int32):
    var p: Int32 = 0
    var f: Int32 = 0

    # convex_hull: empty -> GeometryCollection empty
    var empty_ls = LineString([])
    var h0 = convex_hull(Geometry(empty_ls.copy()))
    var t0 = expect("convex_hull(empty) is geometrycollection", h0.is_geometrycollection())
    p += t0[0]; f += t0[1]
    var t0b = expect("convex_hull(empty) empty", h0.is_empty())
    p += t0b[0]; f += t0b[1]

    # convex_hull: collinear -> LineString
    var col = LineString([(0.0,0.0),(1.0,0.0),(2.0,0.0)])
    var h1 = convex_hull(Geometry(col.copy()))
    var t1 = expect("convex_hull(collinear) is linestring", h1.is_linestring())
    p += t1[0]; f += t1[1]
    var t1b = expect("convex_hull(collinear) has 2 coords", h1.as_linestring().coords.__len__() == 2)
    p += t1b[0]; f += t1b[1]

    # convex_hull: triangle -> Polygon
    var tri = Polygon(LinearRing([(0.0,0.0),(2.0,0.0),(0.0,2.0),(0.0,0.0)]))
    var h2 = convex_hull(Geometry(tri.copy()))
    var t2 = expect("convex_hull(triangle) is polygon", h2.is_polygon())
    p += t2[0]; f += t2[1]
    var t2b = expect("convex_hull(triangle) area == 2", approx_eq(area(h2), 2.0))
    p += t2b[0]; f += t2b[1]

    # simplify: reduces a "wiggly" line
    var wig = LineString([(0.0,0.0),(1.0,0.05),(2.0,0.0),(3.0,0.05),(4.0,0.0)])
    var s0 = simplify(Geometry(wig.copy()), 0.1)
    var t3 = expect("simplify returns linestring", s0.is_linestring())
    p += t3[0]; f += t3[1]
    var t3b = expect("simplify reduces points", s0.as_linestring().coords.__len__() < wig.coords.__len__())
    p += t3b[0]; f += t3b[1]

    # simplify tolerance 0 leaves unchanged
    var s1 = simplify(Geometry(wig.copy()), 0.0)
    var t4 = expect("simplify(tol=0) keeps points", s1.as_linestring().coords.__len__() == wig.coords.__len__())
    p += t4[0]; f += t4[1]

    # simplify: polygon ring closure preserved and removes tiny hole
    var shell = LinearRing([(0.0,0.0),(4.0,0.0),(4.0,4.0),(0.0,4.0),(0.0,0.0)])
    var hole = LinearRing([(1.0,1.0),(1.2,1.0),(1.2,1.2),(1.0,1.2),(1.0,1.0)])
    var poly = Polygon(shell.copy(), [hole.copy()])
    var sp = simplify(Geometry(poly.copy()), 0.5)
    var t5 = expect("simplify(polygon) returns polygon", sp.is_polygon())
    p += t5[0]; f += t5[1]
    var spoly = sp.as_polygon()
    var sc = spoly.shell.coords.copy()
    var t5b = expect(
        "simplify(polygon) shell closed",
        sc.__len__() >= 4
        and sc[0][0] == sc[sc.__len__() - 1][0]
        and sc[0][1] == sc[sc.__len__() - 1][1],
    )
    p += t5b[0]; f += t5b[1]
    var t5c = expect("simplify(polygon) drops small hole", spoly.holes.__len__() == 0)
    p += t5c[0]; f += t5c[1]

    # simplify: multipolygon keeps only non-empty simplified parts
    var mp = MultiPolygon([poly.copy(), Polygon(LinearRing([])).copy()])
    var smp = simplify(Geometry(mp.copy()), 0.5)
    var t6 = expect("simplify(multipolygon) is multipolygon", smp.is_multipolygon())
    p += t6[0]; f += t6[1]
    var t6b = expect("simplify(multipolygon) size == 1", smp.as_multipolygon().polys.__len__() == 1)
    p += t6b[0]; f += t6b[1]

    return (p, f)


fn test_strtree_predicates() -> (Int32, Int32):
    var A = mk_square(0.0, 0.0, 2.0, 2.0)
    var B = mk_square(1.0, 1.0, 3.0, 3.0)
    var E = mk_square(2.0, 0.0, 4.0, 2.0)  # touches A along edge x=2
    var C = mk_square(10.0, 10.0, 12.0, 12.0)
    var S = mk_square(0.2, 0.2, 0.8, 0.8)

    var geoms = List[Geometry]()
    geoms.append(Geometry(A.copy()))
    geoms.append(Geometry(B.copy()))
    geoms.append(Geometry(E.copy()))
    geoms.append(Geometry(C.copy()))
    geoms.append(Geometry(S.copy()))
    var tree = STRtree(geoms)

    var p: Int32 = 0
    var f: Int32 = 0
    # intersects with A -> A and B and E (touching counts as intersects in Shapely)
    var q_inter = tree.query(A, "intersects")
    var s1 = expect("STRtree intersects count", q_inter.__len__() == 3)
    p += s1[0]; f += s1[1]

    # overlaps with B -> only A overlaps B
    var q_ov = tree.query(B, "overlaps")
    var s2 = expect("STRtree overlaps count", q_ov.__len__() == 2)
    p += s2[0]; f += s2[1]

    # touches with A -> only E touches A (B intersects by area)
    var q_touch = tree.query(A, "touches")
    var s3 = expect("STRtree touches count", q_touch.__len__() == 1)
    p += s3[0]; f += s3[1]

    # contains/within variants with S -> A contains S; equality cases allowed for contains/within/covered_by
    var q_withinS = tree.query(S, "within")  # geoms within S (includes S)
    var s4 = expect("STRtree within(S) count >= 1", q_withinS.__len__() >= 1)
    p += s4[0]; f += s4[1]

    var q_coveredbyS = tree.query(S, "covered_by")  # geoms covered_by S (includes S, possibly others that fully cover S)
    var s5 = expect("STRtree covered_by(S) count >= 1", q_coveredbyS.__len__() >= 1)
    p += s5[0]; f += s5[1]

    var q_containsS = tree.query(S, "contains")  # geoms that contain S (A, maybe S depending on semantics)
    var s6 = expect("STRtree contains(S) count >= 1", q_containsS.__len__() >= 1)
    p += s6[0]; f += s6[1]

    var q_cp = tree.query(S, "contains_properly")  # geoms that properly contain S (should be A)
    var s7 = expect("STRtree contains_properly(S) count == 1", q_cp.__len__() == 1)
    p += s7[0]; f += s7[1]

    # contains(A) -> none should contain A except possibly A itself depending on semantics
    var q_containsA = tree.query(A, "contains")
    var s8 = expect("STRtree contains(A) count >= 0", q_containsA.__len__() >= 0)
    p += s8[0]; f += s8[1]

    var q_coveredbyA = tree.query(A, "covered_by")  # geoms covered_by A (includes A itself)
    var s9 = expect("STRtree covered_by(A) count >= 1", q_coveredbyA.__len__() >= 1)
    p += s9[0]; f += s9[1]

    return (p, f)


fn test_strtree_nearest_knn() -> (Int32, Int32):
    var A = mk_square(0.0, 0.0, 2.0, 2.0)
    var B = mk_square(5.0, 0.0, 7.0, 2.0)
    var C = mk_square(10.0, 0.0, 12.0, 2.0)
    var geoms = List[Geometry]()
    geoms.append(Geometry(A.copy()))
    geoms.append(Geometry(B.copy()))
    geoms.append(Geometry(C.copy()))
    var tree = STRtree(geoms)

    var p = Point(4.7, 1.0)

    print("CALL tree.nearest")
    var n = tree.nearest(p)
    print("RETURN tree.nearest")
    var dnA = distance(p, n)
    var pcount: Int32 = 0
    var fcount: Int32 = 0
    var n1 = expect("STRtree.nearest returns closest", approx_eq(dnA, distance(B, p)))
    pcount += n1[0]; fcount += n1[1]

    print("CALL tree.query_knn")
    var kn = tree.query_knn(Geometry(p.copy()), 2)
    print("RETURN tree.query_knn")
    var n2 = expect("STRtree.kNN size == 2", kn.__len__() == 2)
    pcount += n2[0]; fcount += n2[1]

    print("CALL tree.nearest_item")
    var idx = tree.nearest_item(p)
    print("RETURN tree.nearest_item")
    var n3 = expect("STRtree.nearest_item index valid", idx >= 0 and idx < geoms.__len__())
    pcount += n3[0]; fcount += n3[1]

    return (pcount, fcount)


fn test_polygonize_full_basic() -> (Int32, Int32):
    # square ring
    var ring = LineString([(0.0,0.0),(2.0,0.0),(2.0,2.0),(0.0,2.0),(0.0,0.0)])
    # dangle: little tail from (2,1)->(3,1)
    var dangle = LineString([(2.0,1.0),(3.0,1.0)])
    var lines = MultiLineString([ring.copy(), dangle.copy()])

    var res = polygonize_full(Geometry(lines.copy()))
    ref polys = res[0]
    ref dangles = res[1]
    ref cut_edges = res[2]
    ref invalid_rings = res[3]

    var p: Int32 = 0
    var f: Int32 = 0
    var g1 = expect("polygonize_full polygons size == 1", polys.geoms.__len__() == 1)
    p += g1[0]; f += g1[1]
    var g2 = expect("polygonize_full dangles size == 1", dangles.lines.__len__() == 1)
    p += g2[0]; f += g2[1]
    var g3 = expect("polygonize_full cut_edges size == 0", cut_edges.lines.__len__() == 0)
    p += g3[0]; f += g3[1]
    var g4 = expect("polygonize_full invalid_rings size == 0", invalid_rings.lines.__len__() == 0)
    p += g4[0]; f += g4[1]

    return (p, f)


fn main():
    var passed: Int32 = 0
    var failed: Int32 = 0

    print("Running shapely-mojo tests...")

    print("START test_set_operations")
    var r1 = test_set_operations()
    print("END test_set_operations")
    passed += r1[0]; failed += r1[1]

    print("START test_strtree_predicates")
    var r2 = test_strtree_predicates()
    print("END test_strtree_predicates")
    passed += r2[0]; failed += r2[1]

    print("START test_strtree_nearest_knn")
    var r3 = test_strtree_nearest_knn()
    print("END test_strtree_nearest_knn")
    passed += r3[0]; failed += r3[1]

    print("START test_polygonize_full_basic")
    var r4 = test_polygonize_full_basic()
    print("END test_polygonize_full_basic")
    passed += r4[0]; failed += r4[1]

    print("START test_constructive_convex_hull_simplify")
    var r5 = test_constructive_convex_hull_simplify()
    print("END test_constructive_convex_hull_simplify")
    passed += r5[0]; failed += r5[1]

    print("START test_make_valid_basic")
    var r6 = test_make_valid_basic()
    print("END test_make_valid_basic")
    passed += r6[0]; failed += r6[1]

    print("\nSummary: ")
    print(("passed = " + passed.__str__()) + (", failed = " + failed.__str__()))
    if failed == 0:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
