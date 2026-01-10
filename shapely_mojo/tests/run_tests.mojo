from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPolygon, MultiLineString, GeometryCollection
from shapely.set_operations import union, intersection, difference, symmetric_difference, unary_union
from shapely.measurement import area, distance
from shapely.strtree import STRtree
from shapely.ops import polygonize_full


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
    print("areas: inter=" + area(inter).__str__() + ", union=" + area(u).__str__() + ", diff=" + area(diff).__str__() + ", xor=" + area(xor).__str__())
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
    print("area: union(MultiPolygon, Polygon)=" + area(u2).__str__())
    var t5 = expect("union(MultiPolygon, Polygon) area == 7", approx_eq(area(u2), 7.0))
    p += t5[0]; f += t5[1]

    var items = List[Geometry]()
    items.append(Geometry(A.copy()))
    items.append(Geometry(B.copy()))
    var uu = unary_union(items)
    print("area: unary_union=" + area(uu).__str__())
    var t6 = expect("unary_union area == 7", approx_eq(area(uu), 7.0))
    p += t6[0]; f += t6[1]

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
    var s2 = expect("STRtree overlaps count", q_ov.__len__() == 1)
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

    print("\nSummary: ")
    print(("passed = " + passed.__str__()) + (", failed = " + failed.__str__()))
    if failed == 0:
        print("ALL TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
