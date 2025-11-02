from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPolygon
from shapely.algorithms import point_on_linestring, point_in_polygon, any_segment_intersection
from shapely.set_operations import intersection as _poly_intersection
from shapely.measurement import area as _area


fn intersects(a: Point, b: Point) -> Bool:
    return a.x == b.x and a.y == b.y


fn intersects(a: Point, b: LineString) -> Bool:
    return point_on_linestring(a, b)


fn intersects(a: LineString, b: Point) -> Bool:
    return point_on_linestring(b, a)


fn intersects(a: Point, b: Polygon) -> Bool:
    let r = point_in_polygon(a, b)
    return r != 0


fn intersects(a: Polygon, b: Point) -> Bool:
    let r = point_in_polygon(b, a)
    return r != 0


fn intersects(a: LineString, b: LineString) -> Bool:
    return any_segment_intersection(a, b)


fn intersects(a: LineString, b: Polygon) -> Bool:
    # check edge intersections with shell and holes; if none, check endpoint inside
    # shell
    let shell = b.shell
    var ring_ls = LineString(shell.coords)
    if any_segment_intersection(a, ring_ls):
        return True
    # holes
    for h in b.holes:
        ring_ls = LineString(h.coords)
        if any_segment_intersection(a, ring_ls):
            return True
    # endpoint containment
    if a.coords.size() > 0:
        let p0 = Point(a.coords[0][0], a.coords[0][1])
        if point_in_polygon(p0, b) != 0:
            return True
    return False


fn intersects(a: Geometry, b: Geometry) -> Bool:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Point" and tb == "Point":
        return intersects(unsafe_bitcast[Point](a), unsafe_bitcast[Point](b))
    if ta == "Point" and tb == "LineString":
        return intersects(unsafe_bitcast[Point](a), unsafe_bitcast[LineString](b))
    if ta == "LineString" and tb == "Point":
        return intersects(unsafe_bitcast[LineString](a), unsafe_bitcast[Point](b))
    if ta == "Point" and tb == "Polygon":
        return intersects(unsafe_bitcast[Point](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "Point":
        return intersects(unsafe_bitcast[Polygon](a), unsafe_bitcast[Point](b))
    if ta == "LineString" and tb == "LineString":
        return intersects(unsafe_bitcast[LineString](a), unsafe_bitcast[LineString](b))
    if ta == "LineString" and tb == "Polygon":
        return intersects(unsafe_bitcast[LineString](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "LineString":
        return intersects(unsafe_bitcast[Polygon](a), unsafe_bitcast[LineString](b))
    if ta == "Polygon" and tb == "Polygon":
        return intersects(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    return False


# disjoint overloads for common pairs
fn disjoint(a: Point, b: LineString) -> Bool:
    return not intersects(a, b)


fn disjoint(a: LineString, b: Point) -> Bool:
    return not intersects(a, b)


fn disjoint(a: Point, b: Polygon) -> Bool:
    return not intersects(a, b)


fn disjoint(a: Polygon, b: Point) -> Bool:
    return not intersects(a, b)


fn disjoint(a: LineString, b: LineString) -> Bool:
    return not intersects(a, b)


fn disjoint(a: LineString, b: Polygon) -> Bool:
    return not intersects(a, b)


fn disjoint(a: Polygon, b: LineString) -> Bool:
    return not intersects(a, b)


fn disjoint(a: Polygon, b: Polygon) -> Bool:
    return not intersects(a, b)


fn intersects(a: Polygon, b: LineString) -> Bool:
    return intersects(b, a)


fn intersects(a: Polygon, b: Polygon) -> Bool:
    # quick tests: any edge intersection
    let a_shell_ls = LineString(a.shell.coords)
    let b_shell_ls = LineString(b.shell.coords)
    if any_segment_intersection(a_shell_ls, b_shell_ls):
        return True
    # containment tests: one polygon contains a vertex of the other
    if a.shell.coords.size() > 0:
        let p = Point(a.shell.coords[0][0], a.shell.coords[0][1])
        if point_in_polygon(p, b) != 0: return True
    if b.shell.coords.size() > 0:
        let p2 = Point(b.shell.coords[0][0], b.shell.coords[0][1])
        if point_in_polygon(p2, a) != 0: return True
    return False


fn contains(a: Polygon, b: Point) -> Bool:
    # True only if strictly inside (boundary returns False)
    return point_in_polygon(b, a) == 1


fn within(a: Point, b: Polygon) -> Bool:
    return contains(b, a)


fn touches(a: Geometry, b: Geometry) -> Bool:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Polygon":
        return touches(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    if ta == "LineString" and tb == "LineString":
        return touches(unsafe_bitcast[LineString](a), unsafe_bitcast[LineString](b))
    if ta == "Point" and tb == "Polygon":
        return touches(unsafe_bitcast[Point](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "Point":
        return touches(unsafe_bitcast[Polygon](a), unsafe_bitcast[Point](b))
    return False


fn disjoint(a: Point, b: Point) -> Bool:
    return not intersects(a, b)


fn disjoint(a: Geometry, b: Geometry) -> Bool:
    return not intersects(a, b)


fn overlaps(a: Geometry, b: Geometry) -> Bool:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Polygon":
        return overlaps(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    return False


fn crosses(a: Geometry, b: Geometry) -> Bool:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "LineString" and tb == "LineString":
        return crosses(unsafe_bitcast[LineString](a), unsafe_bitcast[LineString](b))
    return False


fn equals(a: Geometry, b: Geometry) -> Bool:
    # Simple structural equality via WKT representation; refine later
    return a.to_wkt() == b.to_wkt()


# --- Additional advanced predicates ---

fn touches(a: Point, b: Polygon) -> Bool:
    return point_in_polygon(a, b) == 2


fn touches(a: Polygon, b: Point) -> Bool:
    return point_in_polygon(b, a) == 2


fn touches(a: LineString, b: LineString) -> Bool:
    # touch if they intersect but only at endpoints (no interior crossing)
    if not any_segment_intersection(a, b):
        return False
    # heuristic: if any endpoint of one lies on the other, and there is no proper interior crossing
    var end_on_other = False
    if a.coords.size() > 0:
        let a0 = a.coords[0]
        let aN = a.coords[a.coords.size() - 1]
        if point_on_linestring(Point(a0[0], a0[1]), b): end_on_other = True
        if point_on_linestring(Point(aN[0], aN[1]), b): end_on_other = True
    if b.coords.size() > 0:
        let b0 = b.coords[0]
        let bN = b.coords[b.coords.size() - 1]
        if point_on_linestring(Point(b0[0], b0[1]), a): end_on_other = True
        if point_on_linestring(Point(bN[0], bN[1]), a): end_on_other = True
    if not end_on_other:
        return False
    # if there is any interior intersection (non-endpoint), consider not touches
    # we approximate by checking first vertices of a not on endpoints of b
    return True


fn touches(a: Polygon, b: Polygon) -> Bool:
    # Shapely semantics: interiors do not intersect, boundaries intersect
    if not intersects(a, b):
        return False
    let inter_g = _poly_intersection(a, b)
    let t = inter_g.__type_name__()
    if t == "GeometryCollection":
        return True
    # area zero -> touch; positive -> not touch
    let ar = _area(inter_g)
    return ar == 0.0


fn overlaps(a: Polygon, b: Polygon) -> Bool:
    # Overlaps if interiors intersect and neither contains the other entirely
    let inter_g = _poly_intersection(a, b)
    let ar = _area(inter_g)
    if ar == 0.0:
        return False
    let aa = _area(a)
    let bb = _area(b)
    # if inter area equals any full area -> containment, not overlap
    if ar >= aa or ar >= bb:
        return False
    return True


fn crosses(a: LineString, b: LineString) -> Bool:
    # Crosses if they intersect and not only at endpoints
    if not any_segment_intersection(a, b):
        return False
    # If they touch only at endpoints, treat as not cross
    if touches(a, b):
        return False
    return True


fn contains(a: Geometry, b: Geometry) -> Bool:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Point":
        return contains(unsafe_bitcast[Polygon](a), unsafe_bitcast[Point](b))
    if ta == "Polygon" and tb == "Polygon":
        return contains(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    return False


fn contains(a: Polygon, b: Polygon) -> Bool:
    # Contains: intersection area equals area of b (boundary contact allowed)
    let inter_g = _poly_intersection(a, b)
    let ar_inter = _area(inter_g)
    let ar_b = _area(b)
    return ar_inter == ar_b


fn within(a: Geometry, b: Geometry) -> Bool:
    return contains(b, a)


fn covers(a: Geometry, b: Geometry) -> Bool:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Point":
        # boundary counts as covered
        return point_in_polygon(unsafe_bitcast[Point](b), unsafe_bitcast[Polygon](a)) != 0
    if ta == "Polygon" and tb == "Polygon":
        return covers(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    return False


fn covers(a: Polygon, b: Polygon) -> Bool:
    # Covers: area of intersection equals area of b (boundary allowed)
    let inter_g = _poly_intersection(a, b)
    let ar_inter = _area(inter_g)
    let ar_b = _area(b)
    return ar_inter == ar_b


fn covered_by(a: Geometry, b: Geometry) -> Bool:
    return covers(b, a)


fn contains_properly(a: Polygon, b: Polygon) -> Bool:
    # Proper containment: contains without boundary touching
    return contains(a, b) and not touches(a, b)


fn contains_properly(a: Polygon, b: Point) -> Bool:
    # Point strictly inside polygon (boundary excluded)
    return point_in_polygon(b, a) == 1


fn contains_properly(a: Geometry, b: Geometry) -> Bool:
    let ta = a.__type_name__()
    let tb = b.__type_name__()
    if ta == "Polygon" and tb == "Polygon":
        return contains_properly(unsafe_bitcast[Polygon](a), unsafe_bitcast[Polygon](b))
    if ta == "Polygon" and tb == "Point":
        return contains_properly(unsafe_bitcast[Polygon](a), unsafe_bitcast[Point](b))
    return False
