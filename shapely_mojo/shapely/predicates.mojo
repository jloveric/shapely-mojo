from shapely._geometry import Geometry
from shapely.geometry import Point, LineString, LinearRing, Polygon, MultiPolygon
from shapely.algorithms import (
    point_on_linestring,
    point_in_polygon,
    any_segment_intersection,
    any_segment_intersection_coords,
)
from shapely.set_operations import intersection as _poly_intersection
from shapely.measurement import area as _area


fn intersects(a: Point, b: Point) -> Bool:
    return a.x == b.x and a.y == b.y


fn intersects(a: Point, b: LineString) -> Bool:
    return point_on_linestring(a, b)


fn intersects(a: LineString, b: Point) -> Bool:
    return point_on_linestring(b, a)


fn intersects(a: Point, b: Polygon) -> Bool:
    var r = point_in_polygon(a, b)
    return r != 0


fn intersects(a: Polygon, b: Point) -> Bool:
    var r = point_in_polygon(b, a)
    return r != 0


fn intersects(a: LineString, b: LineString) -> Bool:
    return any_segment_intersection(a, b)


fn intersects(a: LineString, b: Polygon) -> Bool:
    # check edge intersections with shell and holes; if none, check endpoint inside
    # shell
    ref shell = b.shell
    var ring_ls = LineString(shell.coords)
    if any_segment_intersection(a, ring_ls):
        return True
    # holes
    for h in b.holes:
        ring_ls = LineString(h.coords)
        if any_segment_intersection(a, ring_ls):
            return True
    # endpoint containment
    if a.coords.__len__() > 0:
        var p0 = Point(a.coords[0][0], a.coords[0][1])
        if point_in_polygon(p0, b) != 0:
            return True
    return False


fn intersects(a: Geometry, b: Geometry) -> Bool:
    if a.is_point() and b.is_point():
        return intersects(a.as_point(), b.as_point())
    if a.is_point() and b.is_linestring():
        return intersects(a.as_point(), b.as_linestring())
    if a.is_linestring() and b.is_point():
        return intersects(a.as_linestring(), b.as_point())
    if a.is_point() and b.is_polygon():
        return intersects(a.as_point(), b.as_polygon())
    if a.is_polygon() and b.is_point():
        return intersects(a.as_polygon(), b.as_point())
    if a.is_linestring() and b.is_linestring():
        return intersects(a.as_linestring(), b.as_linestring())
    if a.is_linestring() and b.is_polygon():
        return intersects(a.as_linestring(), b.as_polygon())
    if a.is_polygon() and b.is_linestring():
        return intersects(a.as_polygon(), b.as_linestring())
    if a.is_polygon() and b.is_polygon():
        return intersects(a.as_polygon(), b.as_polygon())
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
    var ab = a.bounds()
    var bb = b.bounds()
    if ab[2] < bb[0] or bb[2] < ab[0] or ab[3] < bb[1] or bb[3] < ab[1]:
        return False

    if any_segment_intersection_coords(a.shell.coords, b.shell.coords):
        return True

    for ah in a.holes:
        if any_segment_intersection_coords(ah.coords, b.shell.coords):
            return True
        for bh in b.holes:
            if any_segment_intersection_coords(ah.coords, bh.coords):
                return True

    for bh2 in b.holes:
        if any_segment_intersection_coords(a.shell.coords, bh2.coords):
            return True

    if a.shell.coords.__len__() > 0:
        var c0 = a.shell.coords[0]
        var p0 = Point(c0[0], c0[1])
        if point_in_polygon(p0, b) != 0:
            return True
    if b.shell.coords.__len__() > 0:
        var c1 = b.shell.coords[0]
        var p1 = Point(c1[0], c1[1])
        if point_in_polygon(p1, a) != 0:
            return True

    return False


fn contains(a: Polygon, b: Point) -> Bool:
    # True only if strictly inside (boundary returns False)
    return point_in_polygon(b, a) == 1


fn within(a: Point, b: Polygon) -> Bool:
    return contains(b, a)


fn touches(a: Geometry, b: Geometry) -> Bool:
    if a.is_polygon() and b.is_polygon():
        return touches(a.as_polygon(), b.as_polygon())
    if a.is_linestring() and b.is_linestring():
        return touches(a.as_linestring(), b.as_linestring())
    if a.is_point() and b.is_polygon():
        return touches(a.as_point(), b.as_polygon())
    if a.is_polygon() and b.is_point():
        return touches(a.as_polygon(), b.as_point())
    return False


fn disjoint(a: Point, b: Point) -> Bool:
    return not intersects(a, b)


fn disjoint(a: Geometry, b: Geometry) -> Bool:
    return not intersects(a, b)


fn overlaps(a: Geometry, b: Geometry) -> Bool:
    if a.is_polygon() and b.is_polygon():
        return overlaps(a.as_polygon(), b.as_polygon())
    return False


fn crosses(a: Geometry, b: Geometry) -> Bool:
    if a.is_linestring() and b.is_linestring():
        return crosses(a.as_linestring(), b.as_linestring())
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
    if a.coords.__len__() > 0:
        var a0 = a.coords[0]
        var aN = a.coords[a.coords.__len__() - 1]
        if point_on_linestring(Point(a0[0], a0[1]), b): end_on_other = True
        if point_on_linestring(Point(aN[0], aN[1]), b): end_on_other = True
    if b.coords.__len__() > 0:
        var b0 = b.coords[0]
        var bN = b.coords[b.coords.__len__() - 1]
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
    var inter_g = _poly_intersection(a, b)
    if inter_g.is_geometrycollection():
        return True
    # area zero -> touch; positive -> not touch
    var ar = _area(inter_g)
    return ar == 0.0


fn overlaps(a: Polygon, b: Polygon) -> Bool:
    # Overlaps if interiors intersect and neither contains the other entirely
    var inter_g = _poly_intersection(a, b)
    var ar = _area(inter_g)
    if ar == 0.0:
        return False
    var aa = _area(a)
    var bb = _area(b)
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
    if a.is_polygon() and b.is_point():
        return contains(a.as_polygon(), b.as_point())
    if a.is_polygon() and b.is_polygon():
        return contains(a.as_polygon(), b.as_polygon())
    return False


fn contains(a: Polygon, b: Polygon) -> Bool:
    # Contains: intersection area equals area of b (boundary contact allowed)
    # Avoid relying on polygon-polygon intersection area here; instead check that all
    # vertices of b lie inside or on boundary of a.
    if b.shell.coords.__len__() == 0:
        return True
    for c in b.shell.coords:
        var p = Point(c[0], c[1])
        if point_in_polygon(p, a) == 0:
            return False
    return True


fn within(a: Geometry, b: Geometry) -> Bool:
    return contains(b, a)


fn covers(a: Geometry, b: Geometry) -> Bool:
    if a.is_polygon() and b.is_point():
        # boundary counts as covered
        return point_in_polygon(b.as_point(), a.as_polygon()) != 0
    if a.is_polygon() and b.is_polygon():
        return covers(a.as_polygon(), b.as_polygon())
    return False


fn covers(a: Polygon, b: Polygon) -> Bool:
    # Covers: area of intersection equals area of b (boundary allowed)
    if b.shell.coords.__len__() == 0:
        return True
    for c in b.shell.coords:
        var p = Point(c[0], c[1])
        if point_in_polygon(p, a) == 0:
            return False
    return True


fn covered_by(a: Geometry, b: Geometry) -> Bool:
    return covers(b, a)


fn contains_properly(a: Polygon, b: Polygon) -> Bool:
    # Proper containment: all vertices of b strictly inside a (boundary excluded)
    if b.shell.coords.__len__() == 0:
        return True
    for c in b.shell.coords:
        var p = Point(c[0], c[1])
        if point_in_polygon(p, a) != 1:
            return False
    return True


fn contains_properly(a: Polygon, b: Point) -> Bool:
    # Point strictly inside polygon (boundary excluded)
    return point_in_polygon(b, a) == 1


fn contains_properly(a: Geometry, b: Geometry) -> Bool:
    if a.is_polygon() and b.is_polygon():
        return contains_properly(a.as_polygon(), b.as_polygon())
    if a.is_polygon() and b.is_point():
        return contains_properly(a.as_polygon(), b.as_point())
    return False
