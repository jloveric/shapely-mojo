from shapely._geometry import Geometry
from shapely.geometry import (
    LinearRing,
    Polygon,
    Point,
    LineString,
    MultiPolygon,
    MultiLineString,
    GeometryCollection,
)
from shapely.measurement import distance as _distance
from shapely.predicates import (
    intersects as _intersects,
    overlaps as _overlaps,
    touches as _touches,
    contains as _contains,
    covers as _covers,
    contains_properly as _contains_properly,
)


struct STRtree:
    var geoms: List[Geometry]
    var boxes: List[
        Tuple[Float64, Float64, Float64, Float64]
    ]  # (minx, miny, maxx, maxy)

    fn __init__(out self, geoms: List[Geometry]):
        self.geoms = geoms.copy()
        self.boxes = List[Tuple[Float64, Float64, Float64, Float64]]()

        fn _bounds_of(g: Geometry) -> (Float64, Float64, Float64, Float64):
            return g.bounds()

        for g in self.geoms:
            self.boxes.append(_bounds_of(g))

    fn _env_intersects(
        self,
        a: Tuple[Float64, Float64, Float64, Float64],
        b: Tuple[Float64, Float64, Float64, Float64],
    ) -> Bool:
        return not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1])

    fn _env_dist2(
        self,
        a: Tuple[Float64, Float64, Float64, Float64],
        b: Tuple[Float64, Float64, Float64, Float64],
    ) -> Float64:
        var dx = 0.0
        if a[2] < b[0]:
            dx = b[0] - a[2]
        elif b[2] < a[0]:
            dx = a[0] - b[2]
        var dy = 0.0
        if a[3] < b[1]:
            dy = b[1] - a[3]
        elif b[3] < a[1]:
            dy = a[1] - b[3]
        return dx * dx + dy * dy

    fn _sqrt_f64(self, x: Float64) -> Float64:
        if x <= 0.0:
            return 0.0
        var r = x
        var i = 0
        while i < 12:
            r = 0.5 * (r + x / r)
            i += 1
        return r

    fn _ceil_div(self, a: Int, b: Int) -> Int:
        var q = a / b
        if a % b != 0:
            q += 1
        return q

    fn _node_union_boxes(
        self, ids: List[Int32]
    ) -> (Float64, Float64, Float64, Float64):
        # Stub retained for compatibility; not used in naive implementation
        return (0.0, 0.0, 0.0, 0.0)

    fn _geom_union_boxes(
        self, ids: List[Int32]
    ) -> (Float64, Float64, Float64, Float64):
        var minx = 1.7976931348623157e308
        var miny = 1.7976931348623157e308
        var maxx = -1.7976931348623157e308
        var maxy = -1.7976931348623157e308
        var i = 0
        while i < ids.__len__():
            var b = self.boxes[ids[i]]
            if b[0] < minx:
                minx = b[0]
            if b[1] < miny:
                miny = b[1]
            if b[2] > maxx:
                maxx = b[2]
            if b[3] > maxy:
                maxy = b[3]
            i += 1
        return (minx, miny, maxx, maxy)

    fn _qsort_geom_by_minx(self, idx: List[Int32], lo: Int, hi: Int):
        var i = lo
        var j = hi
        var pivot = self.boxes[idx[(lo + hi) / 2]][0]
        while i <= j:
            while self.boxes[idx[i]][0] < pivot:
                i += 1
            while self.boxes[idx[j]][0] > pivot:
                j -= 1
            if i <= j:
                var tmp = idx[i]
                idx[i] = idx[j]
                idx[j] = tmp
                i += 1
                j -= 1
        if lo < j:
            self._qsort_geom_by_minx(idx, lo, j)
        if i < hi:
            self._qsort_geom_by_minx(idx, i, hi)

    fn _qsort_geom_slice_by_miny(
        self, idx: List[Int32], start: Int, end_excl: Int
    ):
        fn sort_range(self_ref: STRtree, a: List[Int32], lo: Int, hi: Int):
            var i = lo
            var j = hi
            var pivot = self_ref.boxes[a[(lo + hi) / 2]][1]
            while i <= j:
                while self_ref.boxes[a[i]][1] < pivot:
                    i += 1
                while self_ref.boxes[a[j]][1] > pivot:
                    j -= 1
                if i <= j:
                    var tmp = a[i]
                    a[i] = a[j]
                    a[j] = tmp
                    i += 1
                    j -= 1
            if lo < j:
                sort_range(self_ref, a, lo, j)
            if i < hi:
                sort_range(self_ref, a, i, hi)

        if end_excl - start > 1:
            sort_range(self, idx, start, end_excl - 1)

    fn _build(self, max_children: Int32):
        # Naive implementation does not build a tree
        return

    fn query(self, _target: Geometry) -> List[Geometry]:
        var tb = _target.bounds()
        var out = List[Geometry]()
        var i = 0
        for b in self.boxes:
            if self._env_intersects(b, tb):
                out.append(self.geoms[i].copy())
            i += 1
        return out.copy()

    fn query(self, _target: Polygon) -> List[Geometry]:
        var tb = _target.bounds()
        var out = List[Geometry]()
        var i = 0
        for b in self.boxes:
            if self._env_intersects(b, tb):
                out.append(self.geoms[i].copy())
            i += 1
        return out.copy()

    fn query(self, _target: Geometry, predicate: String) -> List[Geometry]:
        var cands = self.query(_target)
        if predicate == "intersects":
            var out = List[Geometry]()
            for g in cands:
                if _intersects(g, _target):
                    out.append(g.copy())
            return out.copy()
        elif predicate == "touches":
            var out2 = List[Geometry]()
            for g in cands:
                if _touches(g, _target):
                    out2.append(g.copy())
            return out2.copy()
        elif predicate == "overlaps":
            var out3 = List[Geometry]()
            for g in cands:
                if _overlaps(g, _target):
                    out3.append(g.copy())
            return out3.copy()
        elif predicate == "contains":
            var out4 = List[Geometry]()
            for g in cands:
                if _contains(g, _target):
                    out4.append(g.copy())
            return out4.copy()
        elif predicate == "within":
            var out5 = List[Geometry]()
            for g in cands:
                if _contains(_target, g):
                    out5.append(g.copy())
            return out5.copy()
        elif predicate == "covers":
            var out6 = List[Geometry]()
            for g in cands:
                if _covers(g, _target):
                    out6.append(g.copy())
            return out6.copy()
        elif predicate == "covered_by":
            var out7 = List[Geometry]()
            for g in cands:
                if _covers(_target, g):
                    out7.append(g.copy())
            return out7.copy()
        elif predicate == "contains_properly":
            var out10 = List[Geometry]()
            for g in cands:
                if _contains_properly(g, _target):
                    out10.append(g.copy())
            return out10.copy()
        return cands.copy()

    fn query(self, _target: Polygon, predicate: String) -> List[Geometry]:
        var cands = self.query(_target)
        var tgt = Geometry(_target.copy())
        if predicate == "intersects":
            var out = List[Geometry]()
            for g in cands:
                if _intersects(g, tgt):
                    out.append(g.copy())
            return out.copy()
        elif predicate == "touches":
            var out2 = List[Geometry]()
            for g in cands:
                if _touches(g, tgt):
                    out2.append(g.copy())
            return out2.copy()
        elif predicate == "overlaps":
            var out3 = List[Geometry]()
            for g in cands:
                if _overlaps(g, tgt):
                    out3.append(g.copy())
            return out3.copy()
        elif predicate == "contains":
            var out4 = List[Geometry]()
            for g in cands:
                if _contains(g, tgt):
                    out4.append(g.copy())
            return out4.copy()
        elif predicate == "within":
            var out5 = List[Geometry]()
            for g in cands:
                if _contains(tgt, g):
                    out5.append(g.copy())
            return out5.copy()
        elif predicate == "covers":
            var out6 = List[Geometry]()
            for g in cands:
                if _covers(g, tgt):
                    out6.append(g.copy())
            return out6.copy()
        elif predicate == "covered_by":
            var out7 = List[Geometry]()
            for g in cands:
                if _covers(tgt, g):
                    out7.append(g.copy())
            return out7.copy()
        elif predicate == "contains_properly":
            var out10 = List[Geometry]()
            for g in cands:
                if _contains_properly(g, tgt):
                    out10.append(g.copy())
            return out10.copy()
        return cands.copy()

    fn nearest(self, _target: Geometry) -> Geometry:
        if self.boxes.__len__() == 0:
            return Geometry(Polygon(LinearRing(List[Tuple[Float64, Float64]]())))
        var best = 1.7976931348623157e308
        var best_idx = 0
        var i = 0
        for b in self.boxes:
            var gd = _distance(self.geoms[i], _target)
            var d2 = gd * gd
            if d2 < best:
                best = d2
                best_idx = i
                if best == 0.0:
                    return self.geoms[best_idx].copy()
            i += 1
        return self.geoms[best_idx].copy()

    fn nearest(self, _target: Point) -> Geometry:
        if self.boxes.__len__() == 0:
            return Geometry(Polygon(LinearRing(List[Tuple[Float64, Float64]]())))
        var best = 1.7976931348623157e308
        var best_idx = 0
        var i = 0
        for b in self.boxes:
            var gd = _distance(self.geoms[i], _target)
            if gd < best:
                best = gd
                best_idx = i
                if best == 0.0:
                    return self.geoms[best_idx].copy()
            i += 1
        return self.geoms[best_idx].copy()

    fn query_knn(self, _target: Geometry, k: Int32) -> List[Geometry]:
        var out = List[Geometry]()
        if self.boxes.__len__() == 0 or k <= 0:
            return out.copy()
        var idxs = List[Int32]()
        var dists = List[Float64]()
        var i = 0
        for b in self.boxes:
            var gd = _distance(self.geoms[i], _target)
            var d2 = gd * gd
            idxs.append(Int32(i))
            dists.append(d2)
            i += 1
        # selection of k smallest
        var sel = List[Int32]()
        var taken = 0
        while taken < Int(k) and taken < idxs.__len__():
            var mi = taken
            var j = taken + 1
            while j < dists.__len__():
                if dists[j] < dists[mi]:
                    mi = j
                j += 1
            # swap into position
            var td = dists[taken]
            dists[taken] = dists[mi]
            dists[mi] = td
            var ti = idxs[taken]
            idxs[taken] = idxs[mi]
            idxs[mi] = ti
            sel.append(idxs[taken])
            taken += 1
        var r = 0
        while r < sel.__len__():
            out.append(self.geoms[sel[r]].copy())
            r += 1
        return out.copy()

    fn _nearest_idx(self, _target: Geometry) -> (Int32, Float64):
        if self.boxes.__len__() == 0:
            return (-1, 1.7976931348623157e308)
        var best = 1.7976931348623157e308
        var best_idx: Int32 = -1
        var i = 0
        for b in self.boxes:
            var gd = _distance(self.geoms[i], _target)
            if gd < best:
                best = gd
                best_idx = Int32(i)
                if best == 0.0:
                    return (best_idx, best)
            i += 1
        return (best_idx, best)

    fn _nearest_idx(self, _target: Point) -> (Int32, Float64):
        if self.boxes.__len__() == 0:
            return (-1, 1.7976931348623157e308)
        var best = 1.7976931348623157e308
        var best_idx: Int32 = -1
        var i = 0
        for b in self.boxes:
            var gd = _distance(self.geoms[i], _target)
            if gd < best:
                best = gd
                best_idx = Int32(i)
                if best == 0.0:
                    return (best_idx, best)
            i += 1
        return (best_idx, best)

    fn _query_indices(self, _target: Geometry) -> List[Int32]:
        var tb = _target.bounds()
        var out = List[Int32]()
        var i = 0
        for b in self.boxes:
            if self._env_intersects(b, tb):
                out.append(Int32(i))
            i += 1
        return out.copy()

    fn _query_indices(
        self, _target: Geometry, predicate: String
    ) -> List[Int32]:
        var tb = _target.bounds()
        var out = List[Int32]()
        var i = 0
        for b in self.boxes:
            if self._env_intersects(b, tb):
                if predicate == "intersects" and _intersects(
                    self.geoms[i], _target
                ):
                    out.append(Int32(i))
                elif predicate == "touches" and _touches(
                    self.geoms[i], _target
                ):
                    out.append(Int32(i))
                elif predicate == "overlaps" and _overlaps(
                    self.geoms[i], _target
                ):
                    out.append(Int32(i))
                elif predicate == "contains" and _contains(
                    self.geoms[i], _target
                ):
                    out.append(Int32(i))
                elif predicate == "within" and _contains(
                    _target, self.geoms[i]
                ):
                    out.append(Int32(i))
                elif predicate == "covers" and _covers(self.geoms[i], _target):
                    out.append(Int32(i))
                elif predicate == "covered_by" and _covers(
                    _target, self.geoms[i]
                ):
                    out.append(Int32(i))
                elif predicate == "contains_properly" and _contains_properly(
                    self.geoms[i], _target
                ):
                    out.append(Int32(i))
                elif (
                    predicate != "intersects"
                    and predicate != "touches"
                    and predicate != "overlaps"
                ):
                    out.append(Int32(i))
            i += 1
        return out.copy()

    fn nearest_item(self, _target: Geometry) -> Int32:
        var t = self._nearest_idx(_target)
        var idx = t[0]
        return idx

    fn nearest_item(self, _target: Point) -> Int32:
        var t = self._nearest_idx(_target)
        var idx = t[0]
        return idx

    fn query_bulk(self, _targets: List[Geometry]) -> List[Tuple[Int32, Int32]]:
        var pairs = List[Tuple[Int32, Int32]]()
        var ti = 0
        while ti < _targets.__len__():
            var idxs = self._query_indices(_targets[ti])
            var gi = 0
            while gi < idxs.__len__():
                pairs.append((Int32(ti), idxs[gi]))
                gi += 1
            ti += 1
        return pairs.copy()

    fn query_bulk(
        self, _targets: List[Geometry], predicate: String
    ) -> List[Tuple[Int32, Int32]]:
        var pairs = List[Tuple[Int32, Int32]]()
        var ti = 0
        while ti < _targets.__len__():
            var idxs = self._query_indices(_targets[ti], predicate)
            var gi = 0
            while gi < idxs.__len__():
                pairs.append((Int32(ti), idxs[gi]))
                gi += 1
            ti += 1
        return pairs.copy()

    fn nearest_all(self, _targets: List[Geometry]) -> List[Tuple[Int32, Int32]]:
        var out = List[Tuple[Int32, Int32]]()
        var i = 0
        while i < _targets.__len__():
            var t = self._nearest_idx(_targets[i])
            if t[0] != -1:
                out.append((Int32(i), t[0]))
            i += 1
        return out.copy()

    fn nearest_all(
        self,
        _targets: List[Geometry],
        max_distance: Float64,
        return_distance: Bool,
    ) -> List[Tuple[Int32, Int32, Float64]]:
        var out = List[Tuple[Int32, Int32, Float64]]()
        var i = 0
        while i < _targets.__len__():
            var t = self._nearest_idx(_targets[i])
            var idx = t[0]
            var d = t[1]
            if idx != -1:
                if max_distance < 0.0 or d <= max_distance:
                    out.append((Int32(i), idx, d))
            i += 1
        return out.copy()
