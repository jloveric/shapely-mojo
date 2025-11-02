from shapely._geometry import Geometry
from shapely.geometry import LinearRing, Polygon
from shapely.measurement import distance as _distance
from shapely.predicates import intersects as _intersects, overlaps as _overlaps, touches as _touches, contains as _contains, within as _within, covers as _covers, covered_by as _covered_by, disjoint as _disjoint, crosses as _crosses, contains_properly as _contains_properly


struct STRtree:
    var geoms: List[Geometry]
    var boxes: List[Tuple[Float64, Float64, Float64, Float64]]  # (minx, miny, maxx, maxy)

    fn __init__(inout self, geoms: List[Geometry]):
        self.geoms = geoms
        self.boxes = List[Tuple[Float64, Float64, Float64, Float64]]()
        for g in geoms:
            self.boxes.push_back(g.bounds())

    fn _env_intersects(self, a: Tuple[Float64, Float64, Float64, Float64], b: Tuple[Float64, Float64, Float64, Float64]) -> Bool:
        return not (a[2] < b[0] or b[2] < a[0] or a[3] < b[1] or b[3] < a[1])

    fn _env_dist2(self, a: Tuple[Float64, Float64, Float64, Float64], b: Tuple[Float64, Float64, Float64, Float64]) -> Float64:
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
        if x <= 0.0: return 0.0
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

    fn _node_union_boxes(self, ids: List[Int32]) -> (Float64, Float64, Float64, Float64):
        # Stub retained for compatibility; not used in naive implementation
        return (0.0, 0.0, 0.0, 0.0)

    fn _geom_union_boxes(self, ids: List[Int32]) -> (Float64, Float64, Float64, Float64):
        var minx = 1.7976931348623157e308
        var miny = 1.7976931348623157e308
        var maxx = -1.7976931348623157e308
        var maxy = -1.7976931348623157e308
        var i = 0
        while i < ids.size():
            let b = self.boxes[ids[i]]
            if b[0] < minx: minx = b[0]
            if b[1] < miny: miny = b[1]
            if b[2] > maxx: maxx = b[2]
            if b[3] > maxy: maxy = b[3]
            i += 1
        return (minx, miny, maxx, maxy)

    fn _qsort_geom_by_minx(self, inout idx: List[Int32], lo: Int, hi: Int):
        var i = lo
        var j = hi
        let pivot = self.boxes[idx[(lo + hi) / 2]][0]
        while i <= j:
            while self.boxes[idx[i]][0] < pivot:
                i += 1
            while self.boxes[idx[j]][0] > pivot:
                j -= 1
            if i <= j:
                let tmp = idx[i]
                idx[i] = idx[j]
                idx[j] = tmp
                i += 1
                j -= 1
        if lo < j: self._qsort_geom_by_minx(idx, lo, j)
        if i < hi: self._qsort_geom_by_minx(idx, i, hi)

    fn _qsort_geom_slice_by_miny(self, inout idx: List[Int32], start: Int, end_excl: Int):
        fn sort_range(self_ref: STRtree, inout a: List[Int32], lo: Int, hi: Int):
            var i = lo
            var j = hi
            let pivot = self_ref.boxes[a[(lo + hi) / 2]][1]
            while i <= j:
                while self_ref.boxes[a[i]][1] < pivot: i += 1
                while self_ref.boxes[a[j]][1] > pivot: j -= 1
                if i <= j:
                    let tmp = a[i]
                    a[i] = a[j]
                    a[j] = tmp
                    i += 1
                    j -= 1
            if lo < j: sort_range(self_ref, a, lo, j)
            if i < hi: sort_range(self_ref, a, i, hi)
        if end_excl - start > 1:
            sort_range(self, idx, start, end_excl - 1)

    fn _build(inout self, max_children: Int32):
        # Naive implementation does not build a tree
        return

    fn query(self, _target: Geometry) -> List[Geometry]:
        let tb = _target.bounds()
        var out = List[Geometry]()
        var i = 0
        while i < self.geoms.size():
            if self._env_intersects(self.boxes[i], tb):
                out.push_back(self.geoms[i])
            i += 1
        return out


    fn query(self, _target: Geometry, predicate: String) -> List[Geometry]:
        var cands = self.query(_target)
        if predicate == "intersects":
            var out = List[Geometry]()
            for g in cands:
                if _intersects(g, _target): out.push_back(g)
            return out
        elif predicate == "touches":
            var out2 = List[Geometry]()
            for g in cands:
                if _touches(g, _target): out2.push_back(g)
            return out2
        elif predicate == "overlaps":
            var out3 = List[Geometry]()
            for g in cands:
                if _overlaps(g, _target): out3.push_back(g)
            return out3
        elif predicate == "contains":
            var out4 = List[Geometry]()
            for g in cands:
                if _contains(g, _target): out4.push_back(g)
            return out4
        elif predicate == "within":
            var out5 = List[Geometry]()
            for g in cands:
                if _within(g, _target): out5.push_back(g)
            return out5
        elif predicate == "covers":
            var out6 = List[Geometry]()
            for g in cands:
                if _covers(g, _target): out6.push_back(g)
            return out6
        elif predicate == "covered_by":
            var out7 = List[Geometry]()
            for g in cands:
                if _covered_by(g, _target): out7.push_back(g)
            return out7
        elif predicate == "disjoint":
            var out8 = List[Geometry]()
            for g in cands:
                if _disjoint(g, _target): out8.push_back(g)
            return out8
        elif predicate == "crosses":
            var out9 = List[Geometry]()
            for g in cands:
                if _crosses(g, _target): out9.push_back(g)
            return out9
        elif predicate == "contains_properly":
            var out10 = List[Geometry]()
            for g in cands:
                if _contains_properly(g, _target): out10.push_back(g)
            return out10
        return cands

    fn nearest(self, _target: Geometry) -> Geometry:
        if self.geoms.size() == 0:
            return Polygon(LinearRing(List[Tuple[Float64, Float64]]()))
        var best = 1.7976931348623157e308
        var best_idx = 0
        var i = 0
        while i < self.geoms.size():
            let gd = _distance(self.geoms[i], _target)
            let d2 = gd * gd
            if d2 < best:
                best = d2
                best_idx = i
                if best == 0.0:
                    return self.geoms[best_idx]
            i += 1
        return self.geoms[best_idx]


    fn query_knn(self, _target: Geometry, k: Int32) -> List[Geometry]:
        var out = List[Geometry]()
        if self.geoms.size() == 0 or k <= 0:
            return out
        var idxs = List[Int32]()
        var dists = List[Float64]()
        var i = 0
        while i < self.geoms.size():
            let gd = _distance(self.geoms[i], _target)
            let d2 = gd * gd
            idxs.push_back(i as Int32)
            dists.push_back(d2)
            i += 1
        # selection of k smallest
        var sel = List[Int32]()
        var taken = 0
        while taken < (k as Int) and taken < idxs.size():
            var mi = taken
            var j = taken + 1
            while j < dists.size():
                if dists[j] < dists[mi]: mi = j
                j += 1
            # swap into position
            let td = dists[taken]
            dists[taken] = dists[mi]
            dists[mi] = td
            let ti = idxs[taken]
            idxs[taken] = idxs[mi]
            idxs[mi] = ti
            sel.push_back(idxs[taken])
            taken += 1
        var r = 0
        while r < sel.size():
            out.push_back(self.geoms[sel[r]])
            r += 1
        return out


    fn _nearest_idx(self, _target: Geometry) -> (Int32, Float64):
        if self.geoms.size() == 0:
            return (-1, 1.7976931348623157e308)
        var best = 1.7976931348623157e308
        var best_idx: Int32 = 0
        var i = 0
        while i < self.geoms.size():
            let gd = _distance(self.geoms[i], _target)
            if gd < best:
                best = gd
                best_idx = i as Int32
                if best == 0.0:
                    return (best_idx, best)
            i += 1
        return (best_idx, best)


    fn _query_indices(self, _target: Geometry) -> List[Int32]:
        let tb = _target.bounds()
        var out = List[Int32]()
        var i = 0
        while i < self.boxes.size():
            if self._env_intersects(self.boxes[i], tb):
                out.push_back(i as Int32)
            i += 1
        return out


    fn _query_indices(self, _target: Geometry, predicate: String) -> List[Int32]:
        let tb = _target.bounds()
        var out = List[Int32]()
        var i = 0
        while i < self.boxes.size():
            if self._env_intersects(self.boxes[i], tb):
                if predicate == "intersects" and _intersects(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "touches" and _touches(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "overlaps" and _overlaps(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "contains" and _contains(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "within" and _within(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "covers" and _covers(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "covered_by" and _covered_by(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "disjoint" and _disjoint(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "crosses" and _crosses(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate == "contains_properly" and _contains_properly(self.geoms[i], _target):
                    out.push_back(i as Int32)
                elif predicate != "intersects" and predicate != "touches" and predicate != "overlaps":
                    out.push_back(i as Int32)
            i += 1
        return out


    fn nearest_item(self, _target: Geometry) -> Int32:
        let (idx, _d) = self._nearest_idx(_target)
        return idx


    fn query_bulk(self, _targets: List[Geometry]) -> List[Tuple[Int32, Int32]]:
        var pairs = List[Tuple[Int32, Int32]]()
        var ti = 0
        while ti < _targets.size():
            let idxs = self._query_indices(_targets[ti])
            var gi = 0
            while gi < idxs.size():
                pairs.push_back((ti as Int32, idxs[gi]))
                gi += 1
            ti += 1
        return pairs


    fn query_bulk(self, _targets: List[Geometry], predicate: String) -> List[Tuple[Int32, Int32]]:
        var pairs = List[Tuple[Int32, Int32]]()
        var ti = 0
        while ti < _targets.size():
            let idxs = self._query_indices(_targets[ti], predicate)
            var gi = 0
            while gi < idxs.size():
                pairs.push_back((ti as Int32, idxs[gi]))
                gi += 1
            ti += 1
        return pairs


    fn nearest_all(self, _targets: List[Geometry]) -> List[Tuple[Int32, Int32]]:
        var out = List[Tuple[Int32, Int32]]()
        var i = 0
        while i < _targets.size():
            let (idx, _d) = self._nearest_idx(_targets[i])
            if idx != -1:
                out.push_back((i as Int32, idx))
            i += 1
        return out


    fn nearest_all(self, _targets: List[Geometry], max_distance: Float64, return_distance: Bool) -> List[Tuple[Int32, Int32, Float64]]:
        var out = List[Tuple[Int32, Int32, Float64]]()
        var i = 0
        while i < _targets.size():
            let (idx, d) = self._nearest_idx(_targets[i])
            if idx != -1:
                if max_distance < 0.0 or d <= max_distance:
                    out.push_back((i as Int32, idx, d))
            i += 1
        return out
