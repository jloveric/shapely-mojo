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
    var grid: List[List[Int32]]
    var nx: Int32
    var ny: Int32
    var minx: Float64
    var miny: Float64
    var maxx: Float64
    var maxy: Float64
    var cell_w: Float64
    var cell_h: Float64

    fn __init__(out self, geoms: List[Geometry]):
        self.geoms = geoms.copy()
        self.boxes = List[Tuple[Float64, Float64, Float64, Float64]]()
        self.grid = List[List[Int32]]()
        self.nx = 0
        self.ny = 0
        self.minx = 0.0
        self.miny = 0.0
        self.maxx = 0.0
        self.maxy = 0.0
        self.cell_w = 1.0
        self.cell_h = 1.0

        fn _bounds_of(g: Geometry) -> (Float64, Float64, Float64, Float64):
            return g.bounds()

        for g in self.geoms:
            self.boxes.append(_bounds_of(g))

        self._build_grid()

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

    fn _clamp_i32(self, v: Int32, lo: Int32, hi: Int32) -> Int32:
        if v < lo:
            return lo
        if v > hi:
            return hi
        return v

    fn _cell_index(self, ix: Int32, iy: Int32) -> Int:
        return Int(ix) + Int(iy) * Int(self.nx)

    fn _cell_coords_of(self, x: Float64, y: Float64) -> (Int32, Int32):
        if self.nx <= 0 or self.ny <= 0:
            return (0, 0)
        var fx = (x - self.minx) / self.cell_w
        var fy = (y - self.miny) / self.cell_h
        var ix = Int32(Int(fx))
        var iy = Int32(Int(fy))
        ix = self._clamp_i32(ix, 0, self.nx - 1)
        iy = self._clamp_i32(iy, 0, self.ny - 1)
        return (ix, iy)

    fn _cell_range_for_bbox(
        self, b: Tuple[Float64, Float64, Float64, Float64]
    ) -> (Int32, Int32, Int32, Int32):
        if self.nx <= 0 or self.ny <= 0:
            return (0, 0, -1, -1)
        var fx0 = (b[0] - self.minx) / self.cell_w
        var fy0 = (b[1] - self.miny) / self.cell_h
        var fx1 = (b[2] - self.minx) / self.cell_w
        var fy1 = (b[3] - self.miny) / self.cell_h
        var ix0 = Int32(Int(fx0))
        var iy0 = Int32(Int(fy0))
        var ix1 = Int32(Int(fx1))
        var iy1 = Int32(Int(fy1))
        ix0 = self._clamp_i32(ix0, 0, self.nx - 1)
        iy0 = self._clamp_i32(iy0, 0, self.ny - 1)
        ix1 = self._clamp_i32(ix1, 0, self.nx - 1)
        iy1 = self._clamp_i32(iy1, 0, self.ny - 1)
        return (ix0, iy0, ix1, iy1)

    fn _build_grid(mut self):
        self.grid = List[List[Int32]]()
        self.nx = 0
        self.ny = 0
        if self.boxes.__len__() == 0:
            return

        var minx = 1.7976931348623157e308
        var miny = 1.7976931348623157e308
        var maxx = -1.7976931348623157e308
        var maxy = -1.7976931348623157e308
        for b in self.boxes:
            if b[0] < minx:
                minx = b[0]
            if b[1] < miny:
                miny = b[1]
            if b[2] > maxx:
                maxx = b[2]
            if b[3] > maxy:
                maxy = b[3]

        self.minx = minx
        self.miny = miny
        self.maxx = maxx
        self.maxy = maxy

        var n = self.boxes.__len__()
        var side = Int32(Int(self._sqrt_f64(Float64(n))))
        if side < 1:
            side = 1
        if Int(side) * Int(side) < n:
            side += 1

        self.nx = side
        self.ny = side

        var w = self.maxx - self.minx
        var h = self.maxy - self.miny
        if w <= 0.0:
            w = 1.0
        if h <= 0.0:
            h = 1.0

        self.cell_w = w / Float64(self.nx)
        self.cell_h = h / Float64(self.ny)
        if self.cell_w <= 0.0:
            self.cell_w = 1.0
        if self.cell_h <= 0.0:
            self.cell_h = 1.0

        var cell_count = Int(self.nx) * Int(self.ny)
        var ci = 0
        while ci < cell_count:
            self.grid.append(List[Int32]())
            ci += 1

        var i = 0
        while i < self.boxes.__len__():
            var b = self.boxes[i]
            var cr = self._cell_range_for_bbox(b)
            var ix = cr[0]
            while ix <= cr[2]:
                var iy = cr[1]
                while iy <= cr[3]:
                    self.grid[self._cell_index(ix, iy)].append(Int32(i))
                    iy += 1
                ix += 1
            i += 1

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
        var t = self._nearest_idx(_target)
        if t[0] == -1:
            return Geometry(Polygon(LinearRing(List[Tuple[Float64, Float64]]())))
        return self.geoms[Int(t[0])].copy()

    fn nearest(self, _target: Point) -> Geometry:
        if self.boxes.__len__() == 0:
            return Geometry(Polygon(LinearRing(List[Tuple[Float64, Float64]]())))
        var t = self._nearest_idx(_target)
        if t[0] == -1:
            return Geometry(Polygon(LinearRing(List[Tuple[Float64, Float64]]())))
        return self.geoms[Int(t[0])].copy()

    fn query_knn(self, _target: Geometry, k: Int32) -> List[Geometry]:
        var out = List[Geometry]()
        if self.boxes.__len__() == 0 or k <= 0:
            return out.copy()
        var sel = self._knn_indices(_target, k)
        var r = 0
        while r < sel.__len__():
            out.append(self.geoms[Int(sel[r])].copy())
            r += 1
        return out.copy()

    fn _nearest_idx(self, _target: Geometry) -> (Int32, Float64):
        if self.boxes.__len__() == 0:
            return (-1, 1.7976931348623157e308)
        if self.nx <= 0 or self.ny <= 0:
            return self._nearest_idx_fallback(_target)

        var tb = _target.bounds()
        var cx = 0.5 * (tb[0] + tb[2])
        var cy = 0.5 * (tb[1] + tb[3])
        var c = self._cell_coords_of(cx, cy)
        var ix0 = c[0]
        var iy0 = c[1]

        var best = 1.7976931348623157e308
        var best_idx: Int32 = -1
        var ring: Int32 = 0
        var max_ring = self.nx
        if self.ny > max_ring:
            max_ring = self.ny

        while ring <= max_ring:
            var ix = ix0 - ring
            while ix <= ix0 + ring:
                if ix < 0 or ix >= self.nx:
                    ix += 1
                    continue
                var iy = iy0 - ring
                while iy <= iy0 + ring:
                    if iy < 0 or iy >= self.ny:
                        iy += 1
                        continue
                    if ring != 0 and ix != ix0 - ring and ix != ix0 + ring and iy != iy0 - ring and iy != iy0 + ring:
                        iy += 1
                        continue
                    var cidx = self._cell_index(ix, iy)
                    var cj = 0
                    while cj < self.grid[cidx].__len__():
                        var gi = self.grid[cidx][cj]
                        var tgt = _target.copy()
                        var d = _distance(self.geoms[Int(gi)], tgt)
                        if d < best:
                            best = d
                            best_idx = gi
                            if best == 0.0:
                                return (best_idx, best)
                        cj += 1
                    iy += 1
                ix += 1

            if best_idx != -1:
                var bound = Float64(ring) * self.cell_w
                if self.cell_h < self.cell_w:
                    bound = Float64(ring) * self.cell_h
                if bound > best:
                    break
            ring += 1

        return (best_idx, best)

    fn _nearest_idx_fallback(self, _target: Geometry) -> (Int32, Float64):
        var best = 1.7976931348623157e308
        var best_idx: Int32 = -1
        var i = 0
        while i < self.boxes.__len__():
            ref g = self.geoms[i]
            var tgt = _target.copy()
            var gd = _distance(g, tgt)
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
        if self.nx <= 0 or self.ny <= 0:
            return self._nearest_idx_fallback_point(_target)

        var c = self._cell_coords_of(_target.x, _target.y)
        var ix0 = c[0]
        var iy0 = c[1]

        var best = 1.7976931348623157e308
        var best_idx: Int32 = -1
        var ring: Int32 = 0
        var max_ring = self.nx
        if self.ny > max_ring:
            max_ring = self.ny

        while ring <= max_ring:
            var ix = ix0 - ring
            while ix <= ix0 + ring:
                if ix < 0 or ix >= self.nx:
                    ix += 1
                    continue
                var iy = iy0 - ring
                while iy <= iy0 + ring:
                    if iy < 0 or iy >= self.ny:
                        iy += 1
                        continue
                    if ring != 0 and ix != ix0 - ring and ix != ix0 + ring and iy != iy0 - ring and iy != iy0 + ring:
                        iy += 1
                        continue
                    var cidx = self._cell_index(ix, iy)
                    var cj = 0
                    while cj < self.grid[cidx].__len__():
                        var gi = self.grid[cidx][cj]
                        var pt = _target.copy()
                        var d = _distance(self.geoms[Int(gi)], pt)
                        if d < best:
                            best = d
                            best_idx = gi
                            if best == 0.0:
                                return (best_idx, best)
                        cj += 1
                    iy += 1
                ix += 1

            if best_idx != -1:
                var bound = Float64(ring) * self.cell_w
                if self.cell_h < self.cell_w:
                    bound = Float64(ring) * self.cell_h
                if bound > best:
                    break
            ring += 1

        return (best_idx, best)

    fn _nearest_idx_fallback_point(self, _target: Point) -> (Int32, Float64):
        var best = 1.7976931348623157e308
        var best_idx: Int32 = -1
        var i = 0
        while i < self.boxes.__len__():
            ref g = self.geoms[i]
            var pt = _target.copy()
            var gd = _distance(g, pt)
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
        if self.nx <= 0 or self.ny <= 0:
            var i = 0
            for b in self.boxes:
                if self._env_intersects(b, tb):
                    out.append(Int32(i))
                i += 1
            return out.copy()

        var cr = self._cell_range_for_bbox(tb)
        var n = self.boxes.__len__()
        var seen = List[Bool]()
        var si = 0
        while si < n:
            seen.append(False)
            si += 1

        var ix = cr[0]
        while ix <= cr[2]:
            var iy = cr[1]
            while iy <= cr[3]:
                var cidx = self._cell_index(ix, iy)
                var cj = 0
                while cj < self.grid[cidx].__len__():
                    var gi = self.grid[cidx][cj]
                    var ii = Int(gi)
                    if not seen[ii]:
                        seen[ii] = True
                        if self._env_intersects(self.boxes[ii], tb):
                            out.append(gi)
                    cj += 1
                iy += 1
            ix += 1
        return out.copy()

    fn query_items(self, _target: Geometry) -> List[Int32]:
        return self._query_indices(_target)

    fn query_items(self, _target: Geometry, predicate: String) -> List[Int32]:
        return self._query_indices(_target, predicate)

    fn _query_indices(
        self, _target: Geometry, predicate: String
    ) -> List[Int32]:
        var out = List[Int32]()
        var idxs = self._query_indices(_target)
        var j = 0
        while j < idxs.__len__():
            var i = Int(idxs[j])
            if predicate == "intersects" and _intersects(self.geoms[i], _target):
                out.append(idxs[j])
            elif predicate == "touches" and _touches(self.geoms[i], _target):
                out.append(idxs[j])
            elif predicate == "overlaps" and _overlaps(self.geoms[i], _target):
                out.append(idxs[j])
            elif predicate == "contains" and _contains(self.geoms[i], _target):
                out.append(idxs[j])
            elif predicate == "within" and _contains(_target, self.geoms[i]):
                out.append(idxs[j])
            elif predicate == "covers" and _covers(self.geoms[i], _target):
                out.append(idxs[j])
            elif predicate == "covered_by" and _covers(_target, self.geoms[i]):
                out.append(idxs[j])
            elif predicate == "contains_properly" and _contains_properly(
                self.geoms[i], _target
            ):
                out.append(idxs[j])
            elif (
                predicate != "intersects"
                and predicate != "touches"
                and predicate != "overlaps"
            ):
                out.append(idxs[j])
            j += 1
        return out.copy()

    fn _knn_indices(self, _target: Geometry, k: Int32) -> List[Int32]:
        var out = List[Int32]()
        if self.boxes.__len__() == 0 or k <= 0:
            return out.copy()
        if self.nx <= 0 or self.ny <= 0:
            var idxs = List[Int32]()
            var dists = List[Float64]()
            var i = 0
            while i < self.boxes.__len__():
                var tgt = _target.copy()
                var gd = _distance(self.geoms[i], tgt)
                idxs.append(Int32(i))
                dists.append(gd * gd)
                i += 1
            var sel = List[Int32]()
            var taken = 0
            while taken < Int(k) and taken < idxs.__len__():
                var mi = taken
                var j = taken + 1
                while j < dists.__len__():
                    if dists[j] < dists[mi]:
                        mi = j
                    j += 1
                var td = dists[taken]
                dists[taken] = dists[mi]
                dists[mi] = td
                var ti = idxs[taken]
                idxs[taken] = idxs[mi]
                idxs[mi] = ti
                sel.append(idxs[taken])
                taken += 1
            return sel.copy()

        var tb = _target.bounds()
        var cx = 0.5 * (tb[0] + tb[2])
        var cy = 0.5 * (tb[1] + tb[3])
        var c = self._cell_coords_of(cx, cy)
        var ix0 = c[0]
        var iy0 = c[1]

        var best_ids = List[Int32]()
        var best_d2 = List[Float64]()
        var thresh2 = 1.7976931348623157e308

        var ring: Int32 = 0
        var max_ring = self.nx
        if self.ny > max_ring:
            max_ring = self.ny

        while ring <= max_ring:
            var ix = ix0 - ring
            while ix <= ix0 + ring:
                if ix < 0 or ix >= self.nx:
                    ix += 1
                    continue
                var iy = iy0 - ring
                while iy <= iy0 + ring:
                    if iy < 0 or iy >= self.ny:
                        iy += 1
                        continue
                    if ring != 0 and ix != ix0 - ring and ix != ix0 + ring and iy != iy0 - ring and iy != iy0 + ring:
                        iy += 1
                        continue
                    var cidx = self._cell_index(ix, iy)
                    var cj = 0
                    while cj < self.grid[cidx].__len__():
                        var gi = self.grid[cidx][cj]
                        var tgt = _target.copy()
                        var gd = _distance(self.geoms[Int(gi)], tgt)
                        var d2 = gd * gd
                        if best_ids.__len__() < Int(k) or d2 < thresh2:
                            var pos = 0
                            while pos < best_d2.__len__() and best_d2[pos] <= d2:
                                pos += 1
                            best_ids.append(gi)
                            best_d2.append(d2)
                            var s = best_d2.__len__() - 1
                            while s > pos:
                                best_d2[s] = best_d2[s - 1]
                                best_ids[s] = best_ids[s - 1]
                                s -= 1
                            best_d2[pos] = d2
                            best_ids[pos] = gi
                            if best_ids.__len__() > Int(k):
                                var tmp_ids = List[Int32]()
                                var tmp_d2 = List[Float64]()
                                var ti2 = 0
                                while ti2 < Int(k):
                                    tmp_ids.append(best_ids[ti2])
                                    tmp_d2.append(best_d2[ti2])
                                    ti2 += 1
                                best_ids = tmp_ids.copy()
                                best_d2 = tmp_d2.copy()
                            if best_ids.__len__() == Int(k):
                                thresh2 = best_d2[best_d2.__len__() - 1]
                        cj += 1
                    iy += 1
                ix += 1

            if best_ids.__len__() == Int(k):
                var bound = Float64(ring) * self.cell_w
                if self.cell_h < self.cell_w:
                    bound = Float64(ring) * self.cell_h
                if bound * bound > thresh2:
                    break

            ring += 1

        for gi in best_ids:
            out.append(gi)
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
