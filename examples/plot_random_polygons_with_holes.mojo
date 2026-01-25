from python import Python, PythonObject

from shapely._geometry import Geometry
from shapely.geometry import Point, MultiPoint, LinearRing, Polygon
from shapely.constructive import convex_hull, circle
from shapely.algorithms import point_in_polygon


fn _ensure_outputs_dir() raises:
    var os: PythonObject = Python.import_module("os")
    os.makedirs("outputs", exist_ok=True)


struct Rand(Copyable, Movable):
    var state: UInt64

    fn __init__(out self, seed: UInt64):
        self.state = seed

    fn next_u64(mut self) -> UInt64:
        # xorshift64*
        var x = self.state
        x ^= x >> 12
        x ^= x << 25
        x ^= x >> 27
        self.state = x
        return x * 0x2545F4914F6CDD1D

    fn next_f64(mut self) -> Float64:
        # [0, 1)
        var x = self.next_u64()
        var mant = Float64(Int64(x & 0x000F_FFFF_FFFF_FFFF))
        return mant / Float64(Int64(0x0010_0000_0000_0000))

    fn uniform(mut self, lo: Float64, hi: Float64) -> Float64:
        return lo + (hi - lo) * self.next_f64()

    fn randint(mut self, lo: Int32, hi: Int32) -> Int32:
        # inclusive
        var span = UInt64(Int64(hi - lo + 1))
        var v = self.next_u64() % span
        return lo + Int32(Int64(v))


fn _plot_coords(
    plt: PythonObject,
    coords: List[Tuple[Float64, Float64]],
    color: String,
    lw: Int = 2,
    closed: Bool = True,
    alpha: Float64 = 1.0,
) raises:
    var xs = Python.list()
    var ys = Python.list()
    for p in coords:
        xs.append(p[0])
        ys.append(p[1])

    if closed and coords.__len__() > 0:
        var first = coords[0]
        var last = coords[coords.__len__() - 1]
        if first[0] != last[0] or first[1] != last[1]:
            xs.append(first[0])
            ys.append(first[1])

    plt.plot(xs, ys, color=color, linewidth=lw, alpha=alpha)


fn _plot_polygon(plt: PythonObject, poly: Polygon) raises:
    _plot_coords(plt, poly.shell.coords, "tab:blue", lw=2, closed=True, alpha=0.9)
    for h in poly.holes:
        _plot_coords(plt, h.coords, "tab:red", lw=2, closed=True, alpha=0.9)


fn _make_random_hull(mut rng: Rand, npts: Int32) -> Polygon:
    var pts = List[Point]()
    var i: Int32 = 0
    while i < npts:
        var x = rng.uniform(0.0, 10.0)
        var y = rng.uniform(0.0, 10.0)
        pts.append(Point(x, y))
        i += 1

    var hull = convex_hull(Geometry(MultiPoint(pts)))
    if hull.is_polygon():
        return hull.as_polygon()

    # Fallback: guaranteed box
    var shell = List[Tuple[Float64, Float64]]()
    shell.append((0.0, 0.0))
    shell.append((10.0, 0.0))
    shell.append((10.0, 10.0))
    shell.append((0.0, 10.0))
    shell.append((0.0, 0.0))
    return Polygon(LinearRing(shell))


fn _try_place_hole(
    mut rng: Rand,
    shell_poly: Polygon,
    mut centers: List[Tuple[Float64, Float64, Float64]],
    r: Float64,
    quad_segs: Int32,
) -> (LinearRing, Bool):
    var b = shell_poly.bounds()
    var minx = b[0]
    var miny = b[1]
    var maxx = b[2]
    var maxy = b[3]

    var attempts: Int32 = 0
    while attempts < 200:
        var cx = rng.uniform(minx, maxx)
        var cy = rng.uniform(miny, maxy)

        if point_in_polygon(Point(cx, cy), shell_poly) != 1:
            attempts += 1
            continue

        # non-overlap with other holes
        var ok = True
        for c in centers:
            var dx = cx - c[0]
            var dy = cy - c[1]
            var rr = r + c[2] + 0.15
            if dx * dx + dy * dy <= rr * rr:
                ok = False
                break
        if not ok:
            attempts += 1
            continue

        # Ensure sampled circle points are inside the shell
        var disk = circle(cx, cy, r, quad_segs)
        var ring = disk.as_polygon().shell.copy()
        var j: Int32 = 0
        while j < ring.coords.__len__() - 1:
            var p = ring.coords[j]
            if point_in_polygon(Point(p[0], p[1]), shell_poly) == 0:
                ok = False
                break
            j += 1
        if not ok:
            attempts += 1
            continue

        centers.append((cx, cy, r))
        return (ring.copy(), True)

    return (LinearRing(List[Tuple[Float64, Float64]]()), False)


fn _make_random_polygon_with_holes(mut rng: Rand) -> Polygon:
    var shell_poly = _make_random_hull(rng, 25)

    var holes = List[LinearRing]()
    var centers = List[Tuple[Float64, Float64, Float64]]()

    var nholes = rng.randint(1, 4)
    var hi: Int32 = 0
    while hi < nholes:
        var r = rng.uniform(0.35, 1.1)
        var h = _try_place_hole(rng, shell_poly, centers, r, 16)
        if h[1]:
            holes.append(h[0].copy())
            hi += 1
        else:
            # give up early if we can't place more
            break

    return Polygon(shell_poly.shell.copy(), holes)


fn main() raises:
    _ensure_outputs_dir()

    var plt: PythonObject = Python.import_module("matplotlib.pyplot")
    var builtins: PythonObject = Python.import_module("builtins")

    var figsize_list = Python.list()
    figsize_list.append(11)
    figsize_list.append(7)
    var figsize: PythonObject = builtins.tuple(figsize_list)

    var fig: PythonObject = plt.figure(figsize=figsize)

    var seed: UInt64 = 123456789
    var rng = Rand(seed)

    var cols: Int32 = 3
    var rows: Int32 = 2

    var idx: Int32 = 1
    var i: Int32 = 0
    while i < rows * cols:
        var ax = fig.add_subplot(rows, cols, idx)
        var p = _make_random_polygon_with_holes(rng)

        _plot_polygon(plt, p)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(True, linewidth=0.4, alpha=0.4)
        idx += 1
        i += 1

    fig.tight_layout()
    fig.savefig("outputs/random_polygons_with_holes.png", dpi=160)
    print("wrote outputs/random_polygons_with_holes.png")
