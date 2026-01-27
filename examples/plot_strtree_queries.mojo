from python import Python, PythonObject

from shapely._geometry import Geometry
from shapely.creation import box
from shapely.geometry import Point
from shapely.strtree import STRtree


fn _ensure_outputs_dir() raises:
    var os: PythonObject = Python.import_module("os")
    os.makedirs("outputs", exist_ok=True)


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


fn _plot_geom(
    plt: PythonObject,
    geom: Geometry,
    shell_color: String,
    hole_color: String,
    lw: Int = 2,
    alpha: Float64 = 0.9,
) raises:
    if geom.is_polygon():
        var p = geom.as_polygon()
        _plot_coords(plt, p.shell.coords, shell_color, lw=lw, closed=True, alpha=alpha)
        for h in p.holes:
            _plot_coords(plt, h.coords, hole_color, lw=lw, closed=True, alpha=alpha)
    elif geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        for p in mp.polys:
            _plot_coords(plt, p.shell.coords, shell_color, lw=lw, closed=True, alpha=alpha)
            for h in p.holes:
                _plot_coords(plt, h.coords, hole_color, lw=lw, closed=True, alpha=alpha)
    elif geom.is_geometrycollection():
        var gc = geom.as_geometrycollection()
        for g in gc.geoms:
            _plot_geom(plt, g.copy(), shell_color, hole_color, lw=lw, alpha=alpha)


fn main() raises:
    _ensure_outputs_dir()

    var plt: PythonObject = Python.import_module("matplotlib.pyplot")
    var builtins: PythonObject = Python.import_module("builtins")

    # Build a small set of polygons laid out in a grid.
    var geoms = List[Geometry]()
    var x = 0
    while x < 4:
        var y = 0
        while y < 3:
            var xmin = Float64(x) * 3.0
            var ymin = Float64(y) * 3.0
            var xmax = xmin + 2.0
            var ymax = ymin + 2.0
            geoms.append(Geometry(box(xmin, ymin, xmax, ymax)))
            y += 1
        x += 1

    var tree = STRtree(geoms)

    # Query geometry and query point
    var query_poly = Geometry(box(2.5, 1.5, 7.0, 4.5))
    var query_pt = Point(8.2, 1.2)

    var hits = tree.query(query_poly.copy(), "intersects")
    var nearest = tree.nearest(query_pt)
    var knn = tree.query_knn(Geometry(query_pt.copy()), 3)

    print("STRtree size:", geoms.__len__())
    print("query intersects count:", hits.__len__())
    print("kNN count:", knn.__len__())

    # Global plot bounds for consistent scale.
    var minx = 1.0e308
    var miny = 1.0e308
    var maxx = -1.0e308
    var maxy = -1.0e308

    var bounds_geoms = List[Geometry]()
    for g in geoms:
        bounds_geoms.append(g.copy())
    bounds_geoms.append(query_poly.copy())

    for g in bounds_geoms:
        if g.is_empty():
            continue
        var b = g.bounds()
        if b[0] < minx:
            minx = b[0]
        if b[1] < miny:
            miny = b[1]
        if b[2] > maxx:
            maxx = b[2]
        if b[3] > maxy:
            maxy = b[3]

    # Include query point in bounds
    if query_pt.x < minx:
        minx = query_pt.x
    if query_pt.y < miny:
        miny = query_pt.y
    if query_pt.x > maxx:
        maxx = query_pt.x
    if query_pt.y > maxy:
        maxy = query_pt.y

    var pad = 0.5

    var figsize_list = Python.list()
    figsize_list.append(14)
    figsize_list.append(5)
    var figsize: PythonObject = builtins.tuple(figsize_list)

    var fig: PythonObject = plt.figure(figsize=figsize)

    # 1) Inputs overview
    var ax0 = fig.add_subplot(1, 3, 1)
    ax0.set_title("inputs")
    for g in geoms:
        _plot_geom(plt, g, "0.6", "0.6", lw=1, alpha=0.7)
    ax0.set_aspect("equal", adjustable="box")
    ax0.set_xlim(minx - pad, maxx + pad)
    ax0.set_ylim(miny - pad, maxy + pad)
    ax0.grid(True, linewidth=0.4, alpha=0.4)

    # 2) Predicate hits
    var ax1 = fig.add_subplot(1, 3, 2)
    ax1.set_title("query: intersects")
    for g in geoms:
        _plot_geom(plt, g, "0.85", "0.85", lw=1, alpha=0.5)
    _plot_geom(plt, query_poly, "black", "black", lw=2, alpha=0.95)
    for g in hits:
        _plot_geom(plt, g, "tab:blue", "tab:blue", lw=3, alpha=0.95)
    ax1.set_aspect("equal", adjustable="box")
    ax1.set_xlim(minx - pad, maxx + pad)
    ax1.set_ylim(miny - pad, maxy + pad)
    ax1.grid(True, linewidth=0.4, alpha=0.4)

    # 3) Nearest / kNN
    var ax2 = fig.add_subplot(1, 3, 3)
    ax2.set_title("nearest + kNN")
    for g in geoms:
        _plot_geom(plt, g, "0.85", "0.85", lw=1, alpha=0.5)
    for g in knn:
        _plot_geom(plt, g, "tab:green", "tab:green", lw=2, alpha=0.8)
    _plot_geom(plt, nearest, "tab:red", "tab:red", lw=3, alpha=0.95)
    ax2.scatter([query_pt.x], [query_pt.y], color="black", s=50)
    ax2.set_aspect("equal", adjustable="box")
    ax2.set_xlim(minx - pad, maxx + pad)
    ax2.set_ylim(miny - pad, maxy + pad)
    ax2.grid(True, linewidth=0.4, alpha=0.4)

    fig.tight_layout()
    fig.savefig("outputs/strtree_queries.png", dpi=160)
    print("wrote outputs/strtree_queries.png")
