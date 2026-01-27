from python import Python, PythonObject

from shapely._geometry import Geometry
from shapely.creation import box
from shapely.set_operations import union, intersection, difference, symmetric_difference


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

    # Two overlapping squares/rectangles
    var a_poly = box(0.0, 0.0, 2.0, 2.0)
    var b_poly = box(1.0, 0.8, 3.0, 2.6)
    var a = Geometry(a_poly.copy())
    var b = Geometry(b_poly.copy())

    var g_inter = intersection(a.copy(), b.copy())
    var g_union = union(a.copy(), b.copy())
    var g_diff = difference(a.copy(), b.copy())
    var g_xor = symmetric_difference(a.copy(), b.copy())

    var minx = 1.0e308
    var miny = 1.0e308
    var maxx = -1.0e308
    var maxy = -1.0e308

    var bounds_geoms = List[Geometry]()
    bounds_geoms.append(a.copy())
    bounds_geoms.append(b.copy())
    bounds_geoms.append(g_inter.copy())
    bounds_geoms.append(g_union.copy())
    bounds_geoms.append(g_diff.copy())
    bounds_geoms.append(g_xor.copy())

    for gg in bounds_geoms:
        if gg.is_empty():
            continue
        var bb = gg.bounds()
        if bb[0] < minx:
            minx = bb[0]
        if bb[1] < miny:
            miny = bb[1]
        if bb[2] > maxx:
            maxx = bb[2]
        if bb[3] > maxy:
            maxy = bb[3]

    var pad = 0.15

    print("A area:", a.area())
    print("B area:", b.area())
    print("intersection area:", g_inter.area())
    print("union area:", g_union.area())
    print("difference (A-B) area:", g_diff.area())
    print("symmetric_difference area:", g_xor.area())

    var figsize_list = Python.list()
    figsize_list.append(12)
    figsize_list.append(8)
    var figsize: PythonObject = builtins.tuple(figsize_list)
    var fig: PythonObject = plt.figure(figsize=figsize)

    # 2 rows x 3 cols
    var titles = List[String]()
    titles.append("inputs")
    titles.append("intersection")
    titles.append("union")
    titles.append("A - B")
    titles.append("symmetric diff")
    titles.append("")

    var i = 0
    while i < 6:
        var ax = fig.add_subplot(2, 3, i + 1)
        ax.set_title(titles[i])

        if i == 0:
            _plot_geom(plt, a, "black", "black", lw=2, alpha=0.9)
            _plot_geom(plt, b, "tab:orange", "tab:orange", lw=2, alpha=0.9)
        elif i == 1:
            _plot_geom(plt, g_inter, "tab:blue", "tab:red", lw=2, alpha=0.9)
        elif i == 2:
            _plot_geom(plt, g_union, "tab:blue", "tab:red", lw=2, alpha=0.9)
        elif i == 3:
            _plot_geom(plt, g_diff, "tab:blue", "tab:red", lw=2, alpha=0.9)
        elif i == 4:
            _plot_geom(plt, g_xor, "tab:blue", "tab:red", lw=2, alpha=0.9)

        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(minx - pad, maxx + pad)
        ax.set_ylim(miny - pad, maxy + pad)
        ax.grid(True, linewidth=0.4, alpha=0.4)
        i += 1

    fig.tight_layout()
    fig.savefig("outputs/boolean_operations.png", dpi=160)
    print("wrote outputs/boolean_operations.png")
