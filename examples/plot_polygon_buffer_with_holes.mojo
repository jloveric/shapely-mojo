from python import Python, PythonObject

from shapely._geometry import Geometry
from shapely.geometry import Point, LinearRing, Polygon
from shapely.constructive import buffer, circle, JOIN_ROUND, JOIN_BEVEL, JOIN_MITRE
from shapely.validation import make_valid


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


fn _plot_polygon(
    plt: PythonObject,
    poly: Polygon,
    shell_color: String,
    hole_color: String,
    lw: Int = 2,
    alpha: Float64 = 0.9,
) raises:
    _plot_coords(plt, poly.shell.coords, shell_color, lw=lw, closed=True, alpha=alpha)
    for h in poly.holes:
        _plot_coords(plt, h.coords, hole_color, lw=lw, closed=True, alpha=alpha)


fn _plot_geom(
    plt: PythonObject,
    geom: Geometry,
    shell_color: String,
    hole_color: String,
    lw: Int = 2,
    alpha: Float64 = 0.9,
) raises:
    if geom.is_polygon():
        _plot_polygon(plt, geom.as_polygon(), shell_color, hole_color, lw=lw, alpha=alpha)
    elif geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        for p in mp.polys:
            _plot_polygon(plt, p.copy(), shell_color, hole_color, lw=lw, alpha=alpha)


fn _base_polygon_with_holes() -> Polygon:
    # Simple rectangle shell
    var shell = List[Tuple[Float64, Float64]]()
    shell.append((0.0, 0.0))
    shell.append((10.0, 0.0))
    shell.append((10.0, 7.0))
    shell.append((0.0, 7.0))
    shell.append((0.0, 0.0))

    # 3 circular holes
    var holes = List[LinearRing]()

    var h1 = circle(3.0, 2.0, 0.9, 16).as_polygon().shell.copy()
    var h2 = circle(7.0, 2.2, 1.1, 16).as_polygon().shell.copy()
    var h3 = circle(5.0, 5.0, 0.7, 16).as_polygon().shell.copy()

    holes.append(h1.copy())
    holes.append(h2.copy())
    holes.append(h3.copy())

    return Polygon(LinearRing(shell), holes)


fn main() raises:
    _ensure_outputs_dir()

    var plt: PythonObject = Python.import_module("matplotlib.pyplot")
    var builtins: PythonObject = Python.import_module("builtins")

    var figsize_list = Python.list()
    figsize_list.append(12)
    figsize_list.append(8)
    var figsize: PythonObject = builtins.tuple(figsize_list)
    var fig: PythonObject = plt.figure(figsize=figsize)

    var base_poly = _base_polygon_with_holes()
    var base_geom = Geometry(base_poly.copy())

    var joins = List[Tuple[String, Int32]]()
    joins.append(("round", Int32(JOIN_ROUND)))
    joins.append(("bevel", Int32(JOIN_BEVEL)))
    joins.append(("mitre", Int32(JOIN_MITRE)))

    var distances = List[Float64]()
    distances.append(0.6)
    distances.append(1.2)

    var rows = distances.__len__()
    var cols = joins.__len__()

    var r: Int = 0
    while r < rows:
        var c: Int = 0
        while c < cols:
            var idx = r * cols + c + 1
            var ax = fig.add_subplot(rows, cols, idx)
            ax.set_title("d=" + distances[r].__str__() + ", join=" + joins[c][0])

            var d = distances[r]
            var join_style = joins[c][1]

            var buf = buffer(base_geom.copy(), d, 16, 1, join_style, 5.0)
            var fixed = make_valid(buf.copy())

            _plot_geom(plt, base_geom, "black", "black", lw=2, alpha=0.85)
            _plot_geom(plt, buf, "tab:blue", "tab:red", lw=2, alpha=0.85)
            _plot_geom(plt, fixed, "tab:orange", "tab:orange", lw=2, alpha=0.6)

            ax.set_aspect("equal", adjustable="box")
            ax.grid(True, linewidth=0.4, alpha=0.4)

            c += 1
        r += 1

    fig.tight_layout()
    fig.savefig("outputs/polygon_buffer_with_holes.png", dpi=160)
    print("wrote outputs/polygon_buffer_with_holes.png")
