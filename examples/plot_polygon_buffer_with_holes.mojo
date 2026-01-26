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
    var s: Float64 = 2.0
    var hole_scale: Float64 = 1.5

    # Simple rectangle shell
    var shell = List[Tuple[Float64, Float64]]()
    shell.append((0.0 * s, 0.0 * s))
    shell.append((10.0 * s, 0.0 * s))
    shell.append((10.0 * s, 7.0 * s))
    shell.append((0.0 * s, 7.0 * s))
    shell.append((0.0 * s, 0.0 * s))

    # 1 circular hole + 2 polygonal holes
    var holes = List[LinearRing]()

    # Circle hole
    var h1 = circle(3.0 * s, 2.0 * s, 0.75 * s * hole_scale, 16).as_polygon().shell.copy()

    # Diamond (polygon) hole
    var h2c = List[Tuple[Float64, Float64]]()
    var h2cx: Float64 = 7.0 * s
    var h2cy: Float64 = 2.2 * s
    h2c.append((h2cx + (7.0 * s - h2cx) * hole_scale, h2cy + (1.2 * s - h2cy) * hole_scale))
    h2c.append((h2cx + (8.0 * s - h2cx) * hole_scale, h2cy + (2.2 * s - h2cy) * hole_scale))
    h2c.append((h2cx + (7.0 * s - h2cx) * hole_scale, h2cy + (3.2 * s - h2cy) * hole_scale))
    h2c.append((h2cx + (6.0 * s - h2cx) * hole_scale, h2cy + (2.2 * s - h2cy) * hole_scale))
    h2c.append((h2cx + (7.0 * s - h2cx) * hole_scale, h2cy + (1.2 * s - h2cy) * hole_scale))
    var h2 = LinearRing(h2c)

    # Triangle (polygon) hole
    var h3c = List[Tuple[Float64, Float64]]()
    var h3cx: Float64 = 5.1 * s
    var h3cy: Float64 = 5.033333333333333 * s
    h3c.append((h3cx + (4.5 * s - h3cx) * hole_scale, h3cy + (4.6 * s - h3cy) * hole_scale))
    h3c.append((h3cx + (5.7 * s - h3cx) * hole_scale, h3cy + (4.7 * s - h3cy) * hole_scale))
    h3c.append((h3cx + (5.1 * s - h3cx) * hole_scale, h3cy + (5.8 * s - h3cy) * hole_scale))
    h3c.append((h3cx + (4.5 * s - h3cx) * hole_scale, h3cy + (4.6 * s - h3cy) * hole_scale))
    var h3 = LinearRing(h3c)

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
    distances.append(0.9)

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
