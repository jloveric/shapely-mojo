from python import Python, PythonObject

from shapely._geometry import Geometry
from shapely.geometry import LineString, MultiLineString
from shapely.constructive import (
    buffer,
    CAP_ROUND,
    CAP_FLAT,
    CAP_SQUARE,
    JOIN_ROUND,
    JOIN_BEVEL,
    JOIN_MITRE,
)
from shapely.validation import make_valid


fn _ensure_outputs_dir() raises:
    var os: PythonObject = Python.import_module("os")
    os.makedirs("outputs", exist_ok=True)


fn _plot_coords(
    plt: PythonObject,
    coords: List[Tuple[Float64, Float64]],
    color: String,
    lw: Int = 2,
    label: String = "",
    closed: Bool = False,
    alpha: Float64 = 1.0,
) raises:
    # coords: List[Tuple[Float64, Float64]]
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

    if label != "":
        plt.plot(xs, ys, color=color, linewidth=lw, alpha=alpha, label=label)
    else:
        plt.plot(xs, ys, color=color, linewidth=lw, alpha=alpha)


fn _plot_geom(
    plt: PythonObject,
    geom: Geometry,
    color: String,
    lw: Int = 2,
    alpha: Float64 = 0.8,
) raises:
    if geom.is_polygon():
        var poly = geom.as_polygon()
        _plot_coords(plt, poly.shell.coords, color, lw=lw, closed=True, alpha=alpha)
    elif geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        var j = 0
        while j < mp.polys.__len__():
            var p = mp.polys[j].copy()
            _plot_coords(plt, p.shell.coords, color, lw=lw, closed=True, alpha=alpha)
            j += 1


fn main() raises:
    _ensure_outputs_dir()

    var plt: PythonObject = Python.import_module("matplotlib.pyplot")
    var builtins: PythonObject = Python.import_module("builtins")

    var l1 = LineString([(0.0, 0.0), (2.0, 0.5), (4.0, 0.0)])
    var l2 = LineString([(0.0, 2.0), (1.0, 3.0), (2.0, 2.5), (4.0, 3.0)])

    # Default image (kept): round/round
    var buf0_a = buffer(Geometry(l1.copy()), 0.25)
    var buf0_b = buffer(Geometry(l2.copy()), 0.25)
    var figsize_list0 = Python.list()
    figsize_list0.append(7)
    figsize_list0.append(5)
    var figsize0: PythonObject = builtins.tuple(figsize_list0)
    var fig0: PythonObject = plt.figure(figsize=figsize0)
    var ax0 = fig0.add_subplot(1, 1, 1)
    _plot_coords(plt, l1.coords, "black", lw=2, label="input line")
    _plot_coords(plt, l2.coords, "black", lw=2)
    _plot_geom(plt, buf0_a, "tab:blue", lw=2, alpha=0.8)
    _plot_geom(plt, buf0_b, "tab:blue", lw=2, alpha=0.8)
    ax0.set_aspect("equal", adjustable="box")
    ax0.grid(True, linewidth=0.5, alpha=0.5)
    ax0.legend()
    fig0.tight_layout()
    fig0.savefig("outputs/line_buffer.png", dpi=160)

    # Styles grid
    var caps = List[Tuple[String, Int32]]()
    caps.append(("round", Int32(CAP_ROUND)))
    caps.append(("flat", Int32(CAP_FLAT)))
    caps.append(("square", Int32(CAP_SQUARE)))
    var joins = List[Tuple[String, Int32]]()
    joins.append(("round", Int32(JOIN_ROUND)))
    joins.append(("bevel", Int32(JOIN_BEVEL)))
    joins.append(("mitre", Int32(JOIN_MITRE)))

    var figsize_list = Python.list()
    figsize_list.append(11)
    figsize_list.append(9)
    var figsize: PythonObject = builtins.tuple(figsize_list)
    var fig: PythonObject = plt.figure(figsize=figsize)

    var r = 0
    while r < caps.__len__():
        var c = 0
        while c < joins.__len__():
            var idx = r * joins.__len__() + c + 1
            var ax = fig.add_subplot(caps.__len__(), joins.__len__(), idx)

            var cap_label = caps[r][0]
            var join_label = joins[c][0]
            ax.set_title(cap_label + "/" + join_label)

            var buf_a = buffer(Geometry(l1.copy()), 0.25, 8, caps[r][1], joins[c][1])
            var buf_b = buffer(Geometry(l2.copy()), 0.25, 8, caps[r][1], joins[c][1])
            var fixed_a = make_valid(buf_a.copy())
            var fixed_b = make_valid(buf_b.copy())

            _plot_coords(plt, l1.coords, "black", lw=1)
            _plot_coords(plt, l2.coords, "black", lw=1)

            _plot_geom(plt, buf_a, "tab:blue", lw=2, alpha=0.8)
            _plot_geom(plt, buf_b, "tab:blue", lw=2, alpha=0.8)

            _plot_geom(plt, fixed_a, "tab:orange", lw=2, alpha=0.8)
            _plot_geom(plt, fixed_b, "tab:orange", lw=2, alpha=0.8)

            ax.set_aspect("equal", adjustable="box")
            ax.grid(True, linewidth=0.4, alpha=0.4)
            c += 1
        r += 1

    fig.tight_layout()
    fig.savefig("outputs/line_buffer_styles.png", dpi=160)

    print("wrote outputs/line_buffer.png")
    print("wrote outputs/line_buffer_styles.png")
