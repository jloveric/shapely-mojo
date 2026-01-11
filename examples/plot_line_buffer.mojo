from python import Python, PythonObject

from shapely._geometry import Geometry
from shapely.geometry import LineString, MultiLineString
from shapely.constructive import buffer


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


fn main() raises:
    _ensure_outputs_dir()

    var plt: PythonObject = Python.import_module("matplotlib.pyplot")
    var builtins: PythonObject = Python.import_module("builtins")

    var l1 = LineString([(0.0, 0.0), (2.0, 0.5), (4.0, 0.0)])
    var l2 = LineString([(0.0, 2.0), (1.0, 3.0), (2.0, 2.5), (4.0, 3.0)])

    var mls = MultiLineString([l1.copy(), l2.copy()])
    var buf = buffer(Geometry(mls.copy()), 0.25)

    var figsize_list = Python.list()
    figsize_list.append(7)
    figsize_list.append(5)
    var figsize: PythonObject = builtins.tuple(figsize_list)
    var fig: PythonObject = plt.figure(figsize=figsize)
    var ax = fig.add_subplot(1, 1, 1)

    var i = 0
    while i < mls.lines.__len__():
        var ln = mls.lines[i].copy()
        _plot_coords(plt, ln.coords, "black", lw=2, label="input line" if i == 0 else "")
        i += 1

    if buf.is_polygon():
        var poly = buf.as_polygon()
        _plot_coords(plt, poly.shell.coords, "tab:blue", lw=2, label="buffer", closed=True, alpha=0.8)
    elif buf.is_multipolygon():
        var mp = buf.as_multipolygon()
        var j = 0
        while j < mp.polys.__len__():
            var poly2 = mp.polys[j].copy()
            _plot_coords(plt, poly2.shell.coords, "tab:blue", lw=2, label="buffer" if j == 0 else "", closed=True, alpha=0.8)
            j += 1

    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.5, alpha=0.5)
    ax.legend()

    fig.tight_layout()
    fig.savefig("outputs/line_buffer.png", dpi=160)
    print("wrote outputs/line_buffer.png")
