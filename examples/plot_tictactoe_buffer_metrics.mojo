from python import Python, PythonObject

from shapely._geometry import Geometry
from shapely.geometry import LineString, MultiLineString
from shapely.constructive import buffer
from shapely.set_operations import unary_union


fn _ensure_outputs_dir() raises:
    var os: PythonObject = Python.import_module("os")
    os.makedirs("outputs", exist_ok=True)


fn _plot_coords(
    ax: PythonObject,
    coords: List[Tuple[Float64, Float64]],
    color: String,
    lw: Int = 2,
    closed: Bool = False,
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

    ax.plot(xs, ys, color=color, linewidth=lw, alpha=alpha)


fn _plot_geom(ax: PythonObject, geom: Geometry, color: String, lw: Int = 2, alpha: Float64 = 0.9) raises:
    if geom.is_polygon():
        var p = geom.as_polygon()
        _plot_coords(ax, p.shell.coords, color, lw=lw, closed=True, alpha=alpha)
        for h in p.holes:
            _plot_coords(ax, h.coords, color, lw=lw, closed=True, alpha=alpha)
    elif geom.is_multipolygon():
        var mp = geom.as_multipolygon()
        for p in mp.polys:
            _plot_coords(ax, p.shell.coords, color, lw=lw, closed=True, alpha=alpha)
            for h in p.holes:
                _plot_coords(ax, h.coords, color, lw=lw, closed=True, alpha=alpha)


fn _tictactoe_board(size: Float64 = 9.0) -> MultiLineString:
    # Build a 3x3 tic-tac-toe grid as 4 lines.
    # Coordinates are in [0, size] with lines at 1/3 and 2/3.
    var a = size / 3.0
    var b = 2.0 * size / 3.0

    var l0 = LineString([(a, 0.0), (a, size)])
    var l1 = LineString([(b, 0.0), (b, size)])
    var l2 = LineString([(0.0, a), (size, a)])
    var l3 = LineString([(0.0, b), (size, b)])

    return MultiLineString([l0.copy(), l1.copy(), l2.copy(), l3.copy()])


fn main() raises:
    _ensure_outputs_dir()

    var board = _tictactoe_board(9.0)
    var d = 0.35

    # Buffer each line then dissolve overlaps to get a single polygonal footprint.
    var pieces = List[Geometry]()
    for ln in board.lines:
        pieces.append(buffer(Geometry(ln.copy()), d, 16))

    var merged = unary_union(pieces)

    print("buffer distance:", d)
    print("area:", merged.area())
    print("perimeter:", merged.length())

    var plt: PythonObject = Python.import_module("matplotlib.pyplot")
    var fig: PythonObject = plt.figure()
    var ax = fig.add_subplot(1, 1, 1)

    for ln in board.lines:
        _plot_coords(ax, ln.coords, "black", lw=2, closed=False, alpha=0.9)
    _plot_geom(ax, merged, "tab:blue", lw=2, alpha=0.85)

    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, linewidth=0.4, alpha=0.4)
    ax.set_title("tictactoe buffer d=" + d.__str__() + " area=" + merged.area().__str__())

    fig.tight_layout()
    fig.savefig("outputs/tictactoe_buffer.png", dpi=160)
    print("wrote outputs/tictactoe_buffer.png")
