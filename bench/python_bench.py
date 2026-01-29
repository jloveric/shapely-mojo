import sys
from pathlib import Path
import time


def _ensure_real_shapely_import():
    # This repo has a top-level "shapely/" directory (vendored files) which can shadow
    # the PyPI shapely package if you run from the repo root.
    cwd = Path.cwd()

    # Remove cwd from sys.path if it contains a local shapely/ directory.
    if (cwd / "shapely").is_dir():
        sys.path = [p for p in sys.path if p not in ("", str(cwd))]

    import shapely  # noqa: E402

    shapely_file = getattr(shapely, "__file__", "")
    if "site-packages" not in shapely_file and "dist-packages" not in shapely_file:
        raise RuntimeError(
            "Imported a non-site-packages shapely (likely shadowed by repo checkout): "
            f"{shapely_file}"
        )

    return shapely


def _now_ns() -> int:
    return time.perf_counter_ns()


def _print_result(name: str, iters: int, elapsed_ns: int) -> None:
    secs = elapsed_ns / 1e9
    print(f"RESULT\t{name}\t{secs}\t{iters}")


def _run(name: str, warm: int, iters: int, fn) -> None:
    for _ in range(warm):
        fn()
    t0 = _now_ns()
    for _ in range(iters):
        fn()
    t1 = _now_ns()
    _print_result(name, iters, t1 - t0)


def main() -> None:
    shapely = _ensure_real_shapely_import()

    from shapely.geometry import box, Point, Polygon  # type: ignore
    from shapely.strtree import STRtree  # type: ignore

    # buffer
    p = box(0.0, 0.0, 2.0, 2.0)
    _run("buffer_box", warm=200, iters=2000, fn=lambda: p.buffer(0.25, quad_segs=8))

    conv = Polygon(
        [
            (0.0, 0.0),
            (2.0, 0.0),
            (3.0, 1.0),
            (2.0, 2.0),
            (0.0, 2.0),
            (-1.0, 1.0),
            (0.0, 0.0),
        ]
    )
    _run("buffer_convex_poly", warm=200, iters=2000, fn=lambda: conv.buffer(0.25, quad_segs=8))

    conc = Polygon(
        [
            (0.0, 0.0),
            (3.0, 0.0),
            (3.0, 3.0),
            (2.0, 3.0),
            (2.0, 1.0),
            (1.0, 1.0),
            (1.0, 3.0),
            (0.0, 3.0),
            (0.0, 0.0),
        ]
    )
    _run("buffer_concave_poly", warm=200, iters=2000, fn=lambda: conc.buffer(0.25, quad_segs=8))

    # boolean ops
    a = box(0.0, 0.0, 2.0, 2.0)
    b = box(1.0, 0.8, 3.0, 2.6)

    _run("union", warm=500, iters=5000, fn=lambda: a.union(b))
    _run("intersection", warm=500, iters=5000, fn=lambda: a.intersection(b))
    _run("difference", warm=500, iters=5000, fn=lambda: a.difference(b))
    _run(
        "symmetric_difference",
        warm=500,
        iters=5000,
        fn=lambda: a.symmetric_difference(b),
    )

    # STRtree: query + nearest
    items = []
    for x in range(40):
        for y in range(25):
            xmin = float(x) * 1.5
            ymin = float(y) * 1.5
            items.append(box(xmin, ymin, xmin + 1.0, ymin + 1.0))

    tree = STRtree(items)
    qpoly = box(10.0, 10.0, 30.0, 25.0)
    qpt = Point(31.2, 9.7)

    _run(
        "strtree_query_intersects_idx",
        warm=200,
        iters=2000,
        fn=lambda: tree.query(qpoly, predicate="intersects"),
    )

    _run(
        "strtree_query_intersects_geom",
        warm=200,
        iters=2000,
        fn=lambda: tree.geometries.take(tree.query(qpoly, predicate="intersects")),
    )

    _run(
        "strtree_nearest_point_idx",
        warm=200,
        iters=2000,
        fn=lambda: tree.nearest(qpt),
    )

    _run(
        "strtree_nearest_point_geom",
        warm=200,
        iters=2000,
        fn=lambda: tree.geometries.take(tree.nearest(qpt)),
    )

    # Print version last so the RESULT lines are easy to grep.
    print(f"INFO\tpython_shapely\t{getattr(shapely, '__version__', 'unknown')}")


if __name__ == "__main__":
    main()
