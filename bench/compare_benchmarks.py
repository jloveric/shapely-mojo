import json
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple


RE_RESULT = re.compile(r"^RESULT\t(?P<name>[^\t]+)\t(?P<secs>[-+0-9.eE]+)\t(?P<iters>\d+)\s*$")


@dataclass
class BenchResult:
    secs: float
    iters: int


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _run_cmd(cmd: List[str], cwd: Path) -> str:
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            "Command failed\n"
            f"cmd: {' '.join(cmd)}\n"
            f"exit: {proc.returncode}\n"
            f"output:\n{proc.stdout}"
        )
    return proc.stdout


def _parse_results(output: str) -> Dict[str, BenchResult]:
    results: Dict[str, BenchResult] = {}
    for line in output.splitlines():
        m = RE_RESULT.match(line.strip())
        if not m:
            continue
        name = m.group("name")
        secs = float(m.group("secs"))
        iters = int(m.group("iters"))
        results[name] = BenchResult(secs=secs, iters=iters)
    if not results:
        raise RuntimeError(f"No RESULT lines parsed from output:\n{output}")
    return results


def _write_json(path: Path, results: Dict[str, BenchResult]) -> None:
    payload = {
        name: {"secs": r.secs, "iters": r.iters, "secs_per_iter": r.secs / r.iters if r.iters else None}
        for name, r in sorted(results.items())
    }
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def _fmt_float(x: float) -> str:
    if x == 0.0:
        return "0"
    if x < 0.001:
        return f"{x:.6f}"
    if x < 1.0:
        return f"{x:.6f}".rstrip("0").rstrip(".")
    return f"{x:.6f}".rstrip("0").rstrip(".")


def _compare(
    python: Dict[str, BenchResult],
    mojo: Dict[str, BenchResult],
) -> List[Tuple[str, float, float, float]]:
    # returns rows: (name, python_secs, mojo_secs, speedup)
    names = sorted(set(python.keys()) | set(mojo.keys()))
    rows: List[Tuple[str, float, float, float]] = []
    for name in names:
        py = python.get(name)
        mj = mojo.get(name)
        if py is None or mj is None:
            continue
        speedup = py.secs / mj.secs if mj.secs != 0 else float("inf")
        rows.append((name, py.secs, mj.secs, speedup))
    return rows


def _print_table(rows: List[Tuple[str, float, float, float]]) -> None:
    if not rows:
        print("No common benchmark names to compare.")
        return

    name_w = max(len(r[0]) for r in rows)
    print(f"{'name'.ljust(name_w)}  python_s  mojo_s  speedup(py/mojo)")
    for name, py_s, mj_s, sp in rows:
        print(
            f"{name.ljust(name_w)}  {_fmt_float(py_s):>8}  {_fmt_float(mj_s):>6}  {_fmt_float(sp):>15}"
        )


def main() -> int:
    root = _repo_root()
    out_dir = root / "bench" / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    mojo_out = _run_cmd(["uv", "run", "mojo", "run", "-I", "shapely_mojo", "bench/mojo_bench.mojo"], cwd=root)
    mojo_results = _parse_results(mojo_out)

    py_out = _run_cmd(["uv", "run", "python", "bench/python_bench.py"], cwd=root)
    py_results = _parse_results(py_out)

    mojo_json = out_dir / "mojo_results.json"
    py_json = out_dir / "python_results.json"

    _write_json(mojo_json, mojo_results)
    _write_json(py_json, py_results)

    rows = _compare(py_results, mojo_results)
    _print_table(rows)

    print(f"\nWrote: {mojo_json}")
    print(f"Wrote: {py_json}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
