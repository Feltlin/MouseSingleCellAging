"""Run the per-tissue notebooks (02-05c) top to bottom for one or more tissues.

    python src/run_tissue.py Lung
    python src/run_tissue.py Kidney Spleen --from 05a

Each notebook reads TISSUE from the environment. Executed copies, with outputs,
go to notebook/run/{tissue}/ (ignored by git); the notebooks in notebook/ are not touched.
A tissue stops at its first failing notebook, and the error is saved in that copy.
"""
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient
from nbclient.exceptions import CellExecutionError

ROOT = Path(__file__).resolve().parent.parent
NOTEBOOKS = ["02_explore", "03_preprocess", "04_composition", "04b_lung_outlier",
             "05a_classifier", "05b_genes", "05c_clock"]


def run(tissue: str, start: str | None = None, timeout: int = 3600) -> bool:
    out_dir = ROOT / "notebook" / "run" / tissue
    out_dir.mkdir(parents=True, exist_ok=True)
    names = [n for n in NOTEBOOKS if n != "04b_lung_outlier" or tissue == "Lung"]   # 30-M-2 detour: Lung only
    names = names[names.index(start):] if start else names
    os.environ["TISSUE"] = tissue
    for name in names:
        t0 = time.time()
        nb = nbformat.read(ROOT / "notebook" / f"{name}.ipynb", as_version=4)
        client = NotebookClient(nb, timeout=timeout, kernel_name="python3",
                                resources={"metadata": {"path": str(ROOT / "notebook")}})
        try:
            client.execute()
            ok = True
        except CellExecutionError as err:
            print(f"[{tissue}] {name} FAILED:\n{str(err)[-1500:]}", flush=True)
            ok = False
        nbformat.write(nb, out_dir / f"{name}.ipynb")
        print(f"[{tissue}] {name} {'ok' if ok else 'failed'} ({time.time() - t0:.0f}s)", flush=True)
        if not ok:
            return False
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("tissues", nargs="+")
    parser.add_argument("--from", dest="start", choices=NOTEBOOKS, help="first notebook to run")
    args = parser.parse_args()
    results = {t: run(t, args.start) for t in args.tissues}
    print("\nsummary:", {t: "ok" if r else "failed" for t, r in results.items()})
    sys.exit(0 if all(results.values()) else 1)
