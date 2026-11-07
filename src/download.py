"""Download Tabula Muris Senis files from figshare, with an md5 check.

    python src/download.py --list                      # show files and sizes, download nothing
    python src/download.py droplet Kidney Spleen       # selected droplet tissues
    python src/download.py droplet all                 # every droplet tissue
    python src/download.py external                    # eLife gene tables + Kimmel 2019 data

Tissue files go to data/{tissue}/ (the name the notebooks expect, see ag.paths);
external files go to data/external/. Files already present with the right md5 are skipped.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
TISSUE_ARTICLE = 12654728     # "Tabula Muris Senis Data Objects": one h5ad per tissue and method
GENE_ARTICLE = 12827615       # "TMS gene data (final)": eLife 2021 results, Kimmel 2019 data
EXTERNAL = ["global_aging_genes.tsv", "aging_score_model.tsv", "DGE_result_release_sig.zip",
            "Kimmel_GR_2019_data.zip", "190304_maca_bulk.h5ad"]
DROPLET_PREFIX = "tabula-muris-senis-droplet-processed-official-annotations-"
USER_AGENT = "single-cell-aging-download/1.0"


def file_list(article: int) -> list[dict]:
    req = urllib.request.Request(f"https://api.figshare.com/v2/articles/{article}/files?page_size=100",
                                 headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as r:
        return json.load(r)


def droplet_files() -> dict[str, dict]:
    """{tissue: figshare file record} for every droplet tissue."""
    return {f["name"][len(DROPLET_PREFIX):-len(".h5ad")]: f
            for f in file_list(TISSUE_ARTICLE) if f["name"].startswith(DROPLET_PREFIX)}


def md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch(record: dict, dest: Path) -> None:
    if dest.exists() and md5(dest) == record["computed_md5"]:
        print(f"have  {dest.relative_to(ROOT)}")
        return
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    print(f"get   {dest.relative_to(ROOT)} ({record['size'] / 1e6:,.0f} MB)", flush=True)
    req = urllib.request.Request(record["download_url"], headers={"User-Agent": USER_AGENT})
    with urllib.request.urlopen(req) as r, open(tmp, "wb") as out:
        while chunk := r.read(1 << 20):
            out.write(chunk)
    if md5(tmp) != record["computed_md5"]:
        tmp.unlink()
        raise RuntimeError(f"md5 mismatch for {dest.name}; partial file removed")
    tmp.replace(dest)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("what", nargs="?", choices=["droplet", "external"])
    p.add_argument("tissues", nargs="*", help='tissue names, or "all"')
    p.add_argument("--list", action="store_true", help="list files and sizes only")
    args = p.parse_args()

    drop = droplet_files()
    if args.list or args.what is None:
        print("droplet tissues:")
        for t, f in sorted(drop.items(), key=lambda kv: kv[1]["size"]):
            have = (ROOT / "data" / t / f["name"]).exists()
            print(f"  {t:<18} {f['size'] / 1e6:>8,.0f} MB{'   (have)' if have else ''}")
        ext = {f["name"]: f for f in file_list(GENE_ARTICLE)}
        print("external:")
        for n in EXTERNAL:
            print(f"  {n:<32} {ext[n]['size'] / 1e6:>8,.1f} MB")
        return

    if args.what == "droplet":
        names = sorted(drop) if args.tissues == ["all"] else args.tissues
        unknown = [t for t in names if t not in drop]
        if unknown:
            sys.exit(f"unknown tissue(s): {unknown}; choose from {sorted(drop)}")
        for t in names:
            fetch(drop[t], ROOT / "data" / t / drop[t]["name"])
    else:
        ext = {f["name"]: f for f in file_list(GENE_ARTICLE)}
        for n in EXTERNAL:
            fetch(ext[n], ROOT / "data" / "external" / n)


if __name__ == "__main__":
    main()
