#!/usr/bin/env python3
"""
generate_gpd_manifest.py

Walk through
    /gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/G-P-D/
and build a CSV manifest with three columns:

    rsID,PHID,sex

Each .tsv file is assumed to live exactly one directory below G-P-D
and to follow the naming convention

    <PHID>_<sex>.tsv

Examples
--------
Path                                      ->  CSV row created
-----------------------------------------    -------------------------------
.../rs8491/PHD04005_female.tsv           ->  rs8491,PHD04005,female
.../rs63319/PHD10006_all.tsv             ->  rs63319,PHD10006,all
"""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterator, NamedTuple


ROOT = Path(
    "/gpfs/chencao/ysbioinfor/project/proteohubProject/web_out/G-P-D"
).resolve()

#: captures “PHID” and “sex” from filenames like  PHD04005_female.tsv
_FILENAME_RE = re.compile(r"^(?P<phid>[^_]+)_(?P<sex>[^.]+)\.tsv$", re.IGNORECASE)


class Record(NamedTuple):
    rs_id: str
    phid: str
    sex: str


def iter_records(root: Path) -> Iterator[Record]:
    """
    Yield one `Record` for every *.tsv file found beneath *root*.
    """
    for tsv_path in root.rglob("*.tsv"):
        # Skip directories such as root/xxx.tsv that are not in an rsID folder
        if tsv_path.parent == root:
            continue

        rs_id = tsv_path.parent.name
        if not rs_id.startswith("rs"):
            # Not an rsID directory; ignore
            continue

        m = _FILENAME_RE.match(tsv_path.name)
        if m is None:  # malformed file name → skip
            continue

        yield Record(rs_id=rs_id, phid=m.group("phid"), sex=m.group("sex"))


def write_csv(records: Iterator[Record], out_file: Path, overwrite: bool = False) -> None:
    """
    Write *records* to *out_file* in CSV format with header.
    """
    if out_file.exists() and not overwrite:
        raise FileExistsError(f"{out_file} already exists (use --force to overwrite)")

    out_file.parent.mkdir(parents=True, exist_ok=True)

    with out_file.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(("rsID", "PHID", "sex"))
        for rec in sorted(records):  # ensure deterministic order
            writer.writerow(rec)


def parse_cli() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate rsID/PHID/sex manifest for G-P-D results"
    )
    p.add_argument(
        "-r",
        "--root",
        type=Path,
        default=ROOT,
        help=f"Root directory to crawl (default: {ROOT})",
    )
    p.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("gpd_manifest.csv"),
        help="Destination CSV file (default: ./gpd_manifest.csv)",
    )
    p.add_argument(
        "-f",
        "--force",
        action="store_true",
        help="Overwrite output file if it already exists",
    )
    return p.parse_args()


def main() -> None:
    ns = parse_cli()
    recs = list(iter_records(ns.root))
    write_csv(recs, ns.output, overwrite=ns.force)
    print(f"Created {ns.output} with {len(recs):,} rows.")


if __name__ == "__main__":
    main()