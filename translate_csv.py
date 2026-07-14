#!/usr/bin/env python3
"""Translate multi-measurement UV-Vis CSV exports into per-measurement TXT files.

The new spectrophotometer CSV format stores each measurement as two adjacent
columns: wavelength and signal.  The first CSV row contains measurement names,
with blanks between names; the second row contains units such as
"Wavelength (nm),Abs" or "Wavelength (nm),%T".  This script writes one
`wavelength,absorption` TXT file per absorption measurement so the existing
ROOT analysis can keep using the `sample_N.txt` naming convention.

By default, the first two measurement pairs (the first four columns) are
skipped because they are transmission references from the instrument.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Iterable


def slugify(name: str) -> str:
    """Return a filesystem-safe measurement name while preserving useful tags."""
    slug = re.sub(r"\s+", "_", name.strip())
    slug = re.sub(r"[^A-Za-z0-9_.-]+", "_", slug)
    slug = re.sub(r"_+", "_", slug).strip("._-")
    return slug or "measurement"


def split_replicate_name(name: str) -> tuple[str, int | None]:
    """Split equipment names into a base name and optional replicate number.

    The instrument appears to name repeats as `Sample`, `Sample1`, `Sample2`,
    while the analysis expects `Sample_1.txt`, `Sample_2.txt`, `Sample_3.txt`.
    Names that already end in `_N` are left as explicit replicate numbers.
    """
    name = slugify(name)
    explicit = re.match(r"^(?P<base>.+)_(?P<rep>\d+)$", name)
    if explicit:
        return explicit.group("base"), int(explicit.group("rep"))

    implicit = re.match(r"^(?P<base>.*?)(?P<rep>\d+)$", name)
    if implicit and implicit.group("base"):
        return implicit.group("base"), int(implicit.group("rep")) + 1

    return name, None


def normalize_names(raw_names: Iterable[str]) -> list[str]:
    """Convert raw CSV measurement names to `base_N` stems."""
    parsed = [split_replicate_name(name) for name in raw_names]
    counts: dict[str, int] = {}
    normalized: list[str] = []

    for base, replicate in parsed:
        if replicate is None:
            counts[base] = counts.get(base, 0) + 1
            replicate = counts[base]
        else:
            counts[base] = max(counts.get(base, 0), replicate)
        normalized.append(f"{base}_{replicate}")

    return normalized


def convert_csv(input_csv: Path, output_dir: Path, skip_columns: int = 4) -> list[Path]:
    with input_csv.open(newline="") as handle:
        rows = list(csv.reader(handle))

    if len(rows) < 3:
        raise ValueError(f"{input_csv} does not contain the expected two header rows plus data")
    if skip_columns % 2:
        raise ValueError("skip_columns must be even because each measurement uses two columns")

    names_row, units_row, data_rows = rows[0], rows[1], rows[2:]
    output_dir.mkdir(parents=True, exist_ok=True)

    measurement_pairs: list[tuple[int, str]] = []
    for col in range(skip_columns, len(units_row) - 1, 2):
        y_unit = units_row[col + 1].strip().lower()
        if y_unit != "abs":
            continue
        raw_name = names_row[col].strip()
        if not raw_name:
            continue
        measurement_pairs.append((col, raw_name))

    output_stems = normalize_names(name for _, name in measurement_pairs)
    written: list[Path] = []

    for (col, raw_name), stem in zip(measurement_pairs, output_stems):
        out_path = output_dir / f"{stem}.txt"
        points: list[tuple[str, str]] = []
        for row in data_rows:
            if col + 1 >= len(row):
                continue
            wavelength = row[col].strip()
            absorbance = row[col + 1].strip()
            if wavelength and absorbance:
                points.append((wavelength, absorbance))

        with out_path.open("w", newline="") as handle:
            handle.write(f"{raw_name}\n")
            handle.write("Wavelength (nm),Abs\n")
            for wavelength, absorbance in points:
                handle.write(f"{wavelength},{absorbance}\n")
        written.append(out_path)

    return written


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Split a multi-measurement UV-Vis CSV export into analysis-ready TXT files."
    )
    parser.add_argument("input_csv", type=Path, help="CSV export from the new UV-Vis equipment")
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path("."),
        help="directory where TXT files will be written (default: current directory)",
    )
    parser.add_argument(
        "--skip-columns",
        type=int,
        default=4,
        help="number of leading columns to skip (default: 4, i.e. two transmission measurements)",
    )
    args = parser.parse_args()

    written = convert_csv(args.input_csv, args.output_dir, args.skip_columns)
    print(f"Wrote {len(written)} txt files to {args.output_dir}")
    for path in written:
        print(path)


if __name__ == "__main__":
    main()
