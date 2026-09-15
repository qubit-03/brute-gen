#!/usr/bin/env python3
"""
search_cas.py — Find which PDF/XLSX files in a folder contain a given CAS number.

Fast by default: uses poppler's `pdftotext` (C++, very fast) for PDFs and
openpyxl read-only mode for XLSX, and searches files in parallel across CPU
cores. Falls back to slower table-aware extraction (pdfplumber) only for PDFs
where the fast text pass finds nothing — handles PDFs where columns are
laid out in a way plain text extraction garbles.

Usage:
    python search_cas.py <folder> <cas_number> [--exact] [--workers N] [--no-fallback]

Examples:
    python search_cas.py ./chemicals 50-00-0
    python search_cas.py ./chemicals 50-00-0 --exact
    python search_cas.py ./chemicals 50-00-0 --workers 8
    python search_cas.py ./chemicals 50-00-0 --no-fallback   # skip pdfplumber fallback entirely (fastest)
"""

import sys
import re
import shutil
import argparse
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    import openpyxl
except ImportError:
    openpyxl = None

HAVE_PDFTOTEXT = shutil.which("pdftotext") is not None


def normalize(s: str) -> str:
    return re.sub(r"\s+", " ", str(s)).strip()


def cas_matches(text: str, cas: str, exact: bool) -> bool:
    if exact:
        pattern = r"(?<![0-9-])" + re.escape(cas) + r"(?![0-9-])"
        return re.search(pattern, text) is not None
    return cas in text


def search_pdf_fast(path: Path, cas: str, exact: bool):
    """Fast path: pdftotext -layout, piped straight to memory, line-by-line grep."""
    try:
        result = subprocess.run(
            ["pdftotext", "-layout", str(path), "-"],
            capture_output=True, text=True, timeout=60
        )
        text = result.stdout
    except Exception:
        return []
    matches = []
    for line in text.split("\n"):
        if cas_matches(line, cas, exact):
            matches.append((None, normalize(line)))  # no reliable page number from this mode
    return matches


def search_pdf_fallback(path: Path, cas: str, exact: bool):
    """Slow but thorough: pdfplumber, page text + table extraction."""
    import pdfplumber
    matches = []
    try:
        with pdfplumber.open(path) as pdf:
            for i, page in enumerate(pdf.pages, start=1):
                text = page.extract_text() or ""
                for line in text.split("\n"):
                    if cas_matches(line, cas, exact):
                        matches.append((i, normalize(line)))
                for table in page.extract_tables() or []:
                    for row in table:
                        row_text = " | ".join(c for c in row if c)
                        if cas_matches(row_text, cas, exact):
                            matches.append((i, normalize(row_text)))
    except Exception as e:
        return [(None, f"[!] Error reading with pdfplumber: {e}")]
    return matches


def search_xlsx(path: Path, cas: str, exact: bool):
    matches = []
    try:
        wb = openpyxl.load_workbook(path, data_only=True, read_only=True)
        for sheet in wb.worksheets:
            for row_idx, row in enumerate(sheet.iter_rows(values_only=True), start=1):
                row_text = " | ".join(str(c) for c in row if c is not None)
                if cas_matches(row_text, cas, exact):
                    matches.append((sheet.title, row_idx, normalize(row_text)))
        wb.close()
    except Exception as e:
        return [(None, None, f"[!] Error reading {path.name}: {e}")]
    return matches


def process_file(path_str, cas, exact, use_fallback):
    """Runs in a worker process. Returns (path_str, kind, matches)."""
    path = Path(path_str)
    suffix = path.suffix.lower()

    if suffix == ".pdf":
        if HAVE_PDFTOTEXT:
            matches = search_pdf_fast(path, cas, exact)
            if not matches and use_fallback:
                matches = search_pdf_fallback(path, cas, exact)
        else:
            matches = search_pdf_fallback(path, cas, exact)
        return (path_str, "pdf", matches)

    elif suffix == ".xlsx":
        matches = search_xlsx(path, cas, exact)
        return (path_str, "xlsx", matches)

    return (path_str, None, [])


def main():
    parser = argparse.ArgumentParser(description="Search PDF/XLSX files for a CAS number (fast, parallel).")
    parser.add_argument("folder", help="Folder containing PDF/XLSX files")
    parser.add_argument("cas", help="CAS number to search for, e.g. 50-00-0")
    parser.add_argument("--exact", action="store_true",
                         help="Match CAS as a standalone token, not just a substring")
    parser.add_argument("--workers", type=int, default=None,
                         help="Number of parallel worker processes (default: CPU count)")
    parser.add_argument("--no-fallback", action="store_true",
                         help="Skip slow pdfplumber fallback for PDFs with no fast-text match")
    args = parser.parse_args()

    folder = Path(args.folder)
    if not folder.is_dir():
        print(f"Error: '{folder}' is not a valid folder.")
        sys.exit(1)

    if not HAVE_PDFTOTEXT:
        print("Note: 'pdftotext' (poppler-utils) not found — falling back to slower pdfplumber for all PDFs.")
        print("Install it for a big speedup: sudo apt install poppler-utils\n")

    cas = args.cas.strip()
    files = sorted(list(folder.glob("*.pdf")) + list(folder.glob("*.xlsx")))

    if not files:
        print(f"No .pdf or .xlsx files found in {folder}")
        sys.exit(0)

    print(f"Searching {len(files)} file(s) in '{folder}' for CAS '{cas}'"
          f"{' (exact match)' if args.exact else ''}...\n")

    found_any = False
    use_fallback = not args.no_fallback

    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {
            executor.submit(process_file, str(p), cas, args.exact, use_fallback): p
            for p in files
        }
        for future in as_completed(futures):
            path = futures[future]
            try:
                path_str, kind, matches = future.result()
            except Exception as e:
                print(f"  [!] Error processing {path.name}: {e}")
                continue

            if not matches:
                continue

            found_any = True
            print(f"✔ FOUND in: {path.name}")
            for m in matches:
                if kind == "pdf":
                    page_num, line = m
                    loc = f"Page {page_num}" if page_num else "Match"
                    print(f"    {loc}: {line}")
                elif kind == "xlsx":
                    sheet_name, row_num, row = m
                    if sheet_name is None:
                        print(f"    {row}")
                    else:
                        print(f"    Sheet '{sheet_name}' Row {row_num}: {row}")
            print()

    if not found_any:
        print(f"No matches for CAS '{cas}' found in any file.")


if __name__ == "__main__":
    main()