#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import numpy as np

def read_table(path, id_col=0):
    """
    Read a delimited file (CSV/TSV; gz OK) and set the ID column as index.
    Auto-detects delimiter. Leaves numeric strings as numbers when possible.
    """
    try:
        df = pd.read_csv(
            path,
            sep=None,                  # auto-detect delimiter
            engine="python",           # needed for sep=None
            compression="infer",       # supports .gz
        )
    except Exception as e:
        print(f"[ERROR] Failed to read '{path}': {e}", file=sys.stderr)
        sys.exit(2)

    # Normalize column types (strip spaces)
    df.columns = [str(c).strip() for c in df.columns]

    # Set ID column as index (by position or name)
    if isinstance(id_col, int):
        id_name = df.columns[id_col]
    else:
        if id_col not in df.columns:
            print(f"[ERROR] ID column '{id_col}' not found in {path}", file=sys.stderr)
            sys.exit(2)
        id_name = id_col

    # Ensure no duplicate IDs
    if df[id_name].duplicated().any():
        dups = df[id_name][df[id_name].duplicated()].unique()
        print(f"[ERROR] Duplicate IDs in {path}: {list(dups)[:10]}...", file=sys.stderr)
        sys.exit(2)

    df = df.set_index(id_name)

    # Convert remaining columns to numeric where possible
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")

    return df


def compare_tables(df1, df2, rtol=1e-3, atol=1e-3, verbose=False):
    """
    Compare two dataframes with the same schema (index = IDs; columns = samples).
    Returns a dict with comparison results and differences beyond tolerance.
    """
    # Check column sets (ignore order)
    cols1, cols2 = set(df1.columns), set(df2.columns)
    missing_in_2 = sorted(cols1 - cols2)
    missing_in_1 = sorted(cols2 - cols1)

    # Check row (ID) sets (ignore order)
    ids1, ids2 = set(df1.index), set(df2.index)
    ids_missing_in_2 = sorted(ids1 - ids2)
    ids_missing_in_1 = sorted(ids2 - ids1)

    results = {
        "columns_equal": (cols1 == cols2),
        "rows_equal": (ids1 == ids2),
        "missing_columns_in_file2": missing_in_2,
        "missing_columns_in_file1": missing_in_1,
        "missing_ids_in_file2": ids_missing_in_2,
        "missing_ids_in_file1": ids_missing_in_1,
        "mismatch_count": 0,
        "mismatches": None,   # DataFrame of diffs where not close
    }

    if not results["columns_equal"] or not results["rows_equal"]:
        return results  # schema mismatch; stop here

    # Align both by sorted IDs and sorted columns for deterministic comparison
    common_cols = sorted(cols1)  # same as cols2 now
    common_ids = sorted(ids1)    # same as ids2 now
    a = df1.loc[common_ids, common_cols]
    b = df2.loc[common_ids, common_cols]

    # Use numpy.isclose with tolerances
    close_mask = np.isclose(a.values, b.values, rtol=rtol, atol=atol, equal_nan=True)
    # Build DataFrame of differences where not close
    not_close = ~close_mask
    if not_close.any():
        diff = (a - b)
        diff_masked = pd.DataFrame(diff.values, index=a.index, columns=a.columns)
        # Keep only rows/cols where there is a mismatch
        # (filter down to cells failing tolerance)
        where_mismatch = []
        for i, rid in enumerate(diff_masked.index):
            for j, col in enumerate(diff_masked.columns):
                if not_close[i, j]:
                    where_mismatch.append((rid, col, a.iloc[i, j], b.iloc[i, j], diff_masked.iloc[i, j]))

        # Summarize into a DataFrame
        mismatches_df = pd.DataFrame(where_mismatch, columns=["ID", "Column", "Value_file1", "Value_file2", "Diff(file1-file2)"])
        results["mismatch_count"] = len(mismatches_df)
        results["mismatches"] = mismatches_df

    return results


def main():
    p = argparse.ArgumentParser(
        description="Compare two wide tables (same IDs/columns, any order) with numeric tolerances."
    )
    p.add_argument("file1", help="First table (CSV/TSV; .gz OK)")
    p.add_argument("file2", help="Second table (CSV/TSV; .gz OK)")
    p.add_argument("--id-col", default=0,
                   help="ID column (name or 0-based index). Default: 0")
    p.add_argument("--rtol", type=float, default=1e-3,
                   help="Relative tolerance for numeric comparison (np.isclose). Default: 1e-3")
    p.add_argument("--atol", type=float, default=1e-3,
                   help="Absolute tolerance for numeric comparison (np.isclose). Default: 1e-3")
    p.add_argument("--verbose", action="store_true", help="Print extra details")
    p.add_argument("--show-max", type=int, default=50, help="Max mismatches to print. Default: 50")

    args = p.parse_args()

    # Coerce id_col to int if numeric-like
    try:
        id_col = int(args.id_col)
    except ValueError:
        id_col = args.id_col

    df1 = read_table(args.file1, id_col=id_col)
    df2 = read_table(args.file2, id_col=id_col)

    res = compare_tables(df1, df2, rtol=args.rtol, atol=args.atol, verbose=args.verbose)

    # Report schema checks
    if not res["columns_equal"] or not res["rows_equal"]:
        print("SCHEMA MISMATCH")
        if not res["columns_equal"]:
            if res["missing_columns_in_file2"]:
                print(f"  Columns present only in file1: {res['missing_columns_in_file2']}")
            if res["missing_columns_in_file1"]:
                print(f"  Columns present only in file2: {res['missing_columns_in_file1']}")
        if not res["rows_equal"]:
            if res["missing_ids_in_file2"]:
                print(f"  IDs present only in file1: {res['missing_ids_in_file2'][:20]}{' ...' if len(res['missing_ids_in_file2'])>20 else ''}")
            if res["missing_ids_in_file1"]:
                print(f"  IDs present only in file2: {res['missing_ids_in_file1'][:20]}{' ...' if len(res['missing_ids_in_file1'])>20 else ''}")
        sys.exit(1)

    # Report cell-wise comparison
    if res["mismatch_count"] == 0:
        print(f"OK: All {df1.shape[0]} rows and {df1.shape[1]} columns match within rtol={args.rtol}, atol={args.atol}.")
        sys.exit(0)
    else:
        print(f"MISMATCH: {res['mismatch_count']} cell(s) differ beyond rtol={args.rtol}, atol={args.atol}.")
        print(res["mismatches"].head(args.show_max).to_string(index=False))
        if res["mismatch_count"] > args.show_max:
            print(f"... {res['mismatch_count'] - args.show_max} more mismatches not shown.")
        sys.exit(1)


if __name__ == "__main__":
    main()
