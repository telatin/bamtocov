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


def compare_tables(df1, df2, pct_diff=None, abs_diff=None, verbose=False):
    """
    Compare two dataframes with the same schema (index = IDs; columns = samples).
    Returns a dict with comparison results and differences beyond tolerance.

    Tolerance modes:
      - pct_diff: Percentage tolerance relative to the larger absolute value
      - abs_diff: Absolute difference tolerance
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

    # Calculate absolute differences
    abs_diffs = np.abs(a.values - b.values)

    # Determine which values are within tolerance
    if pct_diff is not None:
        # Percentage tolerance: calculate as % of the larger absolute value
        max_vals = np.maximum(np.abs(a.values), np.abs(b.values))
        tolerance = (pct_diff / 100.0) * max_vals
        # Handle NaN: if both are NaN, consider them equal
        both_nan = np.isnan(a.values) & np.isnan(b.values)
        close_mask = (abs_diffs <= tolerance) | both_nan
    elif abs_diff is not None:
        # Absolute tolerance
        both_nan = np.isnan(a.values) & np.isnan(b.values)
        close_mask = (abs_diffs <= abs_diff) | both_nan
    else:
        # Default: exact match (or both NaN)
        both_nan = np.isnan(a.values) & np.isnan(b.values)
        close_mask = (abs_diffs == 0) | both_nan

    # Build DataFrame of differences where not within tolerance
    not_close = ~close_mask
    if not_close.any():
        # Keep only cells where there is a mismatch
        where_mismatch = []
        for i, rid in enumerate(a.index):
            for j, col in enumerate(a.columns):
                if not_close[i, j]:
                    val1 = a.iloc[i, j]
                    val2 = b.iloc[i, j]
                    abs_diff_val = abs_diffs[i, j]
                    where_mismatch.append((rid, col, val1, val2, abs_diff_val))

        # Summarize into a DataFrame
        mismatches_df = pd.DataFrame(where_mismatch, columns=["ID", "Column", "Value_file1", "Value_file2", "AbsDiff"])
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

    # Tolerance options (mutually exclusive)
    tol_group = p.add_mutually_exclusive_group()
    tol_group.add_argument("--pct-diff", type=float, metavar="PCT",
                          help="Percentage tolerance (e.g., 10 for 10%% of the larger absolute value)")
    tol_group.add_argument("--abs-diff", type=float, metavar="VALUE",
                          help="Absolute difference tolerance (e.g., 0.1 to allow ±0.1)")

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

    res = compare_tables(df1, df2, pct_diff=args.pct_diff, abs_diff=args.abs_diff, verbose=args.verbose)

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
        if args.pct_diff is not None:
            tol_str = f"pct-diff={args.pct_diff}%"
        elif args.abs_diff is not None:
            tol_str = f"abs-diff={args.abs_diff}"
        else:
            tol_str = "exact match"
        print(f"OK: All {df1.shape[0]} rows and {df1.shape[1]} columns match within {tol_str}.")
        sys.exit(0)
    else:
        if args.pct_diff is not None:
            tol_str = f"pct-diff={args.pct_diff}%"
        elif args.abs_diff is not None:
            tol_str = f"abs-diff={args.abs_diff}"
        else:
            tol_str = "exact match"
        print(f"MISMATCH: {res['mismatch_count']} cell(s) differ beyond {tol_str}.")
        print(res["mismatches"].head(args.show_max).to_string(index=False))
        if res["mismatch_count"] > args.show_max:
            print(f"... {res['mismatch_count'] - args.show_max} more mismatches not shown.")
        sys.exit(1)


if __name__ == "__main__":
    main()
