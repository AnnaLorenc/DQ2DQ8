#!/usr/bin/env python3
"""Script: aa_kidera_test.py

Summary
-------
Performs multivariate (Hotelling's T²) and univariate (per-Kidera-factor Welch
t-test) comparisons of Kidera-factor profiles across genotype groups, for each
combination of aggregation columns (e.g. IMGT_position × cells, or
IMGT_position × length × cells).

The input is the output of aa_kidera_aggregate.py: one row per
(agg_cols × patient × sample × genotype), with KF1..KF10 columns.

Two comparisons are performed:
  - oo  :  hoDQ2  vs  hoDQ8
  - eo  :  heDQ2DQ8  vs  (hoDQ2 + hoDQ8)

Both comparisons are run within each cell type (implicit if "cells" is among
the aggregation columns).

Output files
------------
1. <output_dir>/kidera_hotelling.tsv
   One row per aggregation-column combination.
   Columns (side-by-side _oo and _eo):
     T2, F, df1, df2, p_value, q_value (BH-FDR), eta_sq, n_a, n_b

2. <output_dir>/kidera_per_factor.tsv
   One row per (aggregation-column combination × Kidera factor).
   Columns (side-by-side _oo and _eo):
     mean_a, mean_b, mean_diff, cohens_d, t_stat, df, p_value, q_value (BH-FDR)
   Plus: p_hotelling_oo, q_hotelling_oo, p_hotelling_eo, q_hotelling_eo

Usage
-----
    python aa_kidera_test.py \\
        --input  sandbox/kidera_aggregated.tsv \\
        --output_dir sandbox/kidera_results \\
        [--agg_columns IMGT_position cells] \\
        [--kf_columns KF1 KF2 KF3 KF4 KF5 KF6 KF7 KF8 KF9 KF10]
"""

import argparse
import math
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from scipy.stats import f as f_dist

try:
    from statsmodels.stats.multitest import multipletests
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False
    print("Warning: statsmodels not found; q-values will be NaN.")


DEFAULT_AGG_COLUMNS = ["IMGT_position", "cells"]
DEFAULT_KF_COLUMNS  = [f"KF{i}" for i in range(1, 11)]


# ---------------------------------------------------------------------------
# Hotelling's T² for two independent groups
# ---------------------------------------------------------------------------

def hotelling_t2(
    X: np.ndarray,
    Y: np.ndarray,
) -> Optional[dict]:
    """Compute Hotelling's T² for two independent multivariate groups.

    Parameters
    ----------
    X : (n_a, p) array for group A
    Y : (n_b, p) array for group B

    Returns
    -------
    dict with keys: T2, F, df1, df2, p_value, eta_sq, n_a, n_b
    or None if the test cannot be computed (insufficient data or singular matrix).

    Notes
    -----
    T² = (n_a·n_b)/(n_a+n_b) · (x̄_a − x̄_b)ᵀ S_pooled⁻¹ (x̄_a − x̄_b)
    F  = (n_a+n_b−p−1) / [p·(n_a+n_b−2)] · T²
    F  ~ F(p, n_a+n_b−p−1)
    η² = T² / (T² + n_a+n_b−2)   [partial eta-squared analog]
    """
    n_a, p = X.shape
    n_b     = Y.shape[0]

    if n_a < 2 or n_b < 2:
        return None
    if n_a + n_b - 2 < p:
        # Degrees of freedom would be non-positive; pooled cov is rank-deficient
        return None

    mean_diff = X.mean(axis=0) - Y.mean(axis=0)

    # Pooled covariance
    S_a = np.cov(X, rowvar=False, ddof=1) if n_a > 1 else np.zeros((p, p))
    S_b = np.cov(Y, rowvar=False, ddof=1) if n_b > 1 else np.zeros((p, p))
    S_pooled = ((n_a - 1) * S_a + (n_b - 1) * S_b) / (n_a + n_b - 2)

    try:
        S_inv = np.linalg.inv(S_pooled)
    except np.linalg.LinAlgError:
        return None

    # Check condition number – if very large, matrix is near-singular
    if np.linalg.cond(S_pooled) > 1e12:
        return None

    T2 = (n_a * n_b / (n_a + n_b)) * float(mean_diff @ S_inv @ mean_diff)

    df1 = p
    df2 = n_a + n_b - p - 1
    if df2 <= 0:
        return None

    F = (n_a + n_b - p - 1) / (p * (n_a + n_b - 2)) * T2
    p_value = float(1.0 - f_dist.cdf(F, df1, df2))
    eta_sq  = T2 / (T2 + n_a + n_b - 2)

    return {
        "T2": T2, "F": F, "df1": df1, "df2": df2,
        "p_value": p_value, "eta_sq": eta_sq,
        "n_a": n_a, "n_b": n_b,
    }


# ---------------------------------------------------------------------------
# Per-factor Welch t-test with Cohen's d
# ---------------------------------------------------------------------------

def welch_t_cohens_d(
    a: np.ndarray,
    b: np.ndarray,
) -> Optional[dict]:
    """Welch two-sample t-test plus Cohen's d for two 1-D arrays.

    Cohen's d uses the pooled SD (Glass's variant would use SD of one group;
    we use pooled to be consistent with the paired two-group framing).

    Returns dict or None if either group has fewer than 2 observations.
    """
    a = a[~np.isnan(a)]
    b = b[~np.isnan(b)]
    if len(a) < 2 or len(b) < 2:
        return None

    mean_a  = float(a.mean())
    mean_b  = float(b.mean())
    diff    = mean_a - mean_b

    var_a   = float(a.var(ddof=1))
    var_b   = float(b.var(ddof=1))
    n_a, n_b = len(a), len(b)

    # Pooled SD for Cohen's d
    pooled_var = ((n_a - 1) * var_a + (n_b - 1) * var_b) / (n_a + n_b - 2)
    pooled_sd  = math.sqrt(pooled_var) if pooled_var > 0 else float("nan")
    cohens_d   = diff / pooled_sd if pooled_sd > 0 else float("nan")

    t_stat, p_value = scipy_stats.ttest_ind(a, b, equal_var=False)

    se_a = var_a / n_a
    se_b = var_b / n_b
    se_sum = se_a + se_b
    denom_df = (se_a ** 2 / (n_a - 1)) + (se_b ** 2 / (n_b - 1))
    if se_sum == 0 or denom_df == 0:
        # Both groups have zero variance – data are constant; test is undefined
        return None
    df_welch = float(se_sum ** 2 / denom_df)

    return {
        "mean_a": mean_a, "mean_b": mean_b, "mean_diff": diff,
        "cohens_d": cohens_d,
        "t_stat": float(t_stat), "df": df_welch,
        "p_value": float(p_value),
    }


# ---------------------------------------------------------------------------
# FDR correction helper
# ---------------------------------------------------------------------------

def fdr_correct(p_values: pd.Series) -> pd.Series:
    """Apply Benjamini–Hochberg FDR correction to a Series of p-values.

    NaN p-values are preserved as NaN in the output.
    """
    if not HAS_STATSMODELS:
        return pd.Series([float("nan")] * len(p_values), index=p_values.index)

    mask = p_values.notna()
    q = pd.Series([float("nan")] * len(p_values), index=p_values.index)
    if mask.sum() == 0:
        return q

    _, q_vals, _, _ = multipletests(p_values[mask].values, method="fdr_bh")
    q[mask] = q_vals
    return q


# ---------------------------------------------------------------------------
# Main testing logic
# ---------------------------------------------------------------------------

def run_hotelling(
    df: pd.DataFrame,
    agg_cols: list[str],
    kf_cols: list[str],
    genotypes_a: list[str],
    genotypes_b: list[str],
) -> pd.DataFrame:
    """Run Hotelling's T² per agg_cols group for one comparison.

    Returns a DataFrame with one row per group and Hotelling stats columns.
    """
    rows = []
    for key, sub in df.groupby(agg_cols):
        key_dict = dict(zip(agg_cols, key if isinstance(key, tuple) else (key,)))

        X = sub[sub["genotype"].isin(genotypes_a)][kf_cols].dropna().values
        Y = sub[sub["genotype"].isin(genotypes_b)][kf_cols].dropna().values

        result = hotelling_t2(X, Y)
        row = {**key_dict}
        if result:
            row.update(result)
        else:
            row.update({k: float("nan") for k in
                        ["T2", "F", "df1", "df2", "p_value", "eta_sq", "n_a", "n_b"]})
        rows.append(row)

    out = pd.DataFrame(rows)
    out["q_value"] = fdr_correct(out["p_value"])
    return out


def run_per_factor(
    df: pd.DataFrame,
    agg_cols: list[str],
    kf_cols: list[str],
    genotypes_a: list[str],
    genotypes_b: list[str],
) -> pd.DataFrame:
    """Run per-KF Welch t-tests per agg_cols group for one comparison.

    Returns a DataFrame with one row per (group × KF).
    """
    rows = []
    for key, sub in df.groupby(agg_cols):
        key_dict = dict(zip(agg_cols, key if isinstance(key, tuple) else (key,)))

        sub_a = sub[sub["genotype"].isin(genotypes_a)]
        sub_b = sub[sub["genotype"].isin(genotypes_b)]

        for kf in kf_cols:
            a = sub_a[kf].dropna().values
            b = sub_b[kf].dropna().values
            result = welch_t_cohens_d(a, b)
            row = {**key_dict, "KF": kf}
            if result:
                row.update(result)
            else:
                row.update({k: float("nan") for k in
                            ["mean_a", "mean_b", "mean_diff", "cohens_d",
                             "t_stat", "df", "p_value"]})
            rows.append(row)

    out = pd.DataFrame(rows)
    out["q_value"] = fdr_correct(out["p_value"])
    return out


def run_tests(
    df: pd.DataFrame,
    agg_cols: list[str],
    kf_cols: list[str],
    output_dir: Path,
) -> None:
    """Orchestrate oo and eo comparisons, write output files."""

    comparisons = {
        "oo": (["hoDQ2"],     ["hoDQ8"],           "hoDQ2", "hoDQ8"),
        "eo": (["heDQ2DQ8"],  ["hoDQ2", "hoDQ8"],  "heDQ2DQ8", "hom"),
    }

    # ---- Hotelling --------------------------------------------------------
    hotelling_frames = {}
    for tag, (ga, gb, _, _) in comparisons.items():
        hotelling_frames[tag] = run_hotelling(df, agg_cols, kf_cols, ga, gb)

    # Merge oo and eo side by side on agg_cols
    h_oo = hotelling_frames["oo"].rename(
        columns={c: f"{c}_oo" for c in hotelling_frames["oo"].columns if c not in agg_cols}
    )
    h_eo = hotelling_frames["eo"].rename(
        columns={c: f"{c}_eo" for c in hotelling_frames["eo"].columns if c not in agg_cols}
    )
    df_hotelling = pd.merge(h_oo, h_eo, on=agg_cols, how="outer").sort_values(agg_cols)

    out1 = output_dir / "kidera_hotelling.tsv"
    df_hotelling.to_csv(out1, sep="\t", index=False)
    print(f"Written Hotelling results ({len(df_hotelling)} rows) → {out1}")

    # ---- Per-factor -------------------------------------------------------
    pf_frames = {}
    for tag, (ga, gb, _, _) in comparisons.items():
        pf_frames[tag] = run_per_factor(df, agg_cols, kf_cols, ga, gb)

    pf_oo = pf_frames["oo"].rename(
        columns={c: f"{c}_oo" for c in pf_frames["oo"].columns
                 if c not in agg_cols and c != "KF"}
    )
    pf_eo = pf_frames["eo"].rename(
        columns={c: f"{c}_eo" for c in pf_frames["eo"].columns
                 if c not in agg_cols and c != "KF"}
    )
    df_pf = pd.merge(pf_oo, pf_eo, on=agg_cols + ["KF"], how="outer")

    # Attach Hotelling p and q values (per group key, broadcast across KFs)
    hot_cols_oo = hotelling_frames["oo"].rename(
        columns={"p_value": "p_hotelling_oo", "q_value": "q_hotelling_oo"}
    )[agg_cols + ["p_hotelling_oo", "q_hotelling_oo"]]
    hot_cols_eo = hotelling_frames["eo"].rename(
        columns={"p_value": "p_hotelling_eo", "q_value": "q_hotelling_eo"}
    )[agg_cols + ["p_hotelling_eo", "q_hotelling_eo"]]

    df_pf = df_pf.merge(hot_cols_oo, on=agg_cols, how="left")
    df_pf = df_pf.merge(hot_cols_eo, on=agg_cols, how="left")
    df_pf = df_pf.sort_values(agg_cols + ["KF"])

    out2 = output_dir / "kidera_per_factor.tsv"
    df_pf.to_csv(out2, sep="\t", index=False)
    print(f"Written per-factor results ({len(df_pf)} rows)  → {out2}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Multivariate (Hotelling's T²) and per-factor (Welch t-test) "
            "comparisons of Kidera-factor profiles across genotype groups."
        )
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to the aggregated Kidera TSV produced by aa_kidera_aggregate.py.",
    )
    parser.add_argument(
        "--output_dir", default="sandbox/kidera_results",
        help="Directory for output TSV files. Default: sandbox/kidera_results.",
    )
    parser.add_argument(
        "--agg_columns", nargs="+",
        default=DEFAULT_AGG_COLUMNS,
        help=(
            "Columns defining each test group (cells is implicit if present). "
            "Default: IMGT_position cells."
        ),
    )
    parser.add_argument(
        "--kf_columns", nargs="+",
        default=DEFAULT_KF_COLUMNS,
        help="Kidera-factor column names. Default: KF1 … KF10.",
    )
    args = parser.parse_args()

    print(f"Reading: {args.input}")
    df = pd.read_csv(args.input, sep="\t")

    needed = args.agg_columns + args.kf_columns + ["genotype"]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        raise ValueError(f"Input is missing columns: {missing}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Aggregation columns : {args.agg_columns}")
    print(f"Kidera factor cols  : {args.kf_columns}")
    genotypes_present = df["genotype"].unique().tolist()
    print(f"Genotypes present   : {sorted(genotypes_present)}")

    run_tests(
        df=df,
        agg_cols=args.agg_columns,
        kf_cols=args.kf_columns,
        output_dir=output_dir,
    )


if __name__ == "__main__":
    main()
