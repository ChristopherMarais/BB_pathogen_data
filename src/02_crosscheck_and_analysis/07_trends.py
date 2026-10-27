"""Step 7 (analysis 6, Jiri's request): publication trend rescaled by the size of the bark and ambrosia beetle literature.

Numerator: her source publications per period (all distinct references of the figshare file, the 4 records without a
taxon included, as in her Figure 4).
Denominators (each with its bias, stated in the output):
  D1 OpenAlex corpus: works whose title or abstract names a bark or ambrosia beetle (topic words, a genus or a binomial
     of the species list), article-like types, all languages; 7 OpenAlex title/abstract searches made for the KG's
     paradigm analysis (species_kg/analysis/paradigms/cache/corpus.parquet). Biases: OpenAlex indexes recent, English,
     DOI-bearing journal work best (pre-1970 works, theses, regional journals under-counted, which inflates the early
     ratios); abstracts are missing for half of pre-1970 works, so old works that name a beetle only in the abstract
     are missed.
  D2 KG reading corpus: distinct literature works (papers, abstracts, scans) from which the KG extracted at least one
     fact about a bark or ambrosia beetle (hosts, distribution, size, fungi, enemies ...), by publication year.
     Biases: the KG's harvest was species-name and taxonomy driven, so taxonomic and faunistic work is over-represented
     and the last 2 years are incomplete in its abstract store.
Rate = numerator / denominator, per 1,000 works, with a Wilson 95 % interval (treating the numerator as a share of the
period's literature). 5-year bins and decades. The 2020s run from 2020 to 2026 (her search ended 24 July 2026).
Test of the latest change: share in 2020-2026 vs 2010-2019, two-proportion z test, and a binomial regression of the
yearly share on year (logit link, 1990-2025) fitted by iteratively reweighted least squares (numpy).
Check with the KG's own pathogen literature: primary sources the KG read that hold an in-scope pathogen / parasite
record (step 6), divided by D2.
Outputs: out/trend_rescaled_decade.csv, out/trend_rescaled_5yr.csv, out/trend_tests.txt
"""
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

A = Path(__file__).resolve().parents[2] / "data"
SNAP = Path(__file__).resolve().parents[2] / "data"
PAR = Path(__file__).resolve().parents[2] / "data"


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def logit_irls(x, k, n, iters=50):
    X = np.column_stack([np.ones_like(x), x])
    b = np.zeros(2)
    for _ in range(iters):
        eta = X @ b
        p = 1 / (1 + np.exp(-eta))
        W = n * p * (1 - p)
        z = eta + (k - n * p) / np.maximum(W, 1e-12)
        XtW = X.T * W
        b_new = np.linalg.solve(XtW @ X, XtW @ z)
        if np.max(np.abs(b_new - b)) < 1e-10:
            b = b_new
            break
        b = b_new
    p = 1 / (1 + np.exp(-(X @ b)))
    cov = np.linalg.inv((X.T * (n * p * (1 - p))) @ X)
    # Pearson overdispersion
    chi2 = np.sum((k - n * p) ** 2 / np.maximum(n * p * (1 - p), 1e-12))
    phi = max(1.0, chi2 / (len(x) - 2))
    return b, np.sqrt(np.diag(cov) * phi), phi


def main():
    F = pd.read_csv(A / "figshare/bark_ambrosia_beetle_pathogen_database_v1.1.csv", encoding="utf-8-sig")
    F["ref"] = F["title"].fillna("").str.replace(r"\s+", " ", regex=True).str.strip()
    N = F.groupby("ref").year.min().reset_index()
    print("her publications (all distinct references incl. records without a taxon):", len(N),
          "| analytical publications:", F[F.pathogens.fillna("").str.strip().ne("") & F.categories.notna()].ref.nunique())
    C = pd.read_parquet(PAR / "corpus.parquet", columns=["id", "year", "language", "type"])
    print("D1 OpenAlex corpus works:", len(C), "| years", C.year.min(), C.year.max())
    con = duckdb.connect(str(SNAP / "kg/kg.duckdb"), read_only=True)
    K = con.sql("""select distinct e.source_id, m.year from edges e join source_meta m on m.source_id = e.source_id
                   where e.method = 'llm_facts+checks'""").df()
    K["year"] = pd.to_numeric(K.year, errors="coerce")
    K = K.dropna(subset=["year"])
    K = K[K.year.between(1750, 2026)]
    print("D2 KG reading corpus works with a year:", len(K))
    KP = pd.read_csv(A / "out/kg_pathogen_sources_primary.csv")
    KP["year"] = pd.to_numeric(KP.year, errors="coerce")
    KP = KP.dropna(subset=["year"])

    def table(binner, label):
        rows = []
        for b, lo, hi in binner:
            n_her = int(N.year.between(lo, hi).sum())
            d1 = int(C.year.between(lo, hi).sum())
            d2 = int(K.year.between(lo, hi).sum())
            kp = int(KP.year.between(lo, hi).sum())
            r = dict(period=b, start=lo, end=hi, her_publications=n_her, openalex_bark_beetle_works=d1,
                     kg_reading_corpus_works=d2, kg_primary_pathogen_sources=kp)
            for nm, den in [("openalex", d1), ("kg", d2)]:
                r[f"per1000_{nm}"] = 1000 * n_her / den if den else np.nan
                lo_, hi_ = wilson(n_her, den)
                r[f"per1000_{nm}_lo"], r[f"per1000_{nm}_hi"] = 1000 * lo_, 1000 * hi_
            r["kg_pathogen_per1000_kg"] = 1000 * kp / d2 if d2 else np.nan
            r["kg_pathogen_per1000_kg_lo"], r["kg_pathogen_per1000_kg_hi"] = [1000 * v for v in wilson(kp, d2)]
            rows.append(r)
        T = pd.DataFrame(rows)
        T.to_csv(A / f"out/trend_rescaled_{label}.csv", index=False)
        pd.set_option("display.width", 250)
        print(T.round(2).to_string())
        return T
    dec = [("pre-1950", 1900, 1949)] + [(f"{d}s", d, d + 9) for d in range(1950, 2020, 10)] + [("2020-2026", 2020, 2026)]
    T10 = table(dec, "decade")
    five = [("pre-1950", 1900, 1949)] + [(f"{y}-{y + 4}", y, y + 4) for y in range(1950, 2020, 5)] + \
        [("2020-2024", 2020, 2024), ("2025-2026", 2025, 2026)]
    T5 = table(five, "5yr")
    # tests
    lines = []
    a = T10.set_index("period")
    for nm, col in [("OpenAlex corpus", "openalex_bark_beetle_works"), ("KG reading corpus", "kg_reading_corpus_works")]:
        k1, n1 = a.loc["2010s", "her_publications"], a.loc["2010s", col]
        k2, n2 = a.loc["2020-2026", "her_publications"], a.loc["2020-2026", col]
        p1, p2 = k1 / n1, k2 / n2
        pp = (k1 + k2) / (n1 + n2)
        z = (p2 - p1) / np.sqrt(pp * (1 - pp) * (1 / n1 + 1 / n2))
        from math import erf, sqrt
        pval = 2 * (1 - 0.5 * (1 + erf(abs(z) / sqrt(2))))
        lines.append(f"{nm}: 2010s {k1}/{n1} = {1000 * p1:.2f} per 1,000; 2020-2026 {k2}/{n2} = {1000 * p2:.2f} per "
                     f"1,000; ratio {p2 / p1:.2f}; two-proportion z = {z:.2f}, p = {pval:.3g}")
        # peak period of the rescaled share
        s = a[f"per1000_{'openalex' if 'OpenAlex' in nm else 'kg'}"]
        lines.append(f"   rescaled share by decade (per 1,000): " + "; ".join(f"{i} {v:.1f}" for i, v in s.items()))
    # yearly logit trend 1990-2025
    yrs = np.arange(1990, 2026)
    k = np.array([int((N.year == y).sum()) for y in yrs])
    for nm, den_s in [("OpenAlex corpus", C.year), ("KG reading corpus", K.year)]:
        n = np.array([int((den_s == y).sum()) for y in yrs])
        m = n > 0
        b, se, phi = logit_irls((yrs[m] - 2000) / 10.0, k[m], n[m])
        lines.append(f"{nm}: yearly share 1990-2025, logit slope per decade = {b[1]:.3f} (SE {se[1]:.3f}, "
                     f"overdispersion {phi:.2f}); odds ratio per decade {np.exp(b[1]):.2f} "
                     f"[{np.exp(b[1] - 1.96 * se[1]):.2f}, {np.exp(b[1] + 1.96 * se[1]):.2f}]")
        b2, se2, phi2 = logit_irls((yrs[m & (yrs >= 2010)] - 2010) / 10.0, k[m & (yrs >= 2010)], n[m & (yrs >= 2010)])
        lines.append(f"   2010-2025 only: odds ratio per decade {np.exp(b2[1]):.2f} [{np.exp(b2[1] - 1.96 * se2[1]):.2f}, "
                     f"{np.exp(b2[1] + 1.96 * se2[1]):.2f}]")
    lines.append("raw counts by year 2005-2026: " + "; ".join(f"{y}:{int((N.year == y).sum())}" for y in range(2005, 2027)))
    lines.append("OpenAlex corpus by year 2015-2026: " + "; ".join(f"{y}:{int((C.year == y).sum())}" for y in range(2015, 2027)))
    (A / "out/trend_tests.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
