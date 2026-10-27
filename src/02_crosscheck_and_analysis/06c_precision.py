"""Step 6c: precision of the candidate list, from manual reading (by the analyst, against her inclusion rules) of
  (1) a random sample of 50 candidate records (unique publication x beetle x organism; drawn 2026-10-23, seed 20261023,
      frozen in data/INTERNAL_sample50_drawn_20261023.csv), and
  (2) a random sample of 40 candidate publications of the focused subset (seed 20261024).
Judgement 'yes' = the publication itself reports the organism as a pathogen / parasite of that beetle (own data or
own observation), fungi / bacteria explicitly as pathogens, primary source; 'borderline' = isolation source or review
of own work; 'no' = reference-list title, citation of another work, passing mention, associate only, wrong pair.
Precision = yes / judged, Wilson 95 % interval; 'lenient' counts borderline as yes. Estimated genuine count =
precision x population (interval from the Wilson bounds).
Output: out/candidate_precision.csv; appends the focused step to out/candidate_funnel.csv.
"""
from pathlib import Path

import numpy as np
import pandas as pd

A = Path(__file__).resolve().parents[2] / "data"


def wilson(k, n, z=1.96):
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return c - h, c + h


def main():
    C = pd.read_csv(A / "out/candidate_records_public.csv")
    L = C[C.stratum.eq("literature read by the KG") & C.beetle_rank.eq("species")]
    U = L.drop_duplicates(["source_id", "beetle", "p_name"])
    J = pd.read_csv(A / "data/manual_judgments_sample50.csv")
    D = pd.read_csv(A / "data/INTERNAL_sample50_drawn_20261023.csv")[["edge_id"]].merge(J, on="edge_id")
    D = D.merge(L[["edge_id", "focused"]], on="edge_id", how="left")
    dropped = D.focused.isna()
    print("sample records no longer candidates (source now matched to her file):", int(dropped.sum()))
    D = D[~dropped]
    rows = []

    def add(name, d, pop, unit):
        n = len(d)
        k = int(d.judgement.eq("yes").sum())
        kl = int(d.judgement.isin(["yes", "borderline"]).sum())
        lo, hi = wilson(k, n)
        llo, lhi = wilson(kl, n)
        rows.append(dict(estimate=name, unit=unit, judged=n, genuine=k, borderline=kl - k, genuine_pct=100 * k / n,
                         ci_lo=100 * lo, ci_hi=100 * hi, lenient_pct=100 * kl / n, lenient_ci_lo=100 * llo,
                         lenient_ci_hi=100 * lhi, population=pop, est_genuine=round(pop * k / n),
                         est_genuine_lo=round(pop * lo), est_genuine_hi=round(pop * hi)))
    add("records, all candidates", D, len(U), "publication x beetle x organism")
    add("records, focused subset", D[D.focused.astype(bool)], int(U.focused.sum()), "publication x beetle x organism")
    add("records, outside focused subset", D[~D.focused.astype(bool)], int((~U.focused).sum()),
        "publication x beetle x organism")
    jp = A / "data/manual_judgments_sample40_publications.csv"
    if jp.exists():
        JP = pd.read_csv(jp)
        add("publications, focused subset", JP, U[U.focused].source_id.nunique(), "publication")
        Fo = U[U.focused]
        print("focused population: publications", Fo.source_id.nunique(), "| new-pair publications",
              Fo[Fo.new_pair].source_id.nunique())
        S2 = pd.read_csv(A / "out/INTERNAL_candidate_sample40_publications.csv")
        cbb_src = set(S2[S2.beetle.eq("Hypothenemus hampei")].source_id)
        JP["cbb"] = JP.source_id.isin(cbb_src)
        print("sample ids matched:", int(JP.source_id.isin(set(S2.source_id)).sum()), "of", len(JP))
        cbb_pop = Fo[Fo.beetle.eq("Hypothenemus hampei")].source_id.nunique()
        add("publications on coffee berry borer, focused subset", JP[JP.cbb], cbb_pop, "publication")
        add("publications on other beetles, focused subset", JP[~JP.cbb], Fo.source_id.nunique() - cbb_pop, "publication")
        R = pd.read_parquet(A / "out/her_records.parquet")
        print("her publications on H. hampei:", R[R.AcceptedBeetleSpeciesName.eq("Hypothenemus hampei")].ref.nunique(),
              "| KG focused candidate publications on H. hampei:", cbb_pop)
        if "new_pair_genuine" in JP.columns:
            add("publications with a genuine new pair, focused subset",
                JP.assign(judgement=np.where(JP.new_pair_genuine.eq("yes"), "yes", "no")),
                U[U.focused].source_id.nunique(), "publication")
    O = pd.DataFrame(rows)
    O.to_csv(A / "out/candidate_precision.csv", index=False)
    pd.set_option("display.width", 250)
    print(O.round(1).to_string())
    print("focused subset keeps", int(D[D.judgement.eq("yes")].focused.astype(bool).sum()), "of",
          int(D.judgement.eq("yes").sum()), "genuine sample records")
    F = pd.read_csv(A / "out/candidate_funnel.csv")
    F = F[~F.step.eq("R6 focused")]
    F = pd.concat([F, pd.DataFrame([dict(step="R6 focused", label="R6 on-topic source, own statement",
                                         sources=L[L.focused].source_id.nunique())])], ignore_index=True)
    F.to_csv(A / "out/candidate_funnel.csv", index=False)


if __name__ == "__main__":
    main()
