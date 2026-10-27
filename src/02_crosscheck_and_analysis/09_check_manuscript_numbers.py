"""Step 9: recompute every number stated in the manuscript text from her own figshare file (v1.1), so that any
disagreement can be fixed before submission. Output: out/manuscript_number_check.csv"""
from pathlib import Path

import pandas as pd

A = Path(__file__).resolve().parents[2] / "data"


def main():
    F = pd.read_csv(A / "figshare/bark_ambrosia_beetle_pathogen_database_v1.1.csv", encoding="utf-8-sig")
    F["ref"] = F["title"].fillna("").str.replace(r"\s+", " ", regex=True).str.strip()
    F["pathogens"] = F.pathogens.fillna("").str.strip()
    R = F[F.pathogens.ne("") & F.categories.notna()].copy()
    TB = pd.read_csv(A / "out/taxonomy_beetles.csv")
    species_level = TB[~TB.accepted_in_dataset.str.contains(r"\bspp?\.|Scolytidae", regex=True)]
    rows = []

    def add(claim, stated, value, note=""):
        rows.append(dict(claim=claim, stated_in_text=stated, recomputed=value,
                         agrees=str(stated).replace(",", "") == str(value).replace(",", ""), note=note))
    add("records in final dataset", 1010, len(R))
    add("rows in figshare file (incl. records without a taxon)", "1,010 + 4", len(F))
    add("primary publications (distinct full references, analytical records)", 280, R.ref.nunique(),
        "distinct references over all 1,014 rows: %d; distinct short citations: %d" % (F.ref.nunique(),
                                                                                       R.source.nunique()))
    add("unique candidate pathogen / parasite taxa (distinct names)", 274, R.pathogens.nunique(),
        "case-insensitive: %d" % R.pathogens.str.lower().nunique())
    add("beetle species (distinct accepted names)", 109, R.AcceptedBeetleSpeciesName.nunique(),
        "of which named to species: %d; not species-level (sp./spp./Scolytidae sp.): %d" % (
            R[~R.AcceptedBeetleSpeciesName.str.contains(r"\bspp?\.|Scolytidae", regex=True)].AcceptedBeetleSpeciesName.nunique(),
            R[R.AcceptedBeetleSpeciesName.str.contains(r"\bspp?\.|Scolytidae", regex=True)].AcceptedBeetleSpeciesName.nunique()))
    add("genera", 27, R.genus.nunique(), "genus values: " + "; ".join(sorted(R.genus.dropna().unique())))
    add("year range", "1914-2026", f"{R.year.min()}-{R.year.max()}")
    for cat, (taxa, recs) in {"Nematode": (121, 220), "Fungi": (59, 407), "Virus": (18, 46)}.items():
        g = R[R.categories.eq(cat)]
        add(f"{cat}: taxa", taxa, g.pathogens.nunique())
        add(f"{cat}: records", recs, len(g))
    for cat in ["Bacteria", "Microsporidia", "Protist"]:
        g = R[R.categories.eq(cat)]
        add(f"{cat}: taxa (not stated in text)", "", g.pathogens.nunique())
        add(f"{cat}: records (not stated in text)", "", len(g))
    v = R["validation type"].fillna("").str.strip()
    add("records with any experimental validation", 238, int(v.ne("").sum()),
        "values: " + "; ".join(f"{k}: {n}" for k, n in v[v.ne("")].value_counts().items()))
    add("share with validation (%)", 23.56, round(100 * v.ne("").mean(), 2))
    top = R.AcceptedBeetleSpeciesName.value_counts()
    add("Ips typographus records", 253, int(top.get("Ips typographus", 0)))
    add("top 10 beetle species", "", "; ".join(f"{k} ({n})" for k, n in top.head(10).items()))
    for cat, stated in {"Bacteria": "Serratia marcescens", "Fungi": "Beauveria bassiana",
                        "Microsporidia": "Chytridiopsis typographi", "Nematode": "Contortylenchus brevicomi",
                        "Protist": "Gregarina typographi", "Virus": "ItEPV"}.items():
        vc = R[R.categories.eq(cat)].pathogens.value_counts()
        add(f"most frequent {cat}", stated, vc.index[0], "top 3: " + "; ".join(f"{k} ({n})" for k, n in vc.head(3).items()))
    df = R[R.AcceptedBeetleSpeciesName.eq("Dendroctonus frontalis")].groupby("ref").year.first()
    add("D. frontalis source publications", 16, len(df))
    add("D. frontalis publications in 1970s-1980s", 12, int(df.between(1970, 1989).sum()))
    dec = R.groupby("ref").year.first().floordiv(10).mul(10).value_counts().sort_index()
    add("publications per decade (her Fig. 4)", "peak 2010s", "; ".join(f"{k}s: {n}" for k, n in dec.items()),
        "peak decade: %ds" % dec.idxmax())
    add("Platypodinae records", 0, int(R.subfamily.eq("Platypodinae").sum()))
    ft = R["full text"].fillna("(blank)").str.strip()
    add("full-text flag values", "", "; ".join(f"{k}: {n}" for k, n in ft.value_counts().items()),
        "spelling variants 'non full text' / 'non full-text' and blanks could be unified")
    er = R["ecological relationship"].fillna("(blank)").str.strip().str.lower().value_counts()
    add("ecological relationship values", "pathogen, parasite, endoparasite, associate",
        "; ".join(f"{k}: {n}" for k, n in er.items()))
    O = pd.DataFrame(rows)
    O.to_csv(A / "out/manuscript_number_check.csv", index=False)
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 200)
    print(O.to_string())


if __name__ == "__main__":
    main()
