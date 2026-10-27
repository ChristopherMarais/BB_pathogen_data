"""Step 8 (analyses 7 and 8): where pathogen surveillance is most needed, and geographic coverage.

Species of concern (KG v2.22):
  introduced  = the species has an 'introduced' status in at least one country stated by a source (KG decision K-14:
                only when the quote or EPPO says so), medium / high confidence;
  regulated   = on an EPPO list (EPPO Global Database categorisation table of the KG);
  next invader = top 50 of the KG risk analysis forecast (species_kg/analysis/risk/out/list_q1_next_invaders.csv:
                ensemble of logistic regression, gradient boosting and a graph neural network trained on knowledge up
                to 2015; AUROC 0.91 for species first stated introduced 2016-2026).
For each: records in her dataset, in-scope pathogen / parasite records in the KG (step 6 rules R1-R4, any source, her
list excluded), and study effort = number of distinct literature works behind any KG fact about the species.
Ranking (short table): species of concern with no pathogen / parasite record in either, ordered by number of
countries where introduced, then regulated, then effort (a well-studied pest with no pathogen record is the clearest
gap; a little-studied one may simply be unstudied).
Geography: her records per country vs the number of introduced and of all species per country in the KG.
Outputs: out/surveillance_species.csv, out/surveillance_top.csv, out/geography_countries.csv
"""
import json
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pycountry

A = Path("/home/gmarais/gmarais/PhD/BarkBeetle_KG/papers/pathogens/analysis")
SNAP = Path("/home/gmarais/gmarais/PhD/BarkBeetle_KG/species_kg/snapshots/kg_v2.22_2026-10-22/species_kg")
RISK = Path("/home/gmarais/gmarais/PhD/BarkBeetle_KG/species_kg/analysis/risk/out")


def main():
    con = duckdb.connect(str(SNAP / "kg/kg.duckdb"), read_only=True)
    SP = con.sql("select species_id, name, genus, tribe, subfamily from species").df()
    PL = con.sql("""select taxon_id, place, country_iso2, statuses, level from places
                    where taxon_rank='species' and predicate='occurs_in'""").df()
    intro = PL[PL.statuses.fillna("").str.contains("introduced") & PL.level.isin(["medium", "high"])]
    n_intro = intro.groupby("taxon_id").country_iso2.nunique()
    print("species introduced somewhere (medium/high):", len(n_intro))
    EC = pd.read_parquet(SNAP / "tables/eppo_categorization.parquet")
    print("EPPO categorisation columns:", list(EC.columns), len(EC))
    sid_col = next((c for c in EC.columns if "species_id" in c or c in ("subject", "taxon_id")), None)
    regulated = set(EC[sid_col].dropna().astype(str)) if sid_col else set()
    if sid_col and not next(iter(regulated), "").startswith("species:"):
        regulated = {f"species:{int(float(x))}" for x in regulated if re.match(r"^\d+(\.0)?$", x)}
    print("EPPO-listed species:", len(regulated))
    NI = pd.read_csv(RISK / "list_q1_next_invaders.csv")
    top50 = set(NI.sort_values("rank_ensemble").head(50).species_id)
    effort = con.sql("""select subject, count(distinct source_id) works from edges
                        where method='llm_facts+checks' and subject like 'species:%' group by 1""").df()
    eff = dict(zip(effort.subject, effort.works))
    R = pd.read_parquet(A / "out/her_records.parquet")
    TB = pd.read_csv(A / "out/taxonomy_beetles.csv")
    sid_of = {r.accepted_in_dataset: str(r.kg_species_id).split(";")[0] for r in TB.itertuples()
              if isinstance(r.kg_species_id, str) and r.kg_species_id}
    R["sid"] = R.AcceptedBeetleSpeciesName.map(sid_of)
    her_n = R.groupby("sid").size()
    her_pubs = R.groupby("sid").ref.nunique()
    K = pd.read_csv(A / "data/kg_inscope_records.csv")
    K = K[K.subject.str.startswith("species:")]
    kg_n = K.groupby("subject").source_id.nunique()
    kg_n_indep = K[~K.in_her].groupby("subject").source_id.nunique()
    S = SP.copy()
    S["introduced_countries"] = S.species_id.map(n_intro).fillna(0).astype(int)
    S["eppo_listed"] = S.species_id.isin(regulated)
    S["next_invader_top50"] = S.species_id.isin(top50)
    S["study_effort_works"] = S.species_id.map(eff).fillna(0).astype(int)
    S["her_records"] = S.species_id.map(her_n).fillna(0).astype(int)
    S["her_publications"] = S.species_id.map(her_pubs).fillna(0).astype(int)
    S["kg_pathogen_sources"] = S.species_id.map(kg_n).fillna(0).astype(int)
    S["kg_pathogen_sources_not_in_her_data"] = S.species_id.map(kg_n_indep).fillna(0).astype(int)
    S["concern"] = (S.introduced_countries > 0) | S.eppo_listed | S.next_invader_top50
    S["no_pathogen_record"] = (S.her_records == 0) & (S.kg_pathogen_sources == 0)
    S.to_csv(A / "out/surveillance_species.csv", index=False)
    c = S[S.concern]
    print("species of concern:", len(c), "| introduced", int((c.introduced_countries > 0).sum()), "| EPPO",
          int(c.eppo_listed.sum()), "| next-invader top 50", int(c.next_invader_top50.sum()))
    for nm, m in [("introduced", c.introduced_countries > 0), ("EPPO-listed", c.eppo_listed),
                  ("next-invader top 50", c.next_invader_top50), ("any concern", c.index == c.index)]:
        g = c[m]
        print(f"{nm}: {len(g)} species | in her dataset {int((g.her_records > 0).sum())} | KG pathogen record only "
              f"{int(((g.her_records == 0) & (g.kg_pathogen_sources > 0)).sum())} | none {int(g.no_pathogen_record.sum())}"
              f" | median effort with records {g[~g.no_pathogen_record].study_effort_works.median():.0f} vs none "
              f"{g[g.no_pathogen_record].study_effort_works.median():.0f}")
    T = c[c.no_pathogen_record].sort_values(["introduced_countries", "eppo_listed", "study_effort_works"],
                                            ascending=False)
    T = T[["name", "tribe", "subfamily", "introduced_countries", "eppo_listed", "next_invader_top50",
           "study_effort_works"]]
    T.head(40).to_csv(A / "out/surveillance_top.csv", index=False)
    pd.set_option("display.width", 220)
    print(T.head(25).to_string())
    # also: top species by introduced countries overall, with their record counts
    print(c.sort_values("introduced_countries", ascending=False).head(25)[
        ["name", "introduced_countries", "eppo_listed", "study_effort_works", "her_records", "kg_pathogen_sources",
         "kg_pathogen_sources_not_in_her_data"]].to_string())
    # next-invader list rows
    print(c[c.next_invader_top50].sort_values("study_effort_works", ascending=False)[
        ["name", "introduced_countries", "study_effort_works", "her_records", "kg_pathogen_sources"]].head(20).to_string())

    # ---------------- geography
    def iso(name):
        if pd.isna(name):
            return ""
        n = str(name or "").strip()
        if not n or n.lower() in ("nan", "none", "<na>", "na"):
            return ""
        fix = {"United States": "US", "USA": "US", "Russia": "RU", "Czech Republic": "CZ", "Czechia": "CZ",
               "South Korea": "KR", "Korea": "KR", "Iran": "IR", "Turkey": "TR", "Türkiye": "TR", "UK": "GB",
               "United Kingdom": "GB", "England": "GB", "Great Britain": "GB", "Scotland": "GB",
               "British Columbia": "CA", "México": "MX", "Mexico": "MX", "Côte d’Ivoire": "CI",
               "Côte d'Ivoire": "CI", "Ivory Coast": "CI",
               "Taiwan": "TW", "Vietnam": "VN",
               "Bosnia and Herzegovina": "BA", "Slovak Republic": "SK", "Moldova": "MD", "Laos": "LA",
               "Former Yugoslavia": "", "Yugoslavia": "", "USSR": "", "Soviet Union": "", "Europe": ""}
        if n in fix:
            return fix[n]
        try:
            return pycountry.countries.lookup(n).alpha_2
        except Exception:
            try:
                return pycountry.countries.search_fuzzy(n)[0].alpha_2
            except Exception:
                return ""
    R["countries"] = R["country or region"].fillna("").map(lambda s: [x.strip() for x in re.split(r"[;,/]| and ", s)
                                                                       if x.strip()])
    X = R.explode("countries")
    X["iso2"] = X.countries.map(iso)
    print("her country strings not mapped:", X[X.iso2.eq("") & X.countries.notna()].countries.value_counts().head(20).to_dict())
    her_c = X[X.iso2.ne("")].groupby("iso2").agg(her_records=("pathogens", "size"),
                                                  her_beetle_species=("sid", "nunique"),
                                                  her_publications=("ref", "nunique"))
    allc = PL[PL.level.isin(["medium", "high"])].groupby("country_iso2").taxon_id.nunique().rename("kg_species")
    intc = intro.groupby("country_iso2").taxon_id.nunique().rename("kg_introduced_species")
    G = pd.concat([allc, intc, her_c], axis=1).fillna(0)
    G.index.name = "iso2"
    G = G.reset_index()
    G["country"] = G.iso2.map(lambda c: getattr(pycountry.countries.get(alpha_2=c), "name", c) if c else "")
    G = G.sort_values("kg_introduced_species", ascending=False)
    G.to_csv(A / "out/geography_countries.csv", index=False)
    print("countries with her records:", int((G.her_records > 0).sum()))
    print(G.head(30).to_string())
    print("countries with >= 10 introduced species and no record in her dataset:")
    print(G[(G.kg_introduced_species >= 10) & (G.her_records == 0)].to_string())


if __name__ == "__main__":
    main()
