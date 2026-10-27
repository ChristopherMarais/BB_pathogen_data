"""Step 11: short tables for the report.
  out/surveillance_ranked_short.csv : species of concern (introduced, EPPO-listed or forecast next invader) with at most
      2 records in her dataset, ranked by number of countries where introduced, then study effort.
  out/taxonomy_name_merge.csv       : her 274 pathogen names grouped into taxa when they differ only by 'cf.', case,
      or are GBIF synonyms of the same accepted species (suggestion: how many distinct taxa the 274 names represent).
"""
import re
from pathlib import Path

import pandas as pd

A = Path(__file__).resolve().parents[2] / "data"


def main():
    S = pd.read_csv(A / "out/surveillance_species.csv")
    c = S[S.concern & (S.her_records <= 2)].sort_values(["introduced_countries", "study_effort_works"], ascending=False)
    cols = ["name", "tribe", "subfamily", "introduced_countries", "eppo_listed", "next_invader_top50",
            "study_effort_works", "her_records", "kg_pathogen_sources", "kg_pathogen_sources_not_in_her_data"]
    c[cols].head(20).to_csv(A / "out/surveillance_ranked_short.csv", index=False)
    pd.set_option("display.width", 250)
    print(c[cols].head(20).to_string())
    T = pd.read_csv(A / "out/taxonomy_pathogens.csv")

    def key(r):
        n = re.sub(r"\b(cf|aff)\.\s*", "", str(r.pathogen_in_dataset)).strip().lower()
        if isinstance(r.gbif_accepted_species, str) and r.gbif_rank in ("SPECIES", "SUBSPECIES", "VARIETY") and \
                r.gbif_match_type in ("EXACT", "FUZZY") and not str(r.flag).startswith("GBIF matched a different spelling: Parat"):
            return "gbif:" + r.gbif_accepted_species.lower()
        return n
    T["taxon_key"] = T.apply(key, axis=1)
    g = T.groupby("taxon_key").agg(names=("pathogen_in_dataset", lambda s: "; ".join(sorted(s))),
                                   n_names=("pathogen_in_dataset", "size"), records=("n_records", "sum"))
    g[g.n_names > 1].to_csv(A / "out/taxonomy_name_merge.csv")
    vague = T.pathogen_in_dataset.str.contains(r"unidentified|unknown|nematodes?$|Nematodan|-like|similar to|"
                                                r"granulosis|polyhedrosis|nonoccluded|un identified|larva|Rod-like|"
                                                r"^gut|^hemolymph|^intestinal", case=False, regex=True)
    print("distinct names", len(T), "| taxon groups after merging cf. / case / GBIF synonyms", g.shape[0],
          "| names that are not taxa (vague descriptions)", int(vague.sum()))
    print("sp.-level names:", int(T.pathogen_in_dataset.str.contains(r"\bspp?\d*\.", regex=True).sum()))
    print(g[g.n_names > 1].to_string())
    print("vague names:", "; ".join(T[vague].pathogen_in_dataset))


if __name__ == "__main__":
    main()
