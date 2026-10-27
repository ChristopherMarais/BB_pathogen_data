"""Profile the inputs (schemas, row counts, value counts) so the analysis scripts use the right columns."""
import glob
import json
from pathlib import Path

import duckdb
import pandas as pd

pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 60)
pd.set_option("display.max_colwidth", 80)
SNAP = Path(__file__).resolve().parents[2] / "data"
BB = Path(__file__).resolve().parents[2] / "data" / "literature_corpus"
ST = Path(__file__).resolve().parents[2] / "data"
A = Path(__file__).resolve().parents[2] / "data"


def show(name, df, n=3, vc=()):
    print(f"\n===== {name}: {df.shape}")
    print(df.dtypes.to_string())
    print(df.head(n).to_string()[:6000])
    for c in vc:
        if c in df.columns:
            print(f"-- value counts {c}")
            print(df[c].astype(str).value_counts().head(30).to_string())


# figshare files
for f in sorted(glob.glob(str(A / "figshare/*"))):
    print("\n##### figshare file", f)
    if f.endswith((".xlsx", ".xls")):
        x = pd.ExcelFile(f)
        for sh in x.sheet_names:
            d = x.parse(sh)
            show(f"{Path(f).name}:{sh}", d, 5)
    elif f.endswith(".csv"):
        show(Path(f).name, pd.read_csv(f), 5)

for t, vc in [("pathogens", ["predicate", "relation", "source_kind", "method", "partner_group", "category"]),
              ("notion_pathogens", ["category", "categories", "ecological_relationship"]),
              ("pathogen_refs", ["in_corpus", "status"]), ("pathogen_refs_fetch", ["status", "route"]),
              ("claims_checked", ["fact_type", "reason", "status"]), ("source_meta", ["kind", "meta_from"])]:
    p = SNAP / "tables" / f"{t}.parquet"
    if p.exists():
        show(t, pd.read_parquet(p), 3, vc)
print(open(SNAP / "tables/notion_pathogen_check.csv").read()[:3000])

con = duckdb.connect(str(SNAP / "kg/kg.duckdb"), read_only=True)
print(con.sql("select table_name, table_type from information_schema.tables order by 1"))
for v in ["enemies", "species", "names", "evidence", "places", "facts", "edges", "source_meta", "nodes"]:
    try:
        d = con.sql(f"select * from {v} limit 3").df()
        n = con.sql(f"select count(*) from {v}").fetchone()[0]
        print(f"\n===== view {v}: {n} rows")
        print(d.dtypes.to_string())
        print(d.to_string()[:5000])
    except Exception as e:
        print(v, "ERR", e)
print(con.sql("select predicate, count(*) from edges group by 1 order by 2 desc"))
print(con.sql("select method, count(*) from edges where predicate in ('has_natural_enemy','has_associated_organism') group by 1 order by 2 desc"))
print(con.sql("select predicate, partner_group, count(*) n from enemies group by 1,2 order by 3 desc limit 60"))

for p, vc in [(BB / "search_2026/candidates.parquet", ["kind", "source"]), (BB / "search_2026/hits_openalex.parquet", ["kind"]),
              (BB / "search_2026/hits_epmc.parquet", ["kind"]), (BB / "search_2026/hits_scopus.parquet", ["kind"]),
              (ST / "dataset/sources.parquet", ["kind", "source_kind", "redistributable", "language"])]:
    show(p.name, pd.read_parquet(p), 3, vc)
for p in sorted(glob.glob(str(ST / "work/sources_*.csv"))):
    d = pd.read_csv(p, nrows=3)
    n = sum(1 for _ in open(p)) - 1
    print(f"\n===== {Path(p).name}: ~{n} rows; cols {list(d.columns)}")
    print(d.head(2).to_string()[:1500])
print(open(BB / "download/processed_dois.txt").read()[:300])
print(pd.read_csv(BB / "download_2026/download_log.csv", nrows=3).to_string())
with open(BB / "abstracts_2026/abstracts.jsonl") as fh:
    for i, line in enumerate(fh):
        if i >= 2:
            break
        print(line[:800])
for d in ["paradigms", "onboarding", "risk", "priorities", "knowledge_structure"]:
    print("\n#####", d)
    for f in sorted(glob.glob(str(SNAP.parent.parent.parent / "analysis" / d / "out" / "**"), recursive=True))[:80]:
        print("  ", f)
