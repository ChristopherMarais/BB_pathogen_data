"""Step 10: publication-quality figures for the pathogen dataset cross-check and trend analysis (Figures 5-7),
styled to match the manuscript's Figures 1-4 (clean sans-serif typography, blue palette, value labels that
match the plotted bar heights and y/x-axis scales, light grid, no top/right spines).
Each figure is saved as PNG (300 dpi) and PDF in analysis/figures/ and figures_for_authors/.
"""
from pathlib import Path
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

A = Path(__file__).resolve().parents[2] / "data"
FIG = Path(__file__).resolve().parents[2] / "figures"
FIG.mkdir(exist_ok=True)
AUTHORS_FIG = A.parent / "figures_for_authors"

BLUE, LIGHT, GREY, ORANGE, INK = "#4C72B0", "#A9BCD0", "#C9CED6", "#DD8452", "#333333"
plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.spines.left": False,
    "axes.edgecolor": "#BBBBBB",
    "axes.grid": True,
    "axes.grid.axis": "y",
    "grid.color": "#E6E6E6",
    "grid.linewidth": 0.8,
    "axes.axisbelow": True,
    "xtick.color": INK,
    "ytick.color": "#444444",
    "axes.titleweight": "bold",
    "axes.titlesize": 10,
    "axes.titlelocation": "left",
})


def wilson_ci(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return max(0.0, 100.0 * (c - h)), min(100.0, 100.0 * (c + h))


def save(fig, name):
    png_path = FIG / f"{name}.png"
    pdf_path = FIG / f"{name}.pdf"
    fig.savefig(png_path, dpi=300, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    plt.close(fig)
    if AUTHORS_FIG.exists() and name in ("Figure5_crosscheck", "Figure6_rescaled_trend", "Figure7_candidates"):
        shutil.copy2(png_path, AUTHORS_FIG / f"{name}.png")
        shutil.copy2(pdf_path, AUTHORS_FIG / f"{name}.pdf")
    print("saved", name)


def fig5():
    S = pd.read_csv(A / "out/recall_summary.csv").set_index("group")
    AG = pd.read_csv(A / "out/agreement_summary.csv")
    fig, axs = plt.subplots(1, 3, figsize=(16.2, 5.2), gridspec_kw={"width_ratios": [1.22, 1.12, 1.26]})

    # ── Panel A: Search recall by decade (including 1920s = 0 to match Figure 4) ──
    ax = axs[0]
    decades = [1910, 1920, 1930, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020]
    lab = [f"{d}s" for d in decades]
    n_vals, k_vals = [], []
    for d in decades:
        g = f"decade {d}"
        if g in S.index:
            n_vals.append(int(S.loc[g, "n_publications"]))
            k_vals.append(int(S.loc[g, "found_kg_search"]))
        else:
            n_vals.append(0)
            k_vals.append(0)
    n = np.array(n_vals)
    k = np.array(k_vals)
    x = np.arange(len(decades))
    ax.bar(x, k, width=0.72, color=BLUE, label="retrieved by automated search (k)")
    ax.bar(x, n - k, bottom=k, width=0.72, color=GREY, label="not retrieved (n - k)", edgecolor="white", linewidth=0.8)
    for i in range(len(x)):
        if n[i] > 0:
            pct = 100.0 * k[i] / n[i]
            ax.text(x[i], n[i] + 1.0, f"{k[i]}/{n[i]}\n({pct:.0f}%)",
                    ha="center", va="bottom", fontsize=7.0, color=INK, linespacing=1.05)
        else:
            ax.text(x[i], 0.8, "0", ha="center", va="bottom", fontsize=7.0, color="#777777")
    ax.set_xticks(x, lab, rotation=40, ha="right", fontsize=8)
    ax.set_ylim(0, 90)
    ax.set_yticks(range(0, 81, 10))
    ax.set_ylabel("Dataset source publications (n = 276)")
    ax.set_title("A  Search recall by decade")
    ax.legend(frameon=False, loc="upper left", fontsize=8)

    # ── Panel B: Recall by DOI status, publication type, and language ──
    ax = axs[1]
    other_langs = ["language Russian", "language Polish", "language Czech"]
    chinese_langs = [i for i in S.index if i.startswith("language Chinese")]

    items = [
        ("All sources", int(S.loc["all", "n_publications"]), int(S.loc["all", "found_kg_search"])),
        ("with article DOI", int(S.loc["has article DOI", "n_publications"]), int(S.loc["has article DOI", "found_kg_search"])),
        ("no article DOI", int(S.loc["no article DOI (none or CAB Abstracts record only)", "n_publications"]),
         int(S.loc["no article DOI (none or CAB Abstracts record only)", "found_kg_search"])),
        ("journal article", int(S.loc["type journal article", "n_publications"]), int(S.loc["type journal article", "found_kg_search"])),
        ("report / proceedings", int(S.loc["type report / proceedings", "n_publications"]), int(S.loc["type report / proceedings", "found_kg_search"])),
        ("book / chapter", int(S.loc["type book / chapter", "n_publications"]), int(S.loc["type book / chapter", "found_kg_search"])),
        ("thesis", int(S.loc["type thesis", "n_publications"]), int(S.loc["type thesis", "found_kg_search"])),
        ("English", int(S.loc["language English", "n_publications"]), int(S.loc["language English", "found_kg_search"])),
        ("German", int(S.loc["language German", "n_publications"]), int(S.loc["language German", "found_kg_search"])),
        ("French", int(S.loc["language French", "n_publications"]), int(S.loc["language French", "found_kg_search"])),
        ("Chinese (inferred)", int(S.loc[chinese_langs, "n_publications"].sum()), int(S.loc[chinese_langs, "found_kg_search"].sum())),
        ("Russian / Polish / Czech", int(S.loc[other_langs, "n_publications"].sum()), int(S.loc[other_langs, "found_kg_search"].sum())),
    ]

    names, vals, lo, hi = [], [], [], []
    for nm, nn, kk in items:
        p = 100.0 * kk / nn if nn else 0.0
        l_ci, h_ci = wilson_ci(kk, nn)
        names.append(f"{nm} ({kk}/{nn})")
        vals.append(p)
        lo.append(l_ci)
        hi.append(h_ci)

    y = np.arange(len(names))[::-1]
    bar_colors = [INK] + [BLUE] * (len(names) - 1)
    ax.barh(y, vals, color=bar_colors, height=0.62, label="recall (%)")
    ax.errorbar(vals, y, xerr=[np.maximum(0, np.array(vals) - lo), np.maximum(0, np.array(hi) - vals)],
                fmt="none", ecolor="#666666", elinewidth=1.0, capsize=2.5, label="Wilson 95% CI")

    # Subtle group dividers between: All (idx 0) | DOI (1..2) | Pub type (3..6) | Language (7..11)
    for split_after in [0, 2, 6]:
        y_div = (y[split_after] + y[split_after + 1]) / 2.0
        ax.axhline(y_div, color="#CCCCCC", linestyle="--", linewidth=0.7)

    # Print exact percentage aligned on the right column so it never collides with CI whiskers
    for yi, v in zip(y, vals):
        ax.text(104, yi, f"{v:.0f}%", va="center", ha="left", fontsize=8, color=INK)

    ax.set_yticks(y, names, fontsize=8)
    ax.set_xlim(0, 118)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Publications retrieved by automated search (%)")
    ax.set_title("B  Recall by identifier, type, and language")
    ax.legend(frameon=True, facecolor="white", edgecolor="none", loc="center right", bbox_to_anchor=(0.86, 0.46), fontsize=7.5)

    # ── Panel C: Agreement with dataset records (bars & labels both show %, with k/n below) ──
    ax = axs[2]
    cats = ["all", "Fungi", "Microsporidia", "Protist", "Nematode", "Bacteria", "Virus"]
    names_c = {
        "all": "all groups",
        "Fungi": "fungi",
        "Microsporidia": "microsporidia",
        "Protist": "protists",
        "Nematode": "nematodes",
        "Bacteria": "bacteria",
        "Virus": "viruses",
    }
    w = 0.42
    x = np.arange(len(cats)) * 1.28
    for j, (subset, col, lbl) in enumerate([
        ("full text", BLUE, "full text processed (n = 138 pairs)"),
        ("abstract only", LIGHT, "abstract only processed (n = 248 pairs)"),
    ]):
        g = AG[AG.subset.eq(subset)].set_index("category")
        v = [float(g.loc[c, "recall_graph_genus"]) if c in g.index else np.nan for c in cats]
        npairs = [int(g.loc[c, "her_pairs"]) if c in g.index else 0 for c in cats]
        krec = [int(round(vi * nn / 100.0)) if (nn and not np.isnan(vi)) else 0 for vi, nn in zip(v, npairs)]
        xpos = x + (j - 0.5) * 0.50
        ax.bar(xpos, v, width=w, color=col, label=lbl, edgecolor="white", linewidth=0.8)
        for xi, vi, kk, nn in zip(xpos, v, krec, npairs):
            if nn:
                dx = (-0.05 if j == 0 else 0.05) if (kk in (1, 12) and nn in (1, 12)) else 0.0
                ax.text(xi + dx, vi + 1.5, f"{vi:.0f}%\n({kk}/{nn})",
                        ha="center", va="bottom", fontsize=6.4, color=INK, linespacing=1.0)
    ax.axvline(0.64, color="#CCCCCC", linestyle="--", linewidth=0.7)
    ax.set_xticks(x, [names_c[c] for c in cats], rotation=30, ha="right", fontsize=8)
    ax.set_ylim(0, 124)
    ax.set_yticks([0, 20, 40, 60, 80, 100])
    ax.set_ylabel("Dataset beetle–organism pairs recovered (%)")
    ax.set_title("C  Extraction recall on shared publications")
    ax.legend(frameon=False, fontsize=7.8, loc="upper right")

    fig.tight_layout()
    save(fig, "Figure5_crosscheck")


def fig6():
    T = pd.read_csv(A / "out/trend_rescaled_decade.csv")
    fig, axs = plt.subplots(1, 2, figsize=(12.6, 4.4))

    # ── Panel A: Raw publication counts per period (n = 280; pre-1950 pooled) ──
    ax = axs[0]
    x = np.arange(len(T))
    cols = [BLUE] * (len(T) - 1) + [LIGHT]
    ax.bar(x, T.her_publications, color=cols, width=0.68)
    for xi, v in zip(x, T.her_publications):
        ax.text(xi, v + 1.2, str(int(v)), ha="center", va="bottom", fontsize=8, color=INK)
    ax.set_xticks(x, T.period, rotation=25, ha="right", fontsize=8)
    ax.set_ylim(0, 90)
    ax.set_yticks(range(0, 81, 10))
    ax.set_ylabel("Number of source publications (n = 280)")
    ax.set_title("A  Source publications per period (dataset)")

    # ── Panel B: Rescaled per 1,000 works with direct point labels ──
    ax = axs[1]
    for nm, col, mk, lbl, dy in [
        ("openalex", BLUE, "o", "per 1,000 bark and ambrosia beetle works in OpenAlex (n = 22,565)", -4.8),
        ("kg", ORANGE, "s", "per 1,000 works processed by the automated search (n = 11,085)", 3.2),
    ]:
        v = T[f"per1000_{nm}"].values
        lo = T[f"per1000_{nm}_lo"].values
        hi = T[f"per1000_{nm}_hi"].values
        ax.fill_between(x, lo, hi, color=col, alpha=0.14, linewidth=0)
        ax.plot(x, v, color=col, marker=mk, markersize=5, linewidth=2, label=lbl)
        for xi, vi in zip(x, v):
            offset = dy if not (nm == "openalex" and xi in (0, 8)) else (3.0 if xi == 0 else -4.8)
            if nm == "kg" and xi == 0:
                offset = -5.2
            ax.annotate(f"{vi:.1f}", (xi, vi), textcoords="offset points", xytext=(0, offset),
                        ha="center", va="bottom" if offset > 0 else "top", fontsize=7.2, color=col, fontweight="bold")
    ax.set_xticks(x, T.period, rotation=25, ha="right", fontsize=8)
    ax.set_ylabel("Dataset source publications per 1,000 works")
    ax.set_title("B  Rescaled by the size of the literature")
    ax.legend(frameon=False, fontsize=7.8, loc="upper right")
    ax.set_ylim(0, 118)
    fig.tight_layout()
    save(fig, "Figure6_rescaled_trend")

    # Supplementary 5-year version
    T5 = pd.read_csv(A / "out/trend_rescaled_5yr.csv")
    fig, ax = plt.subplots(figsize=(10, 3.8))
    x = np.arange(len(T5))
    for nm, col, mk, lbl in [
        ("openalex", BLUE, "o", "per 1,000 works in OpenAlex"),
        ("kg", ORANGE, "s", "per 1,000 works processed by the automated search"),
    ]:
        ax.fill_between(x, T5[f"per1000_{nm}_lo"], T5[f"per1000_{nm}_hi"], color=col, alpha=0.15, linewidth=0)
        ax.plot(x, T5[f"per1000_{nm}"], color=col, marker=mk, markersize=4, linewidth=2, label=lbl)
    ax.set_xticks(x, T5.period, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Dataset source publications per 1,000 works")
    ax.legend(frameon=False, fontsize=8)
    ax.set_ylim(0, None)
    fig.tight_layout()
    save(fig, "FigureS1_rescaled_trend_5yr")

    # Supplementary: automated-search primary sources with a pathogen / parasite record
    fig, ax = plt.subplots(figsize=(10, 3.8))
    x = np.arange(len(T5))
    ax.fill_between(x, T5.kg_pathogen_per1000_kg_lo, T5.kg_pathogen_per1000_kg_hi, color=ORANGE, alpha=0.15, linewidth=0)
    ax.plot(x, T5.kg_pathogen_per1000_kg, color=ORANGE, marker="s", markersize=4, linewidth=2,
            label="Automated-search primary sources with a pathogen / parasite record, per 1,000 works processed")
    ax.fill_between(x, T5.per1000_kg_lo, T5.per1000_kg_hi, color=BLUE, alpha=0.15, linewidth=0)
    ax.plot(x, T5.per1000_kg, color=BLUE, marker="o", markersize=4, linewidth=2,
            label="Dataset source publications, per 1,000 works processed")
    ax.set_xticks(x, T5.period, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Publications per 1,000 works")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_ylim(0, None)
    fig.tight_layout()
    save(fig, "FigureS2_kg_pathogen_literature_share")


def fig7():
    cand_path = A / "out/candidate_records_public.csv"
    if not cand_path.exists():
        cand_path = A / "out/candidate_records_public.csv"
    C = pd.read_csv(cand_path)
    L = C[C.stratum.eq("literature read by the KG") & C.beetle_rank.eq("species")]
    Fo = L[L.focused]

    fig, axs = plt.subplots(1, 2, figsize=(12.8, 4.4), gridspec_kw={"width_ratios": [1.15, 1.1]})

    # ── Panel A: Sequential inclusion filter funnel (showing both sources and records) ──
    ax = axs[0]
    F = pd.read_csv(A / "out/candidate_funnel.csv")
    edge_cols = [c for c in F.columns if c.startswith("edges:")]
    rec_counts = []
    for idx, r in F.iterrows():
        if r["step"] == "R6 focused":
            rec_counts.append(int(len(Fo) + C[C.stratum.eq("literature read by the KG") & ~C.beetle_rank.eq("species") & C.focused].shape[0]))
        else:
            rec_counts.append(int(r[edge_cols].sum()))
    # Note: focused subset has 855 species-level records in candidates_summary.csv
    CS = pd.read_csv(A / "out/candidates_summary.csv").set_index("stratum")
    rec_counts[-1] = int(CS.loc["literature, focused subset", "edges"])

    labels_clean = [
        lbl.replace("in her dataset", "in the dataset").replace("R6 on-topic source, own statement", "R6 focused primary beetle study")
        for lbl in F.label
    ]
    y = np.arange(len(F))[::-1]
    ax.barh(y, F.sources, color=BLUE, height=0.62)
    for yi, s_cnt, r_cnt in zip(y, F.sources, rec_counts):
        ax.text(s_cnt + 35, yi, f"{int(s_cnt):,} ({r_cnt:,} records)", va="center", fontsize=7.8, color=INK)
    ax.set_yticks(y, labels_clean, fontsize=8)
    ax.grid(axis="x")
    ax.grid(axis="y", visible=False)
    ax.set_xlabel("Publications with at least one record")
    ax.set_title("A  From automated-search records to candidates")
    ax.set_xlim(0, F.sources.max() * 1.36)

    # ── Panel B: Focused candidate publications by category (strictly descending) ──
    ax = axs[1]
    cats = ["Fungi", "Nematode", "Bacteria", "Microsporidia", "Protist", "Virus"]
    nice = {
        "Fungi": "fungi",
        "Nematode": "nematodes",
        "Bacteria": "bacteria",
        "Microsporidia": "microsporidia",
        "Protist": "protists",
        "Virus": "viruses",
    }
    newp = [int(Fo[Fo.category.eq(c) & Fo.new_pair].source_id.nunique()) for c in cats]
    totp = [int(Fo[Fo.category.eq(c)].source_id.nunique()) for c in cats]
    oldp = [t - n for t, n in zip(totp, newp)]

    x = np.arange(len(cats))
    ax.bar(x, newp, color=BLUE, width=0.68, label="adds a beetle–organism pair absent from the dataset")
    ax.bar(x, oldp, bottom=newp, color=LIGHT, width=0.68, label="all pairs already in the dataset (new source)",
           edgecolor="white", linewidth=0.8)
    for xi, tot, nw in zip(x, totp, newp):
        ax.text(xi, tot + 3.5, f"{tot}\n({nw} new)", ha="center", va="bottom", fontsize=7.5, color=INK, linespacing=1.05)
    ax.set_xticks(x, [nice[c] for c in cats], rotation=20, ha="right", fontsize=8)
    ax.set_ylim(0, 325)
    ax.set_ylabel("Focused candidate publications (n = 373)")
    ax.set_title("B  Focused candidate publications by category")
    ax.legend(frameon=False, fontsize=7.8, loc="upper right")

    fig.tight_layout()
    save(fig, "Figure7_candidates")


if __name__ == "__main__":
    fig5()
    fig6()
    fig7()
