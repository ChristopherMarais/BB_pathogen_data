"""Candidate passage/page indexing (Aho-Corasick + multilingual pathogen/parasite/fungus regex cues) and structured
language-model extraction (Gemma 4 26B-A4B-it via vLLM, temperature 0, JSON schema + quote alignment) for
bark and ambrosia beetle pathogens and parasites.
"""
import argparse
import json
import os
import sys
from pathlib import Path

import pandas as pd
import regex as re
from rapidfuzz import fuzz

DATA_DIR = Path(__file__).resolve().parents[2] / "data"
W = DATA_DIR / "work"
D = DATA_DIR / "dataset"

# =====================================================================
# PART A: Multilingual Pathogen, Parasite, and Fungal Cue Indexing
# =====================================================================
CUES = {
    "pathogen": re.compile(
        r"\b(?:pathogen\w*|entomopath\w*|parasit\w*|parásit\w*|patógen\w*|pathogène\w*|predat\w*|"
        r"natural enem\w*|ennemis naturels|enemigos naturales|Feinde|Gegenspieler|Krankheit\w*|"
        r"Erreger\w*|infect\w*|Beauveria|Metarhizium|Isaria|Cordyceps|Lecanicillium|microsporid\w*|"
        r"Nosema|Unikaryon|Chytridiopsis|Canningia|gregarin\w*|Gregarina|Mattesia|Menzbieria|"
        r"Malamoeba|nematod\w*|Contortylenchus|Parasitorhabditis|Parasitylenchus|Sphaerularioides|"
        r"Ektaphelenchus|Bursaphelenchus|Cryptaphelenchus|Wolbachia|virus\w*|Bacillus|Serratia)\b",
        re.I,
    ),
    "fungus": re.compile(
        r"\b(?:fung\w*|mycang\w*|symbio\w*|yeasts?|Pilz\w*|champignon\w*|hongo\w*|"
        r"myceli\w*|hyph\w*|conidi\w*|ascospor\w*|moulds?|molds?|levadura\w*|levure\w*)\b",
        re.I,
    ),
}


def cues_of(text: str):
    return [k for k, rx in CUES.items() if rx.search(text or "")]


# =====================================================================
# PART B: Language-Model JSON Extraction & Verbatim Quote Verification
# =====================================================================
TYPES = ["pathogen", "fungus"]
RELATIONS = [
    "pathogen", "parasite", "parasitoid", "predator", "phoretic", "endosymbiont",
    "gut_associate", "mycangial_symbiont", "ambrosia_fungus", "associated",
    "isolated_from", "vectored_pathogen", "not_associated", "other", "none",
]

SCHEMA = {
    "type": "object",
    "properties": {
        "claims": {
            "type": "array",
            "items": {
                "type": "object",
                "properties": {
                    "type": {"type": "string", "enum": TYPES},
                    "taxon": {"type": "string"},
                    "partner": {"type": "string"},
                    "partner_full": {"type": "string"},
                    "partner_group": {
                        "type": "string",
                        "enum": [
                            "fungus", "yeast", "bacterium", "virus", "microsporidium",
                            "protist", "nematode", "mite", "insect", "other", "none",
                        ],
                    },
                    "relation": {"type": "string", "enum": RELATIONS},
                    "setting": {"type": "string", "enum": ["natural", "experimental", "unspecified", "none"]},
                    "quote": {"type": "string"},
                },
                "required": [
                    "type", "taxon", "partner", "partner_full",
                    "partner_group", "relation", "setting", "quote",
                ],
            },
        }
    },
    "required": ["claims"],
}

RULES = """Report every pathogen, parasite, or microbial/fungal association claim for the bark or ambrosia beetle(s) named on the page:
- pathogen: an organism that infects, parasitizes, or causes disease/mortality in the beetle (entomopathogenic fungi such as
  Beauveria or Metarhizium, bacteria, viruses, microsporidia, gregarines and other protists, nematodes, mites,
  parasitoid wasps, predators, endosymbionts, gut bacteria). partner = the organism name as written; partner_full = its full
  scientific name if abbreviated on the page, else ""; partner_group; relation one of: pathogen, parasite, parasitoid,
  predator, phoretic, endosymbiont, gut_associate, other; setting natural (found in or on wild or reared beetles) or
  experimental (infection or biological-control bioassay/trial).
- fungus: a fungus or yeast living with the beetle as a nutritional/mycangial symbiont or carried by it to plants without
  reported harm to the beetle.
The partner must be a named organism (scientific or recognized common name), not a generic word such as "fungi", "bacteria", or "nematodes".
For every claim give `quote`: the shortest exact sentence or phrase of the text that states it, copied verbatim (at most 300 characters).
If the text holds no such claim, return an empty list."""


def norm(s: str) -> str:
    s = re.sub(r"(?:\s*\.){4,}", " ... ", str(s or ""))
    return re.sub(r"\s+", " ", s).strip()


def align(quote: str, text: str):
    """Deterministic Gate (i): Verbatim quote alignment against the source passage/page."""
    q, t = norm(quote), norm(text)
    if not q or not t:
        return 0.0, ""
    if q in t:
        return 100.0, q
    if len(q) < 12:
        m = re.search(r"(?<!\w)" + re.escape(q) + r"(?!\w)", t, re.I)
        return (100.0, m.group(0)) if m else (0.0, "")
    al = fuzz.partial_ratio_alignment(q.lower(), t.lower(), score_cutoff=90)
    if al is None:
        return 0.0, ""
    return float(al.score), t[al.dest_start:al.dest_end]
