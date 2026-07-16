"""Canonical Bangladesh district name harmonization (GADM / BBS / ERA5)."""

from __future__ import annotations

import re
import unicodedata

# Exact BBS spelling in wage and yield CSVs (curly apostrophe U+2019).
BBS_COXS_BAZAR = "Cox\u2019s bazar"

DISTRICT_MAPPING = {
    "Bandarban": "Banderban",
    "Bogra": "Bogura",
    "Brahamanbaria": "Brahmmanbaria",
    "Chittagong": "Chattogram",
    "Comilla": "Cumilla",
    "Cox'S Bazar": BBS_COXS_BAZAR,
    "Cox's Bazar": BBS_COXS_BAZAR,
    "Jessore": "Jashore",
    "Jhalokati": "Jhalokathi",
    "Khagrachhari": "Khagrachari",
    "Maulvibazar": "Maulavi Bazar",
    "Nawabganj": "Chapai Nawabganj",
    "Netrakona": "Netrokona",
    "Pirojpur": "Perojpur",
}


def normalize_apostrophes(name: str) -> str:
    """Replace Unicode apostrophe variants with ASCII for matching."""
    return (
        name.replace("\u2019", "'")
        .replace("\u2018", "'")
        .replace("\u0060", "'")
    )


def district_key(name) -> str | None:
    """Lowercase join key; strips apostrophes and collapses whitespace."""
    if name is None or (isinstance(name, float) and name != name):
        return None
    s = normalize_apostrophes(str(name).strip().lower())
    s = re.sub(r"['\u2019\u2018]", "", s)
    s = re.sub(r"\s+", " ", s)
    return s


def harmonize_district(name) -> str | None:
    """Map GADM / ERA5 names to BBS district labels."""
    if name is None or (isinstance(name, float) and name != name):
        return None
    s = str(name).strip()
    if s in DISTRICT_MAPPING:
        return DISTRICT_MAPPING[s]
    s_norm = normalize_apostrophes(s)
    for src, tgt in DISTRICT_MAPPING.items():
        if normalize_apostrophes(src) == s_norm:
            return tgt
    # Cox's Bazar: any apostrophe / case variant
    if district_key(s) == "coxs bazar":
        return BBS_COXS_BAZAR
    return s
