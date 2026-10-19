"""Motif ve alan taraması: PROSITE kalıpları (BioMarker, analysis paketi).

PROSITE kalıpları düzenli ifadelere (regex) çevrilir ve dizi üzerinde taranır.
Arayüzden bağımsızdır, ağ kullanmaz (çevrimdışı çalışır).

Veri kaynağı:
  - ``analysis/data/prosite.dat`` dosyası varsa tüm PATTERN kayıtları oradan okunur.
    Dosya https://ftp.expasy.org/databases/prosite/ adresinden indirilir.
    PROSITE telif haklıdır, lisans koşullarına uyulmalıdır.
  - Dosya yoksa küçük bir yerleşik örnek set kullanılır (:func:`database_info`
    hangisinin kullanıldığını söyler).
  - Dosya var ama okunamıyorsa sessizce örnek sete düşülmez, hata verilir.

PROSITE kalıp sözdizimi:
  x               herhangi bir aminoasit        ->  .
  [ALT]           listedekilerden biri          ->  [ALT]
  {AM}            listedekiler dışında          ->  [^AM]
  e(n), e(n,m)    tekrar                        ->  e{n}, e{n,m}
  <  ve  >        N-ucu ve C-ucu                ->  ^  ve  $   ([G>] -> (?:G|$))
Çakışan eşleşmeler de bulunur.

PROSITE'ın "SKIP-FLAG=TRUE" işaretli kalıpları çok sık rastlandığı için varsayılan
olarak taramaya girmez (ScanProsite ile aynı), include_frequent=True ile açılır.
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass

from .physchem import _normalize

DEFAULT_DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "prosite.dat")
PROSITE_URL = "https://prosite.expasy.org/{ac}"


class MotifDatabaseError(ValueError):
    """PROSITE dosyası okunamadı. Mesaj arayüzde olduğu gibi gösterilebilir."""


@dataclass
class Pattern:
    accession: str
    name: str
    description: str
    prosite_pattern: str
    regex: str
    search: "re.Pattern"
    frequent: bool


# --------------------------------------------------------------------------
# PROSITE -> regex
# --------------------------------------------------------------------------
_ELEMENT = re.compile(r"(x|[A-Z]|\[[A-Z<>]+\]|\{[A-Z]+\})(?:\((\d+)(?:,(\d+))?\))?")


def prosite_to_regex(pattern: str) -> str:
    """PROSITE kalıbını regex'e çevirir. Geçersiz kalıpta ValueError fırlatır."""
    p = pattern.strip().rstrip(".").strip()
    start_anchor = p.startswith("<")
    if start_anchor:
        p = p[1:]
    end_anchor = p.endswith(">")
    if end_anchor:
        p = p[:-1]
    if not p:
        raise ValueError("Boş kalıp.")

    parts = []
    for el in p.split("-"):
        m = _ELEMENT.fullmatch(el.strip())
        if not m:
            raise ValueError(f"Geçersiz kalıp öğesi: {el!r}")
        core, a, b = m.groups()
        if core == "x":
            body = "."
        elif core.startswith("["):
            inner = core[1:-1]
            has_end, has_start = ">" in inner, "<" in inner
            letters = inner.replace(">", "").replace("<", "")
            if not letters:
                raise ValueError(f"Geçersiz kalıp öğesi: {el!r}")
            alts = [f"[{letters}]"]
            if has_end:
                alts.append("$")
            if has_start:
                alts.append("^")
            body = alts[0] if len(alts) == 1 else "(?:" + "|".join(alts) + ")"
        elif core.startswith("{"):
            body = "[^" + core[1:-1] + "]"
        else:
            body = core
        if a is not None:
            body += "{%s}" % a if b is None else "{%s,%s}" % (a, b)
        parts.append(body)
    return ("^" if start_anchor else "") + "".join(parts) + ("$" if end_anchor else "")


def _compile(regex: str) -> "re.Pattern":
    return re.compile("(?=(" + regex + "))")   # çakışan eşleşmeler için


# --------------------------------------------------------------------------
# prosite.dat ayrıştırma
# --------------------------------------------------------------------------
def _make_pattern(ac, name, desc, pa, frequent) -> Pattern:
    regex = prosite_to_regex(pa)
    return Pattern(ac, name, desc, pa, regex, _compile(regex), frequent)


def parse_prosite_dat(text: str) -> tuple[list[Pattern], int]:
    """prosite.dat metninden PATTERN kayıtlarını okur.

    Döndürür: (kalıplar, atlanan_sayısı). MATRIX (profil) ve RULE kayıtları
    yok sayılır. Çevrilemeyen kalıplar atlanır ve sayılır.
    """
    patterns: list[Pattern] = []
    skipped = 0
    ac = name = ""
    is_pattern = False
    de: list[str] = []
    pa: list[str] = []
    frequent = False

    def finish() -> None:
        nonlocal skipped
        if is_pattern and ac and pa:
            try:
                patterns.append(_make_pattern(ac, name, " ".join(de).rstrip(". ").strip(),
                                              "".join(pa), frequent))
            except ValueError:
                skipped += 1

    for line in text.splitlines():
        tag, body = line[:2], line[5:].strip() if len(line) > 5 else ""
        if line.startswith("//"):
            finish()
            ac = name = ""
            is_pattern, de, pa, frequent = False, [], [], False
        elif tag == "ID":
            name, _, kind = body.partition(";")
            name, is_pattern = name.strip(), kind.strip().rstrip(".") == "PATTERN"
        elif tag == "AC":
            ac = body.rstrip(";").strip()
        elif tag == "DE":
            de.append(body)
        elif tag == "PA":
            pa.append(body)
        elif tag == "CC" and "/SKIP-FLAG=TRUE" in body:
            frequent = True
    return patterns, skipped


# --------------------------------------------------------------------------
# Yerleşik örnek set (prosite.dat yoksa). Gerçek PROSITE kalıplarının küçük bir kısmı.
# --------------------------------------------------------------------------
_BUILTIN_DAT = """\
ID   INSULIN; PATTERN.
AC   PS00262;
DE   Insulin family signature.
PA   C-C-{P}-x(2)-C-[STDNEKPI]-x(3)-[LIVMFS]-x(3)-C.
//
ID   ZINC_FINGER_C2H2_1; PATTERN.
AC   PS00028;
DE   Zinc finger C2H2 type domain signature.
PA   C-x(2,4)-C-x(3)-[LIVMFYWC]-x(8)-H-x(3,5)-H.
//
ID   ATP_GTP_A; PATTERN.
AC   PS00017;
DE   ATP/GTP-binding site motif A (P-loop).
PA   [AG]-x(4)-G-K-[ST].
//
ID   EF_HAND_1; PATTERN.
AC   PS00018;
DE   EF-hand calcium-binding domain.
PA   D-{W}-[DNS]-{ILVFYW}-[DENSTG]-[DNQGHRK]-{GP}-[LIVMC]-[DENQSTAGC]-x(2)-[DE]-[LIVMFYW].
//
ID   ASN_GLYCOSYLATION; PATTERN.
AC   PS00001;
DE   N-glycosylation site.
PA   N-{P}-[ST]-{P}.
CC   /SKIP-FLAG=TRUE;
//
ID   CAMP_PHOSPHO_SITE; PATTERN.
AC   PS00004;
DE   cAMP- and cGMP-dependent protein kinase phosphorylation site.
PA   [RK](2)-x-[ST].
CC   /SKIP-FLAG=TRUE;
//
ID   PKC_PHOSPHO_SITE; PATTERN.
AC   PS00005;
DE   Protein kinase C phosphorylation site.
PA   [ST]-x-[RK].
CC   /SKIP-FLAG=TRUE;
//
ID   CK2_PHOSPHO_SITE; PATTERN.
AC   PS00006;
DE   Casein kinase II phosphorylation site.
PA   [ST]-x(2)-[DE].
CC   /SKIP-FLAG=TRUE;
//
ID   MYRISTYL; PATTERN.
AC   PS00008;
DE   N-myristoylation site.
PA   G-{EDRKHPFYW}-x(2)-[STAGCN]-{P}.
CC   /SKIP-FLAG=TRUE;
//
ID   RGD; PATTERN.
AC   PS00016;
DE   Cell attachment sequence.
PA   R-G-D.
CC   /SKIP-FLAG=TRUE;
//
"""

_CACHE: dict = {}


def _resolve(db_path):
    return os.path.abspath(db_path or DEFAULT_DB_PATH)


def _load(db_path=None):
    path = _resolve(db_path)
    if not os.path.isfile(path):
        key = ("builtin",)
        if key not in _CACHE:
            pats, skipped = parse_prosite_dat(_BUILTIN_DAT)
            _CACHE[key] = (pats, skipped, {"source": "builtin", "path": None, "release": None})
        return _CACHE[key]

    key = (path, os.path.getmtime(path), os.path.getsize(path))
    if key not in _CACHE:
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                text = f.read()
        except OSError:
            raise MotifDatabaseError("prosite.dat dosyası okunamadı. Dosya yolunu ve iznini kontrol edin.")
        pats, skipped = parse_prosite_dat(text)
        if not pats:
            raise MotifDatabaseError(
                "prosite.dat dosyasında kalıp bulunamadı. Dosyayı yeniden indirin: "
                "https://ftp.expasy.org/databases/prosite/")
        rel = re.search(r"Release\s+(\d{4}_\d{2})", text[:20000])
        _CACHE[key] = (pats, skipped, {"source": "prosite.dat", "path": path,
                                       "release": rel.group(1) if rel else None})
    return _CACHE[key]


def load_patterns(db_path=None) -> list[Pattern]:
    """Kullanılacak kalıp listesi (prosite.dat veya yerleşik örnek set)."""
    return _load(db_path)[0]


def database_info(db_path=None) -> dict:
    """Veri kaynağı bilgisi, arayüzde uyarı göstermek için.

    {"source": "prosite.dat" | "builtin", "path": str | None, "count": int,
     "skipped": int, "release": str | None}
    source == "builtin" ise yalnızca küçük bir örnek set kullanılıyordur.
    """
    pats, skipped, meta = _load(db_path)
    return {**meta, "count": len(pats), "skipped": skipped}


# --------------------------------------------------------------------------
# Tarama
# --------------------------------------------------------------------------
def scan_motifs(seq: str, include_frequent: bool = False, db_path=None) -> list[dict]:
    """Dizide PROSITE kalıplarını arar.

    Döndürür (başlangıca göre sıralı, konumlar 1 tabanlı, uçlar dahil):
      [{"id": "PS00262", "name": "INSULIN", "pattern": <regex>,
        "start": 95, "end": 109, "match": "CCTSICSLYQLENYC",
        "description": "Insulin family signature",
        "prosite_pattern": <PROSITE gösterimi>, "frequent": bool, "url": str}, ...]

    include_frequent=False iken PROSITE'ın "çok sık rastlanır" (SKIP-FLAG) kalıpları
    taramaya girmez. Eşleşme yoksa boş liste döner.
    Belirsiz harfler (B, Z, J, X, U, O) dizide kalır, konumlar kaymaz.
    """
    text = _normalize(seq)
    results = []
    for p in load_patterns(db_path):
        if p.frequent and not include_frequent:
            continue
        for m in p.search.finditer(text):
            s, e = m.start(1), m.end(1)
            if e <= s:
                continue
            results.append({
                "id": p.accession, "name": p.name, "pattern": p.regex,
                "start": s + 1, "end": e, "match": text[s:e],
                "description": p.description, "prosite_pattern": p.prosite_pattern,
                "frequent": p.frequent, "url": PROSITE_URL.format(ac=p.accession),
            })
    results.sort(key=lambda d: (d["start"], d["id"]))
    return results
