"""İkincil yapı tahmini: Chou-Fasman algoritması (BioMarker, analysis paketi).

Arayüzden bağımsızdır. Dizi alır, segment listesi döndürür. Hatalı dizide
:class:`physchem.SequenceError` fırlatır (kurallar physchem ile aynıdır).

Algoritmanın üç adımı çıktıda görünür (segmentin ``stage`` alanı):
  1. Çekirdek bölge tespiti
       Heliks:  6 kalıntılık pencerede en az 4 heliks oluşturucu (Pα >= 1,00)
       Yaprak:  5 kalıntılık pencerede en az 3 yaprak oluşturucu (Pβ >= 1,00)
  2. Genişletme
       Çekirdek iki yönde, ilerideki 4 kalıntının ortalaması 1,00'ın altına
       düşene kadar büyütülür. Çakışan bölgeler birleştirilir.
       Heliks için prolin (P) engelleyicidir: çekirdeğe girmez, genişleme prolinde durur.
       Bölgenin uçlarındaki oluşturucu olmayan kalıntılar (değer < 1,00) kırpılır.
  3. Çakışma çözümü
       Heliks ve yaprak çakışırsa çakışan aralıktaki ortalama Pα ve Pβ karşılaştırılır,
       kazanmayan bölgeden o aralık çıkarılır (eşitlikte heliks kazanır).

Bölge kabul koşulları: heliks için uzunluk >= 6, ortalama Pα >= 1,03 ve Pα > Pβ;
yaprak için uzunluk >= 5, ortalama Pβ >= 1,05 ve Pβ > Pα.

Segmente girmeyen kalıntılar "kıvrım/düzensiz" sayılır. Dönü (turn) tahmini yoktur.
Belirsiz harfler (B, Z, J, X, U, O) yapı oluşturucu sayılmaz (değer 0), konumlar kaymaz.

Parametreler: Chou & Fasman (1978) eğilim değerleri.
"""
from __future__ import annotations

from .physchem import SequenceError, _normalize  # noqa: F401

# (Pα, Pβ) - Chou & Fasman 1978
CHOU_FASMAN = {
    "E": (1.51, 0.37), "M": (1.45, 1.05), "A": (1.42, 0.83), "L": (1.21, 1.30),
    "K": (1.16, 0.74), "F": (1.13, 1.38), "Q": (1.11, 1.10), "W": (1.08, 1.37),
    "I": (1.08, 1.60), "V": (1.06, 1.70), "D": (1.01, 0.54), "H": (1.00, 0.87),
    "R": (0.98, 0.93), "T": (0.83, 1.19), "S": (0.77, 0.75), "C": (0.70, 1.19),
    "Y": (0.69, 1.47), "N": (0.67, 0.89), "P": (0.57, 0.55), "G": (0.57, 0.75),
}

HELIX_WINDOW, HELIX_FORMERS, HELIX_MIN_LEN, HELIX_AVG = 6, 4, 6, 1.03
SHEET_WINDOW, SHEET_FORMERS, SHEET_MIN_LEN, SHEET_AVG = 5, 3, 5, 1.05
FORMER_MIN = 1.00     # oluşturucu sayılma eşiği (eşit dahil)
EXTEND_WINDOW = 4     # genişleme sırasında bakılan kalıntı sayısı
_EPS = 1e-9

STAGE_LABELS = {"core": "Çekirdek", "extended": "Genişletilmiş", "resolved": "Çakışma çözüldü"}
TYPE_LABELS = {"h": "Heliks", "s": "Yaprak"}


def _mean(values) -> float:
    values = list(values)
    return sum(values) / len(values)


def _extend(vals, text, s, e, stop_at_pro):
    n = len(vals)
    while e + 1 < n:
        if stop_at_pro and text[e + 1] == "P":
            break
        if _mean(vals[e + 1:e + 1 + EXTEND_WINDOW]) < FORMER_MIN - _EPS:
            break
        e += 1
    while s > 0:
        if stop_at_pro and text[s - 1] == "P":
            break
        if _mean(vals[max(0, s - EXTEND_WINDOW):s]) < FORMER_MIN - _EPS:
            break
        s -= 1
    return s, e


def _merge(regions):
    regions.sort()
    merged = []
    for s, e in regions:
        if merged and s <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    return merged


def _find_regions(vals, text, window, min_formers, helix):
    regions = []
    for i in range(len(vals) - window + 1):
        if sum(v >= FORMER_MIN - _EPS for v in vals[i:i + window]) < min_formers:
            continue
        if helix and "P" in text[i:i + window]:
            continue
        a, b = _extend(vals, text, i, i + window - 1, helix)
        while a < b and vals[a] < FORMER_MIN - _EPS:   # zayıf uçları kırp
            a += 1
        while b > a and vals[b] < FORMER_MIN - _EPS:
            b -= 1
        regions.append([a, b])
    return _merge(regions)


def _accept(s, e, pa, pb, helix):
    length = e - s + 1
    avg_a, avg_b = _mean(pa[s:e + 1]), _mean(pb[s:e + 1])
    if helix:
        return length >= HELIX_MIN_LEN and avg_a >= HELIX_AVG - _EPS and avg_a > avg_b
    return length >= SHEET_MIN_LEN and avg_b >= SHEET_AVG - _EPS and avg_b > avg_a


def _resolve(helices, sheets, pa, pb):
    """Çakışmaları çözer. helices/sheets: [{'s','e','resolved'}] listeleri (0 tabanlı)."""
    while True:
        pair = None
        for h in helices:
            for e in sheets:
                lo, hi = max(h["s"], e["s"]), min(h["e"], e["e"])
                if lo <= hi:
                    pair = (h, e, lo, hi)
                    break
            if pair:
                break
        if not pair:
            return helices, sheets
        h, e, lo, hi = pair
        helix_wins = _mean(pa[lo:hi + 1]) >= _mean(pb[lo:hi + 1])
        winner, loser = (h, e) if helix_wins else (e, h)
        loser_list = sheets if helix_wins else helices
        min_len = SHEET_MIN_LEN if helix_wins else HELIX_MIN_LEN
        winner["resolved"] = True
        loser_list.remove(loser)
        for s, t in ((loser["s"], lo - 1), (hi + 1, loser["e"])):
            if t - s + 1 >= min_len:
                loser_list.append({"s": s, "e": t, "resolved": True})


def predict_secondary_structure(seq: str) -> list[dict]:
    """Chou-Fasman ile ikincil yapı segmentleri.

    Döndürür (başlangıca göre sıralı, 1 tabanlı, uçlar dahil):
      [{"start": int, "end": int, "type": "h" | "s", "stage": str, "length": int}, ...]

    type: "h" heliks, "s" yaprak. Segmentte olmayan kalıntılar kıvrım sayılır.
    stage: "core" (yalnızca çekirdek), "extended" (genişletilmiş),
           "resolved" (heliks-yaprak çakışması çözülürken kesildi veya kazandı).
    """
    text = _normalize(seq)
    pa = [CHOU_FASMAN.get(c, (0.0, 0.0))[0] for c in text]
    pb = [CHOU_FASMAN.get(c, (0.0, 0.0))[1] for c in text]

    helices = [{"s": s, "e": e, "resolved": False}
               for s, e in _find_regions(pa, text, HELIX_WINDOW, HELIX_FORMERS, True)
               if _accept(s, e, pa, pb, True)]
    sheets = [{"s": s, "e": e, "resolved": False}
              for s, e in _find_regions(pb, text, SHEET_WINDOW, SHEET_FORMERS, False)
              if _accept(s, e, pa, pb, False)]
    helices, sheets = _resolve(helices, sheets, pa, pb)

    out = []
    for kind, regs, window in (("h", helices, HELIX_WINDOW), ("s", sheets, SHEET_WINDOW)):
        for r in regs:
            length = r["e"] - r["s"] + 1
            if r["resolved"]:
                stage = "resolved"
            else:
                stage = "core" if length == window else "extended"
            out.append({"start": r["s"] + 1, "end": r["e"] + 1, "type": kind,
                        "stage": stage, "length": length})
    out.sort(key=lambda d: d["start"])
    return out


def secondary_structure_string(seq: str) -> str:
    """Kalıntı başına tek harf: H (heliks), E (yaprak), C (kıvrım). Bant çizimi için."""
    n = len(_normalize(seq))
    states = ["C"] * n
    for seg in predict_secondary_structure(seq):
        mark = "H" if seg["type"] == "h" else "E"
        for i in range(seg["start"] - 1, seg["end"]):
            states[i] = mark
    return "".join(states)


def structure_summary(seq: str) -> dict[str, float]:
    """Yüzdeler: {"helix": .., "sheet": .., "coil": ..}. Toplam 100'dür."""
    s = secondary_structure_string(seq)
    n = len(s)
    return {"helix": s.count("H") * 100.0 / n,
            "sheet": s.count("E") * 100.0 / n,
            "coil": s.count("C") * 100.0 / n}
