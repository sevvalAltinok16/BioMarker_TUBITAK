"""Fizikokimyasal hesaplamalar (BioMarker, analysis paketi).

Bu modül arayüzden bağımsızdır: dizi alır, sözlük/liste döndürür.
Hatalı girdide :class:`SequenceError` fırlatır; mesajlar arayüzde olduğu gibi
gösterilebilsin diye tek cümledir ve kullanıcıya ne yapacağını söyler.

Hesap kaynakları (proje önerisiyle uyumlu):
  - Moleküler ağırlık, GRAVY, stabilite indeksi: Biopython ProteinAnalysis
  - İzoelektrik nokta ve titrasyon eğrisi: pKa değerleriyle özgün hesap
  - Aminoasit bileşimi, hidropati profili: özgün hesap

Dizi kuralları (şartname bölüm 11.2):
  - 20 standart aminoasit kabul edilir.
  - B, Z, J, X, U, O uyarı üretir, hesaplara katılmaz, sayısı "excluded" ile bildirilir.
  - FASTA başlığı, boşluk, satır sonu ve sondaki '*' sessizce temizlenir.
  - Diğer her karakter (rakam, sembol, başka harf) SequenceError üretir.
"""
from __future__ import annotations

from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils.ProtParamData import hw, kd

STANDARD = "ACDEFGHIKLMNPQRSTVWY"
AMBIGUOUS = "BZJXUO"

# pKa seti: EMBOSS değerleri. Bir kaynağa bağlı kalmak için sabit tutuldu,
# istenirse pka= parametresiyle başka set verilebilir. Sistein serbest tiyol
# kabul edilir (disülfid bağları hesaba katılmaz).
PKA_EMBOSS = {
    "Nterm": 8.6, "Cterm": 3.6,
    "K": 10.8, "R": 12.5, "H": 6.5,
    "D": 3.9, "E": 4.1, "C": 8.5, "Y": 10.1,
}
_POSITIVE = ("Nterm", "K", "R", "H")
_NEGATIVE = ("Cterm", "D", "E", "C", "Y")


class SequenceError(ValueError):
    """Geçersiz dizi. Mesaj arayüzde olduğu gibi gösterilebilir."""


# --------------------------------------------------------------------------
# Dizi temizleme
# --------------------------------------------------------------------------
def _normalize(seq: str) -> str:
    """Başlık, boşluk ve sondaki '*' temizlenir, karakterler doğrulanır.

    Belirsiz harfler (B, Z, J, X, U, O) dizide KALIR, konumlar bozulmaz.
    """
    if not isinstance(seq, str):
        raise SequenceError("Dizi metin olmalı. Bir protein dizisi girin.")

    lines = [ln for ln in seq.splitlines() if not ln.lstrip().startswith(">")]
    text = "".join("".join(lines).split()).upper()
    if text.endswith("*"):
        text = text[:-1]

    if not text:
        raise SequenceError("Dizi boş. Bir protein dizisi girin.")

    allowed = set(STANDARD) | set(AMBIGUOUS)
    for pos, ch in enumerate(text, start=1):
        if ch not in allowed:
            raise SequenceError(
                f"Dizide geçersiz karakter var: '{ch}' (konum {pos}). "
                "Bu karakteri silip yeniden deneyin."
            )
    if not any(ch in STANDARD for ch in text):
        raise SequenceError(
            "Dizide hesaplanabilir standart aminoasit yok. Dizinizi kontrol edin."
        )
    return text


def clean_sequence(seq: str) -> tuple[str, int]:
    """Diziyi temizler ve doğrular.

    Döndürür: (temiz_dizi, hariç_tutulan_sayısı)
    temiz_dizi yalnızca 20 standart aminoasidi içerir. Belirsiz harfler
    (B, Z, J, X, U, O) diziden çıkarılır ve sayısı ikinci değerde bildirilir.
    """
    text = _normalize(seq)
    excluded = sum(1 for ch in text if ch in AMBIGUOUS)
    return "".join(ch for ch in text if ch in STANDARD), excluded


def _prepare(seq: str) -> str:
    return clean_sequence(seq)[0]


# --------------------------------------------------------------------------
# Biopython tabanlı hesaplar
# --------------------------------------------------------------------------
def molecular_weight_kda(seq: str) -> float:
    """Moleküler ağırlık (kDa)."""
    return ProteinAnalysis(_prepare(seq)).molecular_weight() / 1000.0


def gravy(seq: str) -> float:
    """GRAVY: Kyte-Doolittle değerlerinin ortalaması."""
    return ProteinAnalysis(_prepare(seq)).gravy()


def instability_index(seq: str) -> float:
    """Stabilite indeksi (Guruprasad). 40 altı kararlı, 40 ve üstü kararsız."""
    return ProteinAnalysis(_prepare(seq)).instability_index()


# --------------------------------------------------------------------------
# Aminoasit bileşimi
# --------------------------------------------------------------------------
def aa_composition(seq: str) -> dict[str, float]:
    """20 standart aminoasit için yüzde bileşimi. Toplam 100'dür."""
    clean = _prepare(seq)
    n = len(clean)
    return {aa: clean.count(aa) * 100.0 / n for aa in STANDARD}


# --------------------------------------------------------------------------
# Net yük, izoelektrik nokta, titrasyon
# --------------------------------------------------------------------------
def _group_counts(clean: str) -> dict[str, int]:
    counts = {aa: clean.count(aa) for aa in "KRHDECY"}
    counts["Nterm"] = 1
    counts["Cterm"] = 1
    return counts


def _charge(counts: dict[str, int], ph: float, pka: dict[str, float]) -> float:
    pos = sum(counts[g] / (1.0 + 10 ** (ph - pka[g])) for g in _POSITIVE)
    neg = sum(counts[g] / (1.0 + 10 ** (pka[g] - ph)) for g in _NEGATIVE)
    return pos - neg


def net_charge(seq: str, ph: float, pka: dict[str, float] | None = None) -> float:
    """Verilen pH'daki net yük (Henderson-Hasselbalch)."""
    pka = pka or PKA_EMBOSS
    return _charge(_group_counts(_prepare(seq)), ph, pka)


def isoelectric_point(seq: str, pka: dict[str, float] | None = None) -> float:
    """Net yükün sıfır olduğu pH (ikiye bölme yöntemiyle)."""
    pka = pka or PKA_EMBOSS
    counts = _group_counts(_prepare(seq))
    lo, hi = 0.0, 14.0
    for _ in range(60):
        mid = (lo + hi) / 2.0
        if _charge(counts, mid, pka) > 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2.0


def titration_curve(seq: str, step: float = 0.1,
                    pka: dict[str, float] | None = None) -> dict:
    """pH 0-14 aralığında net yük eğrisi.

    Döndürür: {"pi": float, "ph": [..], "charge": [..]}
    """
    if step <= 0:
        raise ValueError("step pozitif olmalı.")
    pka = pka or PKA_EMBOSS
    counts = _group_counts(_prepare(seq))
    n = int(round(14.0 / step))
    ph = [round(i * step, 10) for i in range(n + 1)]
    charge = [_charge(counts, p, pka) for p in ph]
    return {"pi": isoelectric_point(seq, pka), "ph": ph, "charge": charge}


# --------------------------------------------------------------------------
# Hidropati profili
# --------------------------------------------------------------------------
# hydrophobic_is_positive: True ise yüksek değer hidrofobik demektir (grafikte
# yukarısı hidrofobik). Hopp-Woods bir hidrofiliklik ölçeğidir, işareti terstir;
# arayüz eksen etiketlerini bu bilgiye göre çevirmelidir.
SCALES = {
    "kyte-doolittle": {"values": kd, "label": "Kyte-Doolittle", "hydrophobic_is_positive": True},
    "hopp-woods": {"values": hw, "label": "Hopp-Woods", "hydrophobic_is_positive": False},
}
_ALIASES = {"kd": "kyte-doolittle", "hw": "hopp-woods"}


def _scale_key(scale: str) -> str:
    key = str(scale).strip().lower().replace("_", "-").replace(" ", "-")
    key = _ALIASES.get(key, key)
    if key not in SCALES:
        raise ValueError(
            f"Bilinmeyen ölçek: '{scale}'. Seçenekler: " + ", ".join(SCALES)
        )
    return key


def scale_info(scale: str = "kyte-doolittle") -> dict:
    """Ölçeğin adı ve yönü: {"label": str, "hydrophobic_is_positive": bool}."""
    info = SCALES[_scale_key(scale)]
    return {"label": info["label"],
            "hydrophobic_is_positive": info["hydrophobic_is_positive"]}


def _check_window(window) -> None:
    if (not isinstance(window, int) or isinstance(window, bool)
            or window < 1 or window % 2 == 0):
        raise ValueError("window tek ve pozitif bir tam sayı olmalı.")


def hydropathy_profile(seq: str, scale: str = "kyte-doolittle",
                       window: int = 9) -> list:
    """Kaydırma penceresiyle hidropati profili.

    Her pencere için ölçek değerlerinin ortalaması alınır. Uzunluk
    len(dizi) - window + 1'dir, dizi pencereden kısaysa boş liste döner.
    Belirsiz harfler (B, Z, J, X, U, O) dizide kalır, konumlar bozulmaz;
    ortalama yalnızca bilinen kalıntılar üzerinden alınır. Penceredeki
    hiçbir kalıntı bilinmiyorsa o nokta None olur.
    """
    _check_window(window)
    table = SCALES[_scale_key(scale)]["values"]
    text = _normalize(seq)
    vals = [table.get(ch) for ch in text]
    out = []
    for i in range(len(vals) - window + 1):
        seg = [v for v in vals[i:i + window] if v is not None]
        out.append(sum(seg) / len(seg) if seg else None)
    return out


def profile_positions(length: int, window: int = 9) -> list[int]:
    """Profil noktalarının dizideki merkez konumları (1 tabanlı), grafiğin x ekseni için."""
    _check_window(window)
    half = window // 2
    return [i + half + 1 for i in range(max(0, length - window + 1))]


# --------------------------------------------------------------------------
# Özet
# --------------------------------------------------------------------------
def analyze(seq: str, window: int = 9) -> dict:
    """Özet kartları ve hidropati profili için tek çağrılık sonuç.

    Anahtarlar: length, mw_kda, pi, gravy, instability, window, profile
    Ek anahtar: excluded (hesaba katılmayan B/Z/J/X/U/O sayısı, uyarı satırı için)

    length, girilen dizinin tam uzunluğudur (hariç tutulanlar dahil).
    profile, Kyte-Doolittle ile hesaplanır ve konumları girilen diziyle hizalıdır.
    Diğer değerler yalnızca standart aminoasitlerle hesaplanır.
    """
    _check_window(window)
    text = _normalize(seq)
    clean = "".join(ch for ch in text if ch in STANDARD)
    pa = ProteinAnalysis(clean)
    return {
        "length": len(text),
        "mw_kda": pa.molecular_weight() / 1000.0,
        "pi": isoelectric_point(clean),
        "gravy": pa.gravy(),
        "instability": pa.instability_index(),
        "window": window,
        "profile": hydropathy_profile(text, "kyte-doolittle", window),
        "excluded": len(text) - len(clean),
    }
