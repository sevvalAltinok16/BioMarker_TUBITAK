"""Hesaplama katmanı için YER TUTUCU.

Arayüz (window.py, ozet_page.py) hesabı buradaki fonksiyonlardan alır.
Hesaplamayı yazan kişi SADECE bu dosyayı değiştirir; arayüze dokunmasına gerek yok.
Aşağıdaki fonksiyonların adları ve döndürdükleri sözlük anahtarları sabit kalmalı.
Şu an hepsi ÖRNEK değer döndürüyor (SAMPLE = True). Gerçek hesap bağlanınca
"sample": False yapılırsa Özet ekranındaki "örnek değerler" uyarısı kendiliğinden kalkar.
"""
import math


class AnalysisError(Exception):
    """Kullanıcıya gösterilecek tek cümlelik hata.
    Örn: raise AnalysisError("Dizi boş. Bir dizi girin.")
    Arayüz mesajı Giriş ekranındaki kırmızı uyarıda gösterir."""


# ---------------------------------------------------------------------------
# Girdi ayrıştırma (arayüz bunu başlık/dizi ayırmak için kullanır)
# ---------------------------------------------------------------------------
def parse_fasta(text: str):
    """İlk kaydı döndürür: (başlık, dizi). Başlıksız ham dizi de kabul edilir."""
    header, parts = "", []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if parts:  # çoklu FASTA: yalnızca ilk kayıt
                break
            header = line[1:].strip()
            continue
        parts.append("".join(line.split()))
    return header, "".join(parts).upper()


def describe(header: str):
    """Başlıktan (erişim no, kayıt adı) çıkarır. Örn. sp|P01308|INS_HUMAN."""
    if not header:
        return "Özel dizi", "Özel dizi"
    if "|" in header:
        bits = header.split("|")
        if len(bits) >= 3:
            return bits[1], bits[2].split()[0]
    first = header.split()[0]
    return first, first


# ---------------------------------------------------------------------------
# 1) Fizikokimyasal analiz  ->  Özet kartları + Hidropati profili
# ---------------------------------------------------------------------------
def analyze(seq: str, window: int = 9) -> dict:
    """
    BURAYA GELECEK (Biopython / ProteinAnalysis):
      - length       : dizi uzunluğu (aa)
      - mw_kda       : moleküler ağırlık, kDa cinsinden
      - pi           : izoelektrik nokta
      - gravy        : ortalama hidropati skoru
      - instability  : stabilite indeksi (40'ın altı "kararlı", üstü "kararsız")
      - profile      : Kyte-Doolittle skorları, `window` pencereli (len(seq)-window+1 değer)
      - standart dışı harf varsa AnalysisError fırlat
    Şimdilik ÖRNEK değer döndürüyor.
    """
    if not seq:
        raise AnalysisError("Dizi boş. Bir dizi girin.")
    window = max(1, min(window, len(seq)))
    n = max(len(seq) - window + 1, 1)
    sample_profile = [2.2 * math.sin(i / 7.0) + 0.8 * math.sin(i / 2.3) for i in range(n)]
    return {
        "sample": True,          # gerçek hesapta False yap
        "length": len(seq),
        "mw_kda": 11.98,         # ÖRNEK
        "pi": 5.39,              # ÖRNEK
        "gravy": 0.12,           # ÖRNEK
        "instability": 36.4,     # ÖRNEK
        "window": window,
        "profile": sample_profile,
    }


# ---------------------------------------------------------------------------
# 2) İkincil yapı (Chou-Fasman)  ->  İkincil yapı sayfası, Özet'teki bant
# ---------------------------------------------------------------------------
def predict_secondary_structure(seq: str):
    """
    BURAYA GELECEK: Chou-Fasman ikincil yapı tahmini.
    Döndürülecek: segment listesi, her biri bir sözlük:
        {"start": 5, "end": 22, "type": "h", "stage": "Genişletilmiş"}
      - start, end : 1 tabanlı, uçlar dahil
      - type       : "h" (heliks) veya "s" (yaprak). Kıvrım ayrıca verilmez,
                     segmentlerin dışında kalan bölgeler kıvrım sayılır.
      - stage      : "Çekirdek", "Genişletilmiş" veya "Çakışma çözüldü"
    Segmentler başlangıca göre sıralı ve birbirini kesmeyen olmalı.
    Şimdilik 110 aa'lık bir dizi için ÖRNEK segmentleri dizi uzunluğuna ölçekleyip döndürüyor.
    """
    n = len(seq)
    scale = n / 110.0
    base = [
        (5, 22, "h", "Genişletilmiş"), (28, 33, "s", "Çekirdek"), (34, 45, "h", "Çekirdek"),
        (50, 57, "s", "Genişletilmiş"), (58, 69, "h", "Çakışma çözüldü"), (79, 86, "s", "Çekirdek"),
    ]
    out, prev = [], 0
    for a, b, t, st in base:
        a2 = max(1, round(a * scale), prev + 1)
        b2 = min(n, round(b * scale))
        if b2 >= a2:
            out.append({"start": a2, "end": b2, "type": t, "stage": st})
            prev = b2
    return out


# ---------------------------------------------------------------------------
# 3) Motif ve alan taraması  ->  Motif sayfası, Özet'teki motif listesi
# ---------------------------------------------------------------------------
def scan_motifs(seq: str):
    """
    BURAYA GELECEK: PROSITE (veya seçilen kaynak) motif/alan taraması.
    Döndürülecek: eşleşme listesi, her eşleşme bir sözlük (aynı motifin birden çok
    eşleşmesi varsa her biri ayrı sözlük olur):
        {"id": "PS00262", "name": "İnsülin ailesi imzası",
         "pattern": "CC[^P].{2}C[STDNEKPI].{3}[LIVMFS].{3}C",
         "start": 95, "end": 109,            # 1 tabanlı, uçlar dahil
         "match": "CCTSICSLYQLENYC",         # eşleşen dizi parçası
         "description": "Kısa Türkçe açıklama."}
    Eşleşme yoksa boş liste döndür.
    Şimdilik ÖRNEK eşleşmeler döndürüyor (dizi yeterince uzunsa).
    """
    d_pkc = "Protein kinaz C tarafından fosforillenebilen serin/treonin içeren bölge."
    base = [
        ("PS00262", "İnsülin ailesi imzası", "CC[^P].{2}C[STDNEKPI].{3}[LIVMFS].{3}C", 95, 109,
         "İnsülin, IGF ve relaksin ailelerinde korunan, disülfid bağlarına katılan sistein kalıntılarını kapsayan imza."),
        ("PS00005", "PKC fosforilasyon bölgesi", "[ST].[RK]", 51, 53, d_pkc),
        ("PS00005", "PKC fosforilasyon bölgesi", "[ST].[RK]", 54, 56, d_pkc),
        ("PS00006", "CK2 fosforilasyon bölgesi", "[ST].{2}[DE]", 54, 57,
         "Kazein kinaz II tarafından fosforillenebilen serin/treonin içeren bölge."),
    ]
    return [
        {"id": i, "name": n, "pattern": pt, "start": a, "end": b,
         "match": seq[a - 1:b], "description": ds}
        for i, n, pt, a, b, ds in base if b <= len(seq)
    ]


# ---------------------------------------------------------------------------
# 4) Homoloji (NCBI BLAST API)  ->  Homoloji sayfası
# ---------------------------------------------------------------------------
# def run_blast(seq, database="Swiss-Prot", program="blastp", evalue=0.001, max_hits=50):
#     """BURAYA GELECEK: NCBI BLAST isteği (uzun sürer, iptal ve zaman aşımı gerekir).
#     Döndürülecek: [{"accession": "P01315", "species": "Sus scrofa", "score": 214,
#                     "evalue": 2e-75, "identity": 96, "coverage": 100,
#                     "query_range": (25, 54), "alignment": (...)}, ...]"""


# ---------------------------------------------------------------------------
# Fizikokimyasal sayfası için ek yer tutucular
# ---------------------------------------------------------------------------
def hydropathy_profile(seq: str, scale: str = "Kyte-Doolittle", window: int = 9):
    """
    BURAYA GELECEK: seçilen ölçeğe ve pencere boyutuna göre hidropati profili.
      scale  : "Kyte-Doolittle" veya "Hopp-Woods"
      window : tek sayı (3, 5, 7, 9 ...)
    Döndürülecek: float listesi, uzunluğu len(seq) - window + 1.
    Şimdilik ÖRNEK dalga döndürüyor.
    """
    n = max(len(seq) - window + 1, 1)
    shift, amp = (0.0, 1.0) if scale == "Kyte-Doolittle" else (1.3, 0.7)
    return [amp * (2.2 * math.sin(i / 7.0 + shift) + 0.8 * math.sin(i / 2.3)) for i in range(n)]


def titration_curve(seq: str):
    """
    BURAYA GELECEK: pH'a göre net yük eğrisi (pH 0-14).
    Döndürülecek: {"pi": float, "ph": [pH değerleri], "charge": [net yük değerleri]}
    Şimdilik ÖRNEK eğri döndürüyor.
    """
    pi = 5.39
    ph = [i * 0.25 for i in range(57)]  # 0 ... 14
    charge = [-9.0 * math.tanh((p - pi) / 3.0) for p in ph]
    return {"pi": pi, "ph": ph, "charge": charge}


def aa_composition(seq: str):
    """
    BURAYA GELECEK: 20 aminoasidin dizideki yüzdesi.
    Döndürülecek: {"A": 8.2, "C": 5.5, ... "Y": 3.6}  (yüzde değerleri)
    Şimdilik ÖRNEK değerler döndürüyor.
    """
    vals = [8.2, 5.5, 1.8, 5.5, 3.6, 8.2, 1.8, 2.7, 2.7, 15.5,
            0.9, 1.8, 4.5, 4.5, 3.6, 4.5, 3.6, 4.5, 0.9, 2.7]
    return dict(zip("ACDEFGHIKLMNPQRSTVWY", vals))