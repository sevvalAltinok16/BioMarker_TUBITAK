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
def is_online() -> bool:
    """
    BURAYA GELECEK: NCBI BLAST sunucusuna ulaşılabiliyor mu?
    Çevrimdışıysa False döndür; arayüz Homoloji sayfasında "İnternet bağlantısı gerekiyor"
    uyarısını gösterir ve düğmeyi kapatır. Hızlı dönmeli (birkaç yüz ms).
    Şimdilik her zaman True.
    """
    return True


def run_blast(seq: str, database: str = "Swiss-Prot", program: str = "blastp",
              evalue: float = 0.001, max_hits: int = 50):
    """
    BURAYA GELECEK: NCBI BLAST isteği. Dakikalar sürebilir.
    Arayüz bu fonksiyonu ARKA PLAN İŞ PARÇACIĞINDA çağırır, yani içinde beklemek
    (bloklamak) sorun değil. İptal arayüz tarafında yapılır.
    Hata durumları:
      - sunucu yanıt vermezse   : raise TimeoutError
      - kullanıcıya gösterilecek tek cümlelik hata: raise AnalysisError("...")
      - diğer her istisna arayüzde genel bir "başarısız oldu" mesajına çevrilir
    Döndürülecek: sonuç listesi (arayüz skora göre sıralar), her sonuç bir sözlük:
        {"accession": "P01315", "species": "Sus scrofa", "score": 214,
         "evalue": 2e-75, "identity": 96, "coverage": 100,
         "alignment": {"query_start": 25, "query_end": 54,
                       "subject_start": 25, "subject_end": 54,
                       "query":   "FVNQHLCGSHLVEALYLVCGERGFFYTPKT",
                       "match":   "|||||||||||||||||||||||||||||| ",
                       "subject": "FVNQHLCGSHLVEALYLVCGERGFFYTPKA"}}
    Eşleşme yoksa boş liste döndür.
    Şimdilik 3 saniye bekleyip ÖRNEK sonuçlar döndürüyor.
    Denemek için (yalnızca geliştirme): BIOMARKER_DEMO_BLAST=timeout ortam değişkeni
    zaman aşımı hatası, =error ise genel hata üretir.
    """
    import os
    import time
    time.sleep(3)  # ÖRNEK: gerçek istek çok daha uzun sürer
    demo = os.environ.get("BIOMARKER_DEMO_BLAST", "")
    if demo == "timeout":
        raise TimeoutError("BLAST sunucusu yanıt vermedi.")
    if demo == "error":
        raise RuntimeError("örnek hata")

    if len(seq) >= 54:
        q, qs = seq[24:54], 25
    else:
        q, qs = seq, 1
    qe = qs + len(q) - 1
    rows = [
        ("P01315", "Sus scrofa", 214, 2e-75, 96, 100),
        ("P01321", "Canis lupus", 214, 2e-75, 96, 100),
        ("P01317", "Bos taurus", 198, 4e-69, 91, 100),
        ("P01325", "Mus musculus", 172, 6e-59, 82, 97),
    ]
    out = []
    for i, (acc, sp, sc, ev, idn, cov) in enumerate(rows):
        subj = list(q)
        for k in range(i + 1):  # ÖRNEK: her satırda biraz farklı hizalama
            pos = len(subj) - 1 - k * 4
            if pos >= 0:
                subj[pos] = "A" if subj[pos] != "A" else "G"
        subj = "".join(subj)
        match = "".join("|" if a == b else " " for a, b in zip(q, subj))
        out.append({
            "accession": acc, "species": sp, "score": sc, "evalue": ev,
            "identity": idn, "coverage": cov,
            "alignment": {"query_start": qs, "query_end": qe, "subject_start": qs,
                          "subject_end": qe, "query": q, "match": match, "subject": subj},
        })
    return out


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


# ---------------------------------------------------------------------------
# PDF rapor  ->  Raporlar sayfası ("PDF rapor al", Ctrl+E)
# ---------------------------------------------------------------------------
def create_pdf_report(info: dict, sections: list, title: str, note: str, path: str):
    """
    BURAYA GELECEK: gerçek PDF raporu (örn. reportlab ile) path konumuna yazılacak.
      info     : sekmenin analiz sonucu. Anahtarlar: accession, entry, label, seq,
                 analysis (analyze sonucu), ss (ikincil yapı), motifs, blast
                 (blast["hits"] listesi, yoksa None), report_hits (Homoloji sayfasında
                 "Rapora ekle" ile işaretlenen erişim numaraları)
      sections : dahil edilecek bölümler, sırayla: "ozet", "fizikokimya", "ikincil",
                 "motif", "homoloji"
      title    : kullanıcının yazdığı başlık
      note     : kullanıcının yazdığı not (boş olabilir)
    Başarısızlıkta tek cümlelik AnalysisError fırlat, arayüz bunu gösterir.
    Şimdilik içinde bir satır yazı olan YER TUTUCU bir PDF yazıyor.
    """
    def ascii_text(t):
        table = str.maketrans("çğıöşüÇĞİÖŞÜ", "cgiosuCGIOSU")
        t = t.translate(table).encode("ascii", "replace").decode("ascii")
        return t.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")

    lines = [ascii_text(title or "BioMarker raporu"), "(yer tutucu rapor)",
             "Bolumler: " + ", ".join(sections)]
    stream = "BT /F1 14 Tf 72 780 Td 18 TL " + " ".join(f"({ln}) Tj T*" for ln in lines) + " ET"
    stream = stream.encode("ascii")
    objs = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Contents 4 0 R "
        b"/Resources << /Font << /F1 5 0 R >> >> >>",
        b"<< /Length %d >>\nstream\n" % len(stream) + stream + b"\nendstream",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
    ]
    out, offsets = b"%PDF-1.4\n", []
    for i, o in enumerate(objs, 1):
        offsets.append(len(out))
        out += b"%d 0 obj\n" % i + o + b"\nendobj\n"
    xref = len(out)
    out += b"xref\n0 %d\n0000000000 65535 f \n" % (len(objs) + 1)
    for off in offsets:
        out += b"%010d 00000 n \n" % off
    out += b"trailer\n<< /Size %d /Root 1 0 R >>\nstartxref\n%d\n%%%%EOF" % (len(objs) + 1, xref)
    with open(path, "wb") as f:
        f.write(out)