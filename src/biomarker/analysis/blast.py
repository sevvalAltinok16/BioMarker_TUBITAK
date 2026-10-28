"""Homoloji taraması: NCBI BLAST (BioMarker, analysis paketi).

Sorgu Biopython ``NCBIWWW.qblast()`` ile NCBI'nin sunucusuna gönderilir, dönen XML
``NCBIXML`` ile ayrıştırılır (proje önerisindeki yöntem). Çevrimiçi çalışır.

Arayüz için üç önemli nokta:
  - BLAST dakikalar sürebilir ve yüzde ilerleme vermez. ``run_blast`` bekleyen bir
    işlemdir: arayüzde arka plan iş parçacığında (QThread) çağrılmalıdır.
  - ``cancel_event`` (threading.Event) kurulursa bekleme bırakılır ve BlastCancelled
    fırlatılır. ``timeout`` dolarsa BlastTimeout fırlatılır. İki durumda da NCBI'ye
    giden sorgu yarıda kesilemez, arka planda kendiliğinden biter ve sonucu atılır.
  - ``on_tick(geçen_saniye)`` her ~1 sn çağrılır, "Geçen süre 00:42" için kullanılabilir.

Hata türleri (hepsi Türkçe, tek cümle, arayüzde olduğu gibi gösterilebilir):
  BlastConnectionError   internet yok
  BlastServerError       sunucu cevap vermedi / hata döndü  (BlastTimeout da bu sınıftan)
  BlastCancelled         kullanıcı iptal etti
  BlastError             diğer sorunlar (beklenmeyen cevap vb.)
Geçersiz dizi: physchem.SequenceError. Geçersiz parametre: ValueError.
"""
from __future__ import annotations

import io
import re
import socket
import threading
import time
import urllib.error

from .physchem import SequenceError, _normalize

DATABASES = {"swissprot": "Swiss-Prot", "nr": "nr"}
PROGRAMS = ("blastp",)
MIN_QUERY_LENGTH = 10          # kendi seçimimiz: bundan kısa dizi anlamlı sonuç vermez
UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
UNIPROT_URL = "https://www.uniprot.org/uniprotkb/{ac}"
NCBI_PROTEIN_URL = "https://www.ncbi.nlm.nih.gov/protein/{ac}"
MY_API_KEY="SANSÜRLEDİM KUSRA BAKMA"



# --------------------------------------------------------------------------
# Hatalar
# --------------------------------------------------------------------------
class BlastError(ValueError):
    """BLAST sorgusu başarısız oldu. Mesaj arayüzde olduğu gibi gösterilebilir."""


class BlastConnectionError(BlastError):
    pass


class BlastServerError(BlastError):
    pass


class BlastTimeout(BlastServerError):
    pass


class BlastCancelled(BlastError):
    pass


MSG_SERVER = "BLAST sunucusu yanıt vermedi. Yeniden deneyin."
MSG_OFFLINE = "İnternet bağlantısı gerekiyor. Bağlantınızı kontrol edip yeniden deneyin."
MSG_BAD_REPLY = "NCBI'den beklenmeyen bir cevap geldi. Biraz sonra yeniden deneyin."
MSG_FAILED = "NCBI BLAST sorgusu başarısız oldu. Biraz sonra yeniden deneyin."
MSG_CANCEL = "BLAST sorgusu iptal edildi."


# --------------------------------------------------------------------------
# XML ayrıştırma
# --------------------------------------------------------------------------
def _clean_title(hit_def: str) -> tuple[str, str]:
    """(ad, tür) çıkarır. Birleşik başlıklarda yalnızca ilk başlık kullanılır."""
    title = hit_def.split(" >")[0].strip()
    organism = ""
    m = re.search(r"\[([^\]]+)\]\s*$", title)
    if m:
        organism = m.group(1).strip()
        title = title[:m.start()].strip()
    full = re.match(r"RecName:\s*Full=([^;]+)", title)
    if full:
        title = full.group(1).strip()
    return title.rstrip(";").strip(), organism


def _union_length(intervals: list[tuple[int, int]]) -> int:
    total, cur_s, cur_e = 0, None, None
    for s, e in sorted(intervals):
        if cur_e is None or s > cur_e + 1:
            if cur_e is not None:
                total += cur_e - cur_s + 1
            cur_s, cur_e = s, e
        else:
            cur_e = max(cur_e, e)
    if cur_e is not None:
        total += cur_e - cur_s + 1
    return total


def _hit_url(accession: str) -> str:
    acc = accession.split(".")[0]
    if UNIPROT_RE.match(acc):
        return UNIPROT_URL.format(ac=acc)
    return NCBI_PROTEIN_URL.format(ac=accession)


def parse_blast_xml(xml_text: str, evalue: float | None = None,
                    max_results: int | None = None) -> dict:
    """BLAST XML çıktısını sözlüğe çevirir (ağ kullanmaz).

    Döndürür:
      {"program", "database", "version", "query_length", "total_hits",
       "hits": [ {accession, name, organism, description, length,
                  bit_score, score, evalue, identity_percent, coverage_percent,
                  url, alignment: {...}} ]}
    Hit'ler bit skoruna göre azalan sıralıdır. Hizalama, hit'in en iyi HSP'sine aittir.
    evalue verilirse bundan büyük E-değerli hit'ler atılır.
    """
    from Bio.Blast import NCBIXML

    if "<BlastOutput" not in xml_text:
        raise BlastError(MSG_BAD_REPLY)
    try:
        record = NCBIXML.read(io.StringIO(xml_text))
    except Exception:
        raise BlastError(MSG_BAD_REPLY)

    qlen = int(record.query_letters or 0)
    hits = []
    for aln in record.alignments:
        if not aln.hsps:
            continue
        best = max(aln.hsps, key=lambda h: h.bits)
        ev = min(h.expect for h in aln.hsps)
        if evalue is not None and ev > evalue:
            continue
        name, organism = _clean_title(aln.hit_def or "")
        accession = (aln.accession or aln.hit_id or "").strip()
        cover = _union_length([(h.query_start, h.query_end) for h in aln.hsps])
        hits.append({
            "accession": accession,
            "name": name,
            "organism": organism,
            "description": (aln.hit_def or "").split(" >")[0].strip(),
            "length": int(aln.length),
            "bit_score": float(best.bits),
            "score": float(best.score),
            "evalue": float(ev),
            "identity_percent": best.identities * 100.0 / best.align_length,
            "coverage_percent": min(100.0, cover * 100.0 / qlen) if qlen else 0.0,
            "url": _hit_url(accession),
            "alignment": {
                "query_start": int(best.query_start), "query_end": int(best.query_end),
                "subject_start": int(best.sbjct_start), "subject_end": int(best.sbjct_end),
                "query_seq": best.query, "subject_seq": best.sbjct, "midline": best.match,
                "identities": int(best.identities), "positives": int(best.positives or 0),
                "gaps": int(best.gaps or 0), "length": int(best.align_length),
            },
        })
    hits.sort(key=lambda h: (-h["bit_score"], h["evalue"]))
    total = len(hits)
    if max_results is not None:
        hits = hits[:max_results]
    return {
        "program": record.application.lower() if record.application else "blastp",
        "database": record.database or "",
        "version": record.version or "",
        "query_length": qlen,
        "total_hits": total,
        "hits": hits,
    }


def format_alignment(alignment: dict, width: int = 60) -> str:
    """Hizalamayı arayüzde eş aralıklı fontla gösterilecek metne çevirir.

    Sorgu / eşleşme çizgisi / Hedef üç satırı, ``width`` kalıntılık bloklar halinde,
    her satırın başında ve sonunda konum numarasıyla.
    """
    q, mid, s = alignment["query_seq"], alignment["midline"], alignment["subject_seq"]
    qpos, spos = alignment["query_start"], alignment["subject_start"]
    num_w = len(str(max(alignment["query_end"], alignment["subject_end"])))
    lines = []
    for i in range(0, len(q), width):
        qc, mc, sc = q[i:i + width], mid[i:i + width], s[i:i + width]
        q_end = qpos + sum(1 for c in qc if c != "-") - 1
        s_end = spos + sum(1 for c in sc if c != "-") - 1
        lines.append(f"Sorgu  {qpos:>{num_w}}  {qc}  {q_end}")
        lines.append(f"       {'':>{num_w}}  {mc}")
        lines.append(f"Hedef  {spos:>{num_w}}  {sc}  {s_end}")
        lines.append("")
        qpos, spos = q_end + 1, s_end + 1
    return "\n".join(lines).rstrip("\n")


# --------------------------------------------------------------------------
# Sorgu çalıştırma
# --------------------------------------------------------------------------
def _read_handle(handle) -> str:
    data = handle.read()
    return data.decode("utf-8", errors="replace") if isinstance(data, bytes) else str(data)


def _translate_error(exc: Exception) -> BlastError:
    if isinstance(exc, urllib.error.HTTPError):
        return BlastServerError(MSG_SERVER)
    if isinstance(exc, (urllib.error.URLError, socket.timeout, TimeoutError,
                        ConnectionError, OSError)):
        return BlastConnectionError(MSG_OFFLINE)
    return BlastError(MSG_FAILED)


def run_blast(seq: str, database: str = "swissprot", evalue: float = 0.001,
              max_results: int = 50, timeout: float = 1200, email: str | None = None,
              api_key: str | None = None,  # <-- 1. API Key parametresi buraya eklendi
              cancel_event: "threading.Event | None" = None, on_tick=None,
              tick_interval: float = 1.0, qblast=None) -> dict:
    if not api_key:
        api_key = MY_API_KEY
    """NCBI BLAST (blastp) sorgusu çalıştırır ve sonucu :func:`parse_blast_xml` biçiminde döndürür.

    Ek alanlar: "elapsed_seconds", "evalue_threshold".
    ``qblast`` yalnızca testler içindir (sahte fonksiyon verilebilir).
    """
    if database not in DATABASES:
        raise ValueError("Bilinmeyen veritabanı. Seçenekler: " + ", ".join(DATABASES))
    if not (isinstance(evalue, (int, float)) and not isinstance(evalue, bool) and evalue > 0):
        raise ValueError("E-değeri eşiği pozitif bir sayı olmalı.")
    if not isinstance(max_results, int) or isinstance(max_results, bool) or max_results < 1:
        raise ValueError("En çok sonuç sayısı 1 veya daha büyük bir tam sayı olmalı.")
    if timeout <= 0:
        raise ValueError("Zaman aşımı pozitif olmalı.")
    if email is not None and "@" not in email:
        raise ValueError("E-posta adresi geçerli değil.")

    query = _normalize(seq)
    if len(query) < MIN_QUERY_LENGTH:
        raise SequenceError(
            f"Dizi BLAST için çok kısa. En az {MIN_QUERY_LENGTH} aminoasit girin.")

    if qblast is None:
        from Bio.Blast import NCBIWWW
        if email:
            NCBIWWW.email = email
        if api_key:
            NCBIWWW.api_key = api_key  # <-- 2. NCBI'a API Key tanıtılıyor
        qblast = NCBIWWW.qblast

    if cancel_event is not None and cancel_event.is_set():
        raise BlastCancelled(MSG_CANCEL)

    box: dict = {}

    def worker() -> None:
        try:
            kwargs = {
                "expect": evalue,
                "hitlist_size": max_results,
                "format_type": "XML"
            }
            
            # BURAYA BİR BİLGİLENDİRME EKLEYELİM:
            print(f"NCBI'a istek gönderiliyor (Veritabanı: {database}). Lütfen bekleyin...")
            
            handle = qblast("blastp", database, query, **kwargs)
            
            print("NCBI'dan yanıt başarıyla alındı, işleniyor...")
            box["xml"] = _read_handle(handle)
        except Exception as exc:  
            print(f"NCBI HATASI DETAYI: {exc}") 
            box["error"] = exc

    t0 = time.monotonic()
    thread = threading.Thread(target=worker, daemon=True, name="biomarker-blast")
    thread.start()
    last_tick = 0.0
    while thread.is_alive():
        thread.join(0.1)
        elapsed = time.monotonic() - t0
        if cancel_event is not None and cancel_event.is_set():
            raise BlastCancelled(MSG_CANCEL)
        if thread.is_alive() and elapsed > timeout:
            raise BlastTimeout(MSG_SERVER)
        if on_tick is not None and elapsed - last_tick >= tick_interval:
            last_tick = elapsed
            on_tick(elapsed)

    if "error" in box:
        raise _translate_error(box["error"])
    result = parse_blast_xml(box["xml"], evalue=evalue, max_results=max_results)
    result["elapsed_seconds"] = time.monotonic() - t0
    result["evalue_threshold"] = evalue
    return result
