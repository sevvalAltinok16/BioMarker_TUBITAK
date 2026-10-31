"""Homoloji taraması: BLAST (EMBL-EBI NCBI BLAST+ servisi) (BioMarker, analysis paketi).

NCBI'nin kendi BLAST sunucusu bu proje için güvenilir çalışmadığından, aynı NCBI BLAST+
programını çalıştıran EMBL-EBI REST servisi kullanılır (iş gönder, durumu sor, sonucu al).
Dönen BLAST XML'i Biopython ``NCBIXML`` ile ayrıştırılır. Çevrimiçi çalışır.

Akış: ``run`` (POST) -> iş kimliği -> ``status`` (GET, birkaç saniyede bir) -> FINISHED ->
``result/<id>/xml`` (GET). Tüm istekler bu modülün kendi kodundadır, bu yüzden zaman aşımı
ve iptal gerçekten çalışır: bekleme döngüsü bırakılır ve başka istek gönderilmez.
(EBI iş iptal uç noktası sunmaz, sunucudaki iş kendiliğinden biter.)

Arayüz için notlar:
  - ``run_blast`` bekleyen bir işlemdir, arka plan iş parçacığında (QThread) çağrılmalıdır.
  - ``cancel_event`` (threading.Event) kurulursa BlastCancelled fırlar.
  - ``on_tick(geçen_saniye)`` her ~1 sn çağrılır ("Geçen süre 00:42" için). Geri çağrıdaki
    hata sorguyu bozmaz.
  - EBI servisi geçerli bir e-posta adresi ister (Ayarlar'dan alınmalı).
  - Veritabanları: "swissprot" (UniProtKB/Swiss-Prot) ve "uniprotkb" (Swiss-Prot + TrEMBL).
    NCBI'nin "nr" veritabanı bu serviste yoktur.
  - En çok sonuç sayısı sunucunun varsayılanıyla (50) sınırlıdır.
  - E-değeri eşiği sonuçlara istemci tarafında uygulanır.

Hata türleri (hepsi Türkçe, tek cümle; ``.detail`` alanında sunucunun ham cevabı bulunur):
  BlastConnectionError   internet yok
  BlastServerError       sunucu cevap vermedi / işi tamamlayamadı  (BlastTimeout da bu sınıftan)
  BlastCancelled         kullanıcı iptal etti
  BlastError             diğer sorunlar (e-posta eksik, beklenmeyen cevap vb.)
Geçersiz dizi: physchem.SequenceError. Geçersiz parametre: ValueError.
"""
from __future__ import annotations

import io
import re
import socket
import time
import urllib.error
import urllib.parse
import urllib.request

from .physchem import SequenceError, _normalize

BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
DATABASES = {
    "swissprot": {"ebi": "uniprotkb_swissprot", "label": "Swiss-Prot"},
    "uniprotkb": {"ebi": "uniprotkb", "label": "UniProtKB"},
}
MIN_QUERY_LENGTH = 10          # kendi seçimimiz: bundan kısa dizi anlamlı sonuç vermez
JOB_ID_RE = re.compile(r"^ncbiblast-\S+$")
UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
UNIPROT_URL = "https://www.uniprot.org/uniprotkb/{ac}"
NCBI_PROTEIN_URL = "https://www.ncbi.nlm.nih.gov/protein/{ac}"
USER_AGENT = "BioMarker/1.0 (akademik proje)"


# --------------------------------------------------------------------------
# Hatalar
# --------------------------------------------------------------------------
class BlastError(ValueError):
    """BLAST sorgusu başarısız oldu. Mesaj arayüzde olduğu gibi gösterilebilir."""

    def __init__(self, message: str, detail: str | None = None):
        super().__init__(message)
        self.detail = detail


class BlastConnectionError(BlastError):
    pass


class BlastServerError(BlastError):
    pass


class BlastTimeout(BlastServerError):
    pass


class BlastCancelled(BlastError):
    pass


MSG_SERVER = "BLAST sunucusu yanıt vermedi. Yeniden deneyin."
MSG_JOB_FAILED = "BLAST sunucusu sorguyu işleyemedi. Dizinizi kontrol edip yeniden deneyin."
MSG_JOB_LOST = "BLAST işi sunucuda bulunamadı. Yeniden deneyin."
MSG_OFFLINE = "İnternet bağlantısı gerekiyor. Bağlantınızı kontrol edip yeniden deneyin."
MSG_BAD_REPLY = "BLAST sunucusundan beklenmeyen bir cevap geldi. Biraz sonra yeniden deneyin."
MSG_REJECTED = "BLAST sunucusu isteği reddetti. E-posta adresini ve diziyi kontrol edin."
MSG_RATE = "Çok fazla istek gönderildi. Biraz bekleyip yeniden deneyin."
MSG_EMAIL = "BLAST için e-posta adresi gerekiyor. Ayarlar'dan girin."
MSG_CANCEL = "BLAST sorgusu iptal edildi."


# --------------------------------------------------------------------------
# XML ayrıştırma
# --------------------------------------------------------------------------
def _clean_title(hit_def: str) -> tuple[str, str]:
    """(ad, tür) çıkarır. UniProt ("Ad OS=Tür OX=..."), NCBI ("Ad [Tür]") ve
    Swiss-Prot ("RecName: Full=Ad; ... [Tür]") biçimlerini tanır."""
    title = (hit_def or "").split(" >")[0].strip()
    organism = ""
    m = re.search(r"\sOS=(.+?)(?=\s+(?:OX|GN|PE|SV)=|$)", title)
    if m:
        organism = m.group(1).strip()
        title = title[:m.start()].strip()
    else:
        m = re.search(r"\[([^\]]+)\]\s*$", title)
        if m:
            organism = m.group(1).strip()
            title = title[:m.start()].strip()
    full = re.match(r"RecName:\s*Full=([^;]+)", title)
    if full:
        title = full.group(1).strip()
    return title.rstrip(";").strip(), organism


def _accession(aln) -> str:
    """Erişim numarası: Hit_accession, yoksa 'sp|P01315|INS_PIG' biçimli kimlikten."""
    for raw in ((aln.accession or ""), (aln.hit_id or "")):
        raw = raw.strip()
        if not raw:
            continue
        parts = raw.split("|")
        if len(parts) >= 2 and parts[0] in ("sp", "tr", "ref", "gb", "emb", "dbj", "pdb"):
            return parts[1].strip()
        return parts[0].strip()
    return ""


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
        raise BlastError(MSG_BAD_REPLY, detail=xml_text[:500])
    try:
        record = NCBIXML.read(io.StringIO(xml_text))
    except Exception as exc:  # noqa: BLE001
        raise BlastError(MSG_BAD_REPLY, detail=f"{type(exc).__name__}: {exc}")

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
        accession = _accession(aln)
        cover = _union_length([(h.query_start, h.query_end) for h in aln.hsps])
        ident = (best.identities * 100.0 / best.align_length) if best.align_length else 0.0
        hits.append({
            "accession": accession,
            "name": name,
            "organism": organism,
            "description": (aln.hit_def or "").split(" >")[0].strip(),
            "length": int(aln.length),
            "bit_score": float(best.bits),
            "score": float(best.score),
            "evalue": float(ev),
            "identity_percent": ident,
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
# EBI REST istemcisi
# --------------------------------------------------------------------------
def _http(url: str, data: dict | None, timeout: float, opener) -> str:
    """Tek bir HTTP isteği. Hataları Türkçe BlastError türlerine çevirir."""
    body = urllib.parse.urlencode(data).encode("utf-8") if data is not None else None
    req = urllib.request.Request(url, data=body, headers={"User-Agent": USER_AGENT})
    try:
        with opener(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace").strip()
    except urllib.error.HTTPError as e:
        try:
            detail = e.read().decode("utf-8", errors="replace")[:500]
        except Exception:  # noqa: BLE001
            detail = None
        detail = f"HTTP {e.code} {e.reason}" + (f": {detail}" if detail else "")
        if e.code in (400, 403, 404, 415, 422):
            raise BlastError(MSG_REJECTED, detail=detail)
        if e.code == 429:
            raise BlastError(MSG_RATE, detail=detail)
        raise BlastServerError(MSG_SERVER, detail=detail)
    except (urllib.error.URLError, socket.timeout, TimeoutError, ConnectionError, OSError) as e:
        raise BlastConnectionError(MSG_OFFLINE, detail=f"{type(e).__name__}: {e}")


def _check_params(database, evalue, max_results, timeout, email, poll_interval):
    if database not in DATABASES:
        raise ValueError("Bilinmeyen veritabanı. Seçenekler: " + ", ".join(DATABASES))
    if not (isinstance(evalue, (int, float)) and not isinstance(evalue, bool) and evalue > 0):
        raise ValueError("E-değeri eşiği pozitif bir sayı olmalı.")
    if not isinstance(max_results, int) or isinstance(max_results, bool) or max_results < 1:
        raise ValueError("En çok sonuç sayısı 1 veya daha büyük bir tam sayı olmalı.")
    if timeout <= 0 or poll_interval <= 0:
        raise ValueError("Zaman aşımı ve sorgulama aralığı pozitif olmalı.")
    if not email or "@" not in email:
        raise BlastError(MSG_EMAIL)


def run_blast(seq: str, database: str = "swissprot", evalue: float = 0.001,
              max_results: int = 50, timeout: float = 300, email: str | None = None,
              cancel_event=None, on_tick=None, tick_interval: float = 1.0,
              poll_interval: float = 3.0, opener=None) -> dict:
    """EMBL-EBI üzerinden blastp sorgusu çalıştırır, sonucu :func:`parse_blast_xml` biçiminde döndürür.

    Ek alanlar: "job_id", "provider", "database_label", "elapsed_seconds", "evalue_threshold".
    ``opener`` yalnızca testler içindir (urllib.request.urlopen yerine sahte açıcı).
    """
    _check_params(database, evalue, max_results, timeout, email, poll_interval)
    query = _normalize(seq)
    if len(query) < MIN_QUERY_LENGTH:
        raise SequenceError(
            f"Dizi BLAST için çok kısa. En az {MIN_QUERY_LENGTH} aminoasit girin.")
    if cancel_event is not None and cancel_event.is_set():
        raise BlastCancelled(MSG_CANCEL)

    opener = opener or urllib.request.urlopen
    http_timeout = min(60.0, max(5.0, float(timeout)))
    t0 = time.monotonic()

    def guard() -> float:
        elapsed = time.monotonic() - t0
        if cancel_event is not None and cancel_event.is_set():
            raise BlastCancelled(MSG_CANCEL)
        if elapsed > timeout:
            raise BlastTimeout(MSG_SERVER, detail=f"{timeout:g} sn içinde sonuç gelmedi.")
        return elapsed

    job_id = _http(BASE_URL + "/run", {
        "email": email, "program": "blastp", "stype": "protein",
        "database": DATABASES[database]["ebi"], "sequence": ">BioMarker\n" + query,
    }, http_timeout, opener)
    if not JOB_ID_RE.match(job_id):
        raise BlastError(MSG_BAD_REPLY, detail=job_id[:500])

    next_poll, last_tick = 0.0, 0.0
    while True:
        elapsed = guard()
        if on_tick is not None and elapsed - last_tick >= tick_interval:
            last_tick = elapsed
            try:
                on_tick(elapsed)
            except Exception:  # noqa: BLE001 - arayüz geri çağrısı sorguyu bozmasın
                pass
        if elapsed >= next_poll:
            status = _http(f"{BASE_URL}/status/{job_id}", None, http_timeout, opener).upper()
            if status == "FINISHED":
                break
            if status in ("ERROR", "FAILURE"):
                raise BlastServerError(MSG_JOB_FAILED, detail=f"Durum: {status}")
            if status == "NOT_FOUND":
                raise BlastServerError(MSG_JOB_LOST, detail=f"Durum: {status}")
            if status not in ("QUEUED", "RUNNING", "PENDING"):
                raise BlastError(MSG_BAD_REPLY, detail=f"Beklenmeyen durum: {status[:200]}")
            next_poll = time.monotonic() - t0 + poll_interval
        time.sleep(0.05)

    guard()
    xml = _http(f"{BASE_URL}/result/{job_id}/xml", None, http_timeout, opener)
    result = parse_blast_xml(xml, evalue=evalue, max_results=max_results)
    result.update({
        "job_id": job_id, "provider": "EMBL-EBI",
        "database_label": DATABASES[database]["label"],
        "elapsed_seconds": time.monotonic() - t0, "evalue_threshold": evalue,
    })
    return result
