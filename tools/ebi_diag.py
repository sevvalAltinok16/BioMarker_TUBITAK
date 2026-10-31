"""EMBL-EBI BLAST tanı aracı: gerçek servisin cevabını ekrana yazar ve XML'i dosyaya kaydeder.

Kullanım (proje kökünden; BioMarker kodundan bağımsızdır, yalnızca Biopython ister):
    python tools/ebi_diag.py --email senin@epostan.com

Ne yapar: insülin dizisini EBI'nin ncbiblast servisine gönderir, durumu birkaç saniyede bir
yazar, bitince sonucu ebi_sonuc.xml dosyasına kaydeder ve XML'in biçimini (BLAST XML mi)
ile kaç hit bulunduğunu söyler. Bir şey ters giderse sunucunun ham cevabını gösterir.
"""
from __future__ import annotations

import argparse
import io
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"
INS = ("MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
       "RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN")


def stamp(msg, t0):
    print(f"[{time.monotonic() - t0:6.1f} sn] {msg}", flush=True)


def call(url, data=None, timeout=60):
    body = urllib.parse.urlencode(data).encode() if data is not None else None
    req = urllib.request.Request(url, data=body, headers={"User-Agent": "BioMarker-diag/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.status, r.read().decode("utf-8", errors="replace").strip()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--email", required=True)
    ap.add_argument("--db", default="uniprotkb_swissprot", help="uniprotkb_swissprot veya uniprotkb")
    ap.add_argument("--minutes", type=float, default=10)
    ap.add_argument("--seq", default=INS)
    ap.add_argument("--out", default="ebi_sonuc.xml")
    a = ap.parse_args()

    t0 = time.monotonic()
    stamp(f"Gönderiliyor: blastp / {a.db} / {len(a.seq)} aa", t0)
    try:
        code, job = call(BASE + "/run", {"email": a.email, "program": "blastp", "stype": "protein",
                                         "database": a.db, "sequence": ">BioMarker\n" + a.seq})
    except urllib.error.HTTPError as e:
        stamp(f"GÖNDERİM HATASI: HTTP {e.code} {e.reason}", t0)
        print(e.read().decode("utf-8", errors="replace")[:1000])
        return 1
    except Exception as e:  # noqa: BLE001
        stamp(f"GÖNDERİM HATASI (bağlantı): {type(e).__name__}: {e}", t0)
        return 1
    stamp(f"HTTP {code}. İş kimliği: {job}", t0)
    if not job.startswith("ncbiblast-"):
        print("Beklenen biçimde bir iş kimliği değil. Sunucunun cevabı:\n", job[:1000])
        return 1

    deadline = t0 + a.minutes * 60
    while time.monotonic() < deadline:
        time.sleep(3)
        try:
            _, status = call(f"{BASE}/status/{job}")
        except Exception as e:  # noqa: BLE001
            stamp(f"Durum sorgusu hatası: {type(e).__name__}: {e}", t0)
            continue
        stamp(f"Durum: {status}", t0)
        if status == "FINISHED":
            break
        if status in ("ERROR", "FAILURE", "NOT_FOUND"):
            return 2
    else:
        stamp("Süre doldu, iş bitmedi.", t0)
        return 3

    _, xml = call(f"{BASE}/result/{job}/xml", timeout=120)
    with open(a.out, "w", encoding="utf-8") as f:
        f.write(xml)
    stamp(f"Sonuç kaydedildi: {a.out} ({len(xml)} karakter)", t0)
    print("İlk 300 karakter:\n" + xml[:300].replace("\n", " "))
    if "<BlastOutput" in xml:
        try:
            from Bio.Blast import NCBIXML
            rec = NCBIXML.read(io.StringIO(xml))
            first = rec.alignments[0] if rec.alignments else None
            stamp(f"BLAST XML biçimi TAMAM. {len(rec.alignments)} hit"
                  + (f", ilk hit: {first.accession} | {first.hit_def[:60]}" if first else ""), t0)
        except Exception as e:  # noqa: BLE001
            stamp(f"BLAST XML ama ayrıştırılamadı: {type(e).__name__}: {e}", t0)
    else:
        stamp("Bu BLAST XML değil (başka bir biçim). Dosyayı bana gönder, ayrıştırıcıyı ona göre uyarlarım.", t0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
