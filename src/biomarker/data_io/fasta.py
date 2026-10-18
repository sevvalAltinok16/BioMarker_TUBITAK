"""FASTA metni ve dosyası okuma. Ağ kullanmaz, arayüzden bağımsızdır."""
from __future__ import annotations

import os
import re
from dataclasses import dataclass

from .errors import InputError

MAX_FILE_BYTES = 20 * 1024 * 1024  # 20 MB


@dataclass
class FastaRecord:
    header: str      # '>' işareti olmadan
    sequence: str    # boşluksuz, ham (doğrulama physchem.clean_sequence'te)


def parse_fasta(text: str) -> list[FastaRecord]:
    """Metindeki tüm FASTA kayıtlarını döndürür. Dizisi boş kayıtlar atlanır.

    Başlıksız ham dizi de kabul edilir (header boş olur).
    """
    records: list[FastaRecord] = []
    header: str | None = None
    parts: list[str] = []

    def flush() -> None:
        seq = "".join(parts)
        if seq:
            records.append(FastaRecord(header or "", seq))

    for raw in text.replace("\ufeff", "").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            flush()
            header, parts = line[1:].strip(), []
        else:
            parts.append("".join(line.split()))
    flush()
    return records


def parse_header(header: str) -> tuple[str | None, str]:
    """Başlıktan (erişim no, ad) çıkarır.

    UniProt:  'sp|P01308|INS_HUMAN Insulin OS=Homo sapiens ...' -> ('P01308', 'Insulin')
    NCBI:     'NP_000198.1 insulin preproprotein [Homo sapiens]' -> ('NP_000198.1', 'insulin preproprotein')
    Diğer:    ilk sözcük erişim no, kalanı ad.
    """
    header = header.strip()
    if not header:
        return None, ""
    first, _, rest = header.partition(" ")
    pipes = first.split("|")
    if len(pipes) >= 3 and pipes[0] in ("sp", "tr"):
        accession = pipes[1]
    else:
        accession = first
    name = rest.split(" OS=")[0]
    name = re.sub(r"\s*\[[^\]]*\]\s*$", "", name).strip()
    return accession or None, name


def read_fasta_file(path: str) -> list[FastaRecord]:
    """Dosyadan FASTA kayıtlarını okur."""
    if not os.path.isfile(path):
        raise InputError("Dosya bulunamadı. Dosya yolunu kontrol edin.")
    if path.lower().endswith(".pdb"):
        raise InputError("PDB dosyasından dizi okunamıyor. FASTA dosyası kullanın.")
    if os.path.getsize(path) > MAX_FILE_BYTES:
        raise InputError("Dosya çok büyük. Tek dizili bir FASTA dosyası kullanın.")
    with open(path, "rb") as f:
        data = f.read()
    if b"\x00" in data:
        raise InputError("Bu dosya metin dosyası değil. FASTA veya .txt dosyası kullanın.")
    records = parse_fasta(data.decode("utf-8", errors="replace"))
    if not records:
        raise InputError("Dosyada dizi bulunamadı. Dosyanın içeriğini kontrol edin.")
    return records
