"""BioMarker veri girişi: FASTA metni/dosyası, UniProt ve NCBI kimliği."""
from __future__ import annotations

import re
from dataclasses import dataclass

from .errors import InputError
from .fasta import FastaRecord, parse_fasta, parse_header, read_fasta_file
from .remote import fetch_ncbi, fetch_uniprot

__all__ = [
    "InputError", "LoadedSequence", "FastaRecord", "load_input", "load_file",
    "parse_fasta", "parse_header", "fetch_uniprot", "fetch_ncbi", "classify",
]

UNIPROT_RE = re.compile(
    r"^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?$"
)
NCBI_RE = re.compile(r"^(?:[A-Z]{2}_\d{6,9}|[A-Z]{3}\d{5,7})(?:\.\d+)?$")


@dataclass
class LoadedSequence:
    sequence: str                 # ham dizi (doğrulama physchem'de yapılır)
    header: str                   # FASTA başlığı ('>' olmadan), olmayabilir
    accession: str | None         # erişim no (sekme etiketi için)
    name: str                     # protein adı, olmayabilir
    source: str                   # "text" | "file" | "uniprot" | "ncbi"
    record_count: int             # girdide kaç dizi vardı
    warning: str | None = None    # ör. "3 dizi bulundu, ilki analiz edildi."


def classify(text: str) -> str:
    """Girdinin türünü söyler: 'fasta', 'uniprot', 'ncbi' veya 'sequence'."""
    t = text.strip()
    if t.startswith(">"):
        return "fasta"
    token = t.upper()
    if len(t.split()) == 1:
        if UNIPROT_RE.match(token):
            return "uniprot"
        if NCBI_RE.match(token):
            return "ncbi"
    return "sequence"


def _from_records(records: list[FastaRecord], source: str, accession=None) -> LoadedSequence:
    first = records[0]
    acc, name = parse_header(first.header)
    n = len(records)
    return LoadedSequence(
        sequence=first.sequence,
        header=first.header,
        accession=accession or acc,
        name=name,
        source=source,
        record_count=n,
        warning=f"{n} dizi bulundu, ilki analiz edildi." if n > 1 else None,
    )


def load_input(text: str, *, email: str | None = None, opener=None) -> LoadedSequence:
    """Giriş kutusuna yazılan metni diziye çevirir.

    Metin FASTA, ham dizi, UniProt kimliği veya NCBI kimliği olabilir.
    Hata durumunda InputError fırlatır.
    """
    if not isinstance(text, str) or not text.strip():
        raise InputError("Girdi boş. Bir dizi, FASTA veya kimlik girin.")

    kind = classify(text)
    if kind == "uniprot":
        fasta = fetch_uniprot(text, opener=opener)
        return _from_records(parse_fasta(fasta), "uniprot", accession=text.strip().upper())
    if kind == "ncbi":
        fasta = fetch_ncbi(text, email or "", opener=opener)
        return _from_records(parse_fasta(fasta), "ncbi", accession=text.strip().upper())

    records = parse_fasta(text)
    if not records:
        raise InputError("Girdide dizi bulunamadı. Dizi veya kimlik girin.")
    return _from_records(records, "text")


def load_file(path: str) -> LoadedSequence:
    """FASTA/.txt dosyasından dizi okur. Birden çok dizi varsa ilkini alır."""
    return _from_records(read_fasta_file(path), "file")
