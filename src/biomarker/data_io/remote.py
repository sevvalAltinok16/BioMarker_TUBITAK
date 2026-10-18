"""UniProt ve NCBI'den dizi getirme.

Standart kütüphane (urllib) kullanılır, ek paket gerekmez. Testlerde ağa çıkmamak
için `opener` parametresiyle sahte bir açıcı verilebilir.
"""
from __future__ import annotations

import urllib.error
import urllib.parse
import urllib.request

from .errors import InputError

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{acc}.fasta"
NCBI_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    "?db=protein&id={acc}&rettype=fasta&retmode=text&tool=BioMarker&email={email}"
)
USER_AGENT = "BioMarker/1.0 (akademik proje)"


def _get(url: str, source: str, ident: str, opener=None, timeout: float = 15) -> str:
    opener = opener or urllib.request.urlopen
    req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
    try:
        with opener(req, timeout=timeout) as resp:
            return resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as e:
        if e.code in (400, 404, 410):
            raise InputError(f"{source} içinde '{ident}' kimliği bulunamadı. Kimliği kontrol edin.")
        raise InputError(f"{source} sunucusu şu an yanıt vermiyor. Biraz sonra yeniden deneyin.")
    except (urllib.error.URLError, OSError):
        raise InputError("İnternet bağlantısı kurulamadı. Bağlantınızı kontrol edip yeniden deneyin.")


def _check_fasta(text: str, source: str, ident: str) -> str:
    if not text.lstrip().startswith(">"):
        raise InputError(f"{source} içinde '{ident}' kimliği bulunamadı. Kimliği kontrol edin.")
    return text


def fetch_uniprot(accession: str, opener=None, timeout: float = 15) -> str:
    """UniProt erişim numarasının FASTA metnini döndürür (ör. P01308, P04637-2)."""
    acc = accession.strip().upper()
    url = UNIPROT_URL.format(acc=urllib.parse.quote(acc))
    return _check_fasta(_get(url, "UniProt", acc, opener, timeout), "UniProt", acc)


def fetch_ncbi(accession: str, email: str, opener=None, timeout: float = 15) -> str:
    """NCBI protein erişim numarasının FASTA metnini döndürür (ör. NP_000198.1)."""
    if not email or "@" not in email:
        raise InputError("NCBI sorguları için e-posta adresi gerekiyor. Ayarlar'dan girin.")
    acc = accession.strip().upper()
    url = NCBI_URL.format(acc=urllib.parse.quote(acc), email=urllib.parse.quote(email))
    return _check_fasta(_get(url, "NCBI", acc, opener, timeout), "NCBI", acc)
