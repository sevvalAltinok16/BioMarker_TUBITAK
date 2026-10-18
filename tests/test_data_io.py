"""data_io testleri (ağa çıkmaz, sahte açıcı kullanır). Çalıştırma: python -m pytest tests -q"""
import pathlib
import sys
import urllib.error

import pytest

# Proje kökündeki src klasörünü otomatik ekler.
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from biomarker import data_io as dio  # noqa: E402
from biomarker.analysis import physchem as pc  # noqa: E402

INS = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    "RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)
UNIPROT_FASTA = (
    ">sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606 GN=INS PE=1 SV=1\n"
    + INS[:60] + "\n" + INS[60:] + "\n"
)
NCBI_FASTA = ">NP_000198.1 insulin preproprotein [Homo sapiens]\n" + INS + "\n"


class FakeResp:
    def __init__(self, text):
        self._data = text.encode("utf-8")

    def read(self):
        return self._data

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


def fake_opener(text):
    def opener(req, timeout=None):
        fake_opener.last_url = req.full_url
        return FakeResp(text)
    return opener


def failing_opener(exc):
    def opener(req, timeout=None):
        raise exc
    return opener


# ---------- sınıflandırma ----------
@pytest.mark.parametrize("text,kind", [
    ("P01308", "uniprot"), ("p01308", "uniprot"), ("Q9Y6K9", "uniprot"),
    ("A0A024R161", "uniprot"), ("P04637-2", "uniprot"),
    ("NP_000198.1", "ncbi"), ("XP_011520136", "ncbi"), ("AAA59172.1", "ncbi"),
    (">x\nACD", "fasta"), ("MALWMRLLPLL", "sequence"), ("ACDEFGHIKL MNPQ", "sequence"),
])
def test_classify(text, kind):
    assert dio.classify(text) == kind


# ---------- FASTA metni ----------
def test_tek_dizi_basliksiz():
    r = dio.load_input("  mal wmr\nllp  ")
    assert r.sequence == "malwmrllp" and r.source == "text" and r.record_count == 1
    assert r.warning is None


def test_cok_dizili_fasta_ilkini_alir_ve_uyarir():
    r = dio.load_input(">a ilk\nACDE\n>b ikinci\nFGHI\n>c\nKLMN")
    assert r.sequence == "ACDE"
    assert r.record_count == 3
    assert r.warning == "3 dizi bulundu, ilki analiz edildi."


def test_bos_kayitlar_atlanir():
    r = dio.load_input(">bos\n>dolu\nACDE")
    assert r.sequence == "ACDE" and r.record_count == 1


def test_windows_satir_sonu_ve_bom():
    r = dio.load_input("\ufeff>x\r\nACD\r\nEFG\r\n")
    assert r.sequence == "ACDEFG"


@pytest.mark.parametrize("bad", ["", "   \n  ", None, ">sadece baslik"])
def test_bos_girdi_hata(bad):
    with pytest.raises(dio.InputError):
        dio.load_input(bad)


# ---------- başlık çözümleme ----------
def test_baslik_uniprot():
    assert dio.parse_header("sp|P01308|INS_HUMAN Insulin OS=Homo sapiens OX=9606") == ("P01308", "Insulin")


def test_baslik_ncbi():
    assert dio.parse_header("NP_000198.1 insulin preproprotein [Homo sapiens]") == (
        "NP_000198.1", "insulin preproprotein")


def test_baslik_genel_ve_bos():
    assert dio.parse_header("seq1 benim proteinim") == ("seq1", "benim proteinim")
    assert dio.parse_header("") == (None, "")


# ---------- UniProt ----------
def test_uniprot_basarili():
    r = dio.load_input("P01308", opener=fake_opener(UNIPROT_FASTA))
    assert r.sequence == INS and r.accession == "P01308"
    assert r.name == "Insulin" and r.source == "uniprot"
    assert fake_opener.last_url == "https://rest.uniprot.org/uniprotkb/P01308.fasta"


def test_uniprot_kucuk_harf_buyutulur():
    dio.load_input("p01308", opener=fake_opener(UNIPROT_FASTA))
    assert fake_opener.last_url.endswith("/P01308.fasta")


def test_uniprot_bulunamadi_404():
    err = urllib.error.HTTPError("u", 404, "nf", {}, None)
    with pytest.raises(dio.InputError, match="bulunamadı"):
        dio.load_input("P99999", opener=failing_opener(err))


def test_uniprot_bos_cevap_bulunamadi():
    with pytest.raises(dio.InputError, match="bulunamadı"):
        dio.load_input("P99999", opener=fake_opener(""))


def test_sunucu_hatasi_500():
    err = urllib.error.HTTPError("u", 503, "x", {}, None)
    with pytest.raises(dio.InputError, match="yanıt vermiyor"):
        dio.load_input("P01308", opener=failing_opener(err))


def test_internet_yok():
    with pytest.raises(dio.InputError, match="İnternet"):
        dio.load_input("P01308", opener=failing_opener(urllib.error.URLError("dns")))


def test_zaman_asimi():
    with pytest.raises(dio.InputError, match="İnternet"):
        dio.load_input("P01308", opener=failing_opener(TimeoutError()))


# ---------- NCBI ----------
def test_ncbi_basarili():
    r = dio.load_input("NP_000198.1", email="a@b.com", opener=fake_opener(NCBI_FASTA))
    assert r.sequence == INS and r.accession == "NP_000198.1" and r.source == "ncbi"
    assert "db=protein" in fake_opener.last_url and "email=a%40b.com" in fake_opener.last_url


@pytest.mark.parametrize("email", [None, "", "mail-yok"])
def test_ncbi_eposta_gerekli(email):
    with pytest.raises(dio.InputError, match="e-posta"):
        dio.load_input("NP_000198.1", email=email, opener=fake_opener(NCBI_FASTA))


def test_ncbi_bulunamadi():
    with pytest.raises(dio.InputError, match="bulunamadı"):
        dio.load_input("NP_000000.1", email="a@b.com", opener=fake_opener("Error: ..."))


# ---------- dosya ----------
def test_dosya_okuma(tmp_path):
    p = tmp_path / "a.fasta"
    p.write_text(">x protein\nACDE\nFGHI\n>y\nKLMN\n", encoding="utf-8")
    r = dio.load_file(str(p))
    assert r.sequence == "ACDEFGHI" and r.source == "file"
    assert r.warning == "2 dizi bulundu, ilki analiz edildi."


def test_dosya_yok():
    with pytest.raises(dio.InputError, match="bulunamadı"):
        dio.load_file("/yok/boyle/dosya.fasta")


def test_dosya_bos(tmp_path):
    p = tmp_path / "bos.fasta"
    p.write_text("")
    with pytest.raises(dio.InputError, match="dizi bulunamadı"):
        dio.load_file(str(p))


def test_pdb_dosyasi_reddedilir(tmp_path):
    p = tmp_path / "x.pdb"
    p.write_text("ATOM")
    with pytest.raises(dio.InputError, match="PDB"):
        dio.load_file(str(p))


def test_ikili_dosya_reddedilir(tmp_path):
    p = tmp_path / "x.bin"
    p.write_bytes(b"PK\x03\x04\x00\x00")
    with pytest.raises(dio.InputError, match="metin dosyası değil"):
        dio.load_file(str(p))


# ---------- analizle birlikte ----------
def test_kimlikten_analize_uctan_uca():
    r = dio.load_input("P01308", opener=fake_opener(UNIPROT_FASTA))
    out = pc.analyze(r.sequence)
    assert out["length"] == 110
    assert out["mw_kda"] == pytest.approx(11.98, abs=0.01)


def test_gecersiz_karakterli_fasta_analizde_hata_verir():
    r = dio.load_input(">x\nACD1EF")
    with pytest.raises(pc.SequenceError):
        pc.analyze(r.sequence)
