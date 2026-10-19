"""motifs (PROSITE) testleri. Gerçek prosite.dat'a bağlı değildir, kendi örnek dosyasını kurar.
Çalıştırma: python -m pytest tests -q"""
import os
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from biomarker.analysis import motifs as mo  # noqa: E402
from biomarker.analysis import physchem as pc  # noqa: E402

INS = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    "RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)

DAT = """\
CC   Ornek PROSITE dosyasi. Release 2099_01
ID   INSULIN; PATTERN.
AC   PS00262;
DE   Insulin family signature.
PA   C-C-{P}-x(2)-C-[STDNEKPI]-x(3)-[LIVMFS]-x(3)-C.
//
ID   PKC_PHOSPHO_SITE; PATTERN.
AC   PS00005;
DE   Protein kinase C phosphorylation site.
PA   [ST]-x-[RK].
CC   /TAXO-RANGE=??EPV; /SKIP-FLAG=TRUE;
//
ID   CK2_PHOSPHO_SITE; PATTERN.
AC   PS00006;
DE   Casein kinase II phosphorylation site.
PA   [ST]-x(2)-[DE].
CC   /SKIP-FLAG=TRUE;
//
ID   BIR_REPEAT_2; MATRIX.
AC   PS50143;
DE   BIR repeat profile.
MA   /GENERAL_SPEC: ALPHABET='ABCDEFGHIKLMNPQRSTVWXYZ';
//
ID   LONG_ONE; PATTERN.
AC   PS99999;
DE   Cok uzun
DE   ve iki satirlik aciklama.
PA   A-x-
PA   [ST]-G.
//
ID   BROKEN; PATTERN.
AC   PS00000;
DE   Bozuk kalip.
PA   Z-1-?.
//
"""


@pytest.fixture()
def db(tmp_path):
    p = tmp_path / "prosite.dat"
    p.write_text(DAT, encoding="utf-8")
    return str(p)


# ---------- PROSITE -> regex ----------
@pytest.mark.parametrize("pa,rx", [
    ("N-{P}-[ST]-{P}.", "N[^P][ST][^P]"),
    ("[RK](2)-x-[ST]", "[RK]{2}.[ST]"),
    ("C-x(2,4)-C", "C.{2,4}C"),
    ("<M-x-A", "^M.A"),
    ("A-x(2)>", "A.{2}$"),
    ("A-[G>]", "A(?:[G]|$)"),
    ("R-G-D.", "RGD"),
    ("C-C-{P}-x(2)-C-[STDNEKPI]-x(3)-[LIVMFS]-x(3)-C.",
     "CC[^P].{2}C[STDNEKPI].{3}[LIVMFS].{3}C"),
])
def test_prosite_to_regex(pa, rx):
    assert mo.prosite_to_regex(pa) == rx


@pytest.mark.parametrize("bad", ["", ".", "Z-1", "A-(2)", "A--C", "ab-c", "[]-A", "A-[x]"])
def test_gecersiz_kalip_hata_verir(bad):
    with pytest.raises(ValueError):
        mo.prosite_to_regex(bad)


# ---------- prosite.dat ayrıştırma ----------
def test_ayristirma_kayitlar_ve_atlananlar():
    pats, skipped = mo.parse_prosite_dat(DAT)
    ids = [p.accession for p in pats]
    assert ids == ["PS00262", "PS00005", "PS00006", "PS99999"]   # MATRIX yok sayıldı
    assert skipped == 1                                          # BROKEN çevrilemedi


def test_ayristirma_skip_flag_ve_aciklama():
    pats, _ = mo.parse_prosite_dat(DAT)
    by = {p.accession: p for p in pats}
    assert by["PS00005"].frequent and by["PS00006"].frequent
    assert not by["PS00262"].frequent
    assert by["PS00262"].description == "Insulin family signature"
    assert by["PS99999"].description == "Cok uzun ve iki satirlik aciklama"
    assert by["PS99999"].prosite_pattern == "A-x-[ST]-G."        # çok satırlı PA birleşti
    assert by["PS00262"].name == "INSULIN"


# ---------- tarama ----------
def test_insulin_imzasi_bulunur(db):
    got = mo.scan_motifs(INS, db_path=db)
    assert len(got) == 1
    m = got[0]
    assert (m["id"], m["start"], m["end"], m["match"]) == ("PS00262", 95, 109, "CCTSICSLYQLENYC")
    assert m["name"] == "INSULIN" and m["description"] == "Insulin family signature"
    assert m["pattern"] == "CC[^P].{2}C[STDNEKPI].{3}[LIVMFS].{3}C"
    assert m["url"] == "https://prosite.expasy.org/PS00262" and m["frequent"] is False


def test_sik_rastlanan_kaliplar_varsayilan_kapali_ve_acilabilir(db):
    got = mo.scan_motifs(INS, include_frequent=True, db_path=db)
    triples = [(m["id"], m["start"], m["end"]) for m in got]
    assert ("PS00005", 51, 53) in triples and ("PS00005", 54, 56) in triples
    assert ("PS00006", 54, 57) in triples and ("PS00262", 95, 109) in triples
    assert all(m["frequent"] for m in got if m["id"] in ("PS00005", "PS00006"))


def test_sonuclar_baslangica_gore_sirali(db):
    got = mo.scan_motifs(INS, include_frequent=True, db_path=db)
    assert [m["start"] for m in got] == sorted(m["start"] for m in got)


def test_eslesme_yoksa_bos_liste(db):
    assert mo.scan_motifs("GGGGGGGGGG", include_frequent=True, db_path=db) == []


def test_cakisan_eslesmeler(tmp_path):
    p = tmp_path / "x.dat"
    p.write_text("ID   AA; PATTERN.\nAC   PS00001;\nDE   t.\nPA   A-A.\n//\n")
    got = mo.scan_motifs("AAA", db_path=str(p))
    assert [(m["start"], m["end"]) for m in got] == [(1, 2), (2, 3)]


def test_demirlemeler(tmp_path):
    p = tmp_path / "x.dat"
    p.write_text("ID   N1; PATTERN.\nAC   PS00001;\nDE   t.\nPA   <M-x-A.\n//\n"
                 "ID   C1; PATTERN.\nAC   PS00002;\nDE   t.\nPA   K-[LG>].\n//\n")
    n = mo.scan_motifs("MKAMKA", db_path=str(p))
    assert [(m["id"], m["start"], m["end"]) for m in n] == [("PS00001", 1, 3)]  # yalnız N-ucu
    c = mo.scan_motifs("KAKKL", db_path=str(p))
    assert [(m["id"], m["start"], m["end"]) for m in c] == [("PS00002", 4, 5)]   # K + L
    e = mo.scan_motifs("KAKKK", db_path=str(p))
    assert [(m["id"], m["start"], m["end"]) for m in e] == [("PS00002", 5, 5)]   # K + C-ucu


def test_belirsiz_harf_konumlari_bozmaz(tmp_path):
    p = tmp_path / "x.dat"
    p.write_text("ID   AA; PATTERN.\nAC   PS00001;\nDE   t.\nPA   A-x-C.\n//\n")
    got = mo.scan_motifs("GGAXCGG", db_path=str(p))
    assert [(m["start"], m["end"], m["match"]) for m in got] == [(3, 5, "AXC")]


def test_fasta_ve_kucuk_harf_kabul(db):
    assert mo.scan_motifs(">x\n" + INS.lower(), db_path=db)[0]["id"] == "PS00262"


@pytest.mark.parametrize("bad", ["", "ACD1", None])
def test_gecersiz_dizi_hata_verir(db, bad):
    with pytest.raises(pc.SequenceError):
        mo.scan_motifs(bad, db_path=db)


# ---------- veri tabanı bilgisi ve hatalar ----------
def test_database_info_dosya(db):
    info = mo.database_info(db)
    assert info["source"] == "prosite.dat" and info["count"] == 4
    assert info["skipped"] == 1 and info["release"] == "2099_01"
    assert os.path.samefile(info["path"], db)


def test_dosya_yoksa_yerlesik_ornek_set(tmp_path):
    missing = str(tmp_path / "yok.dat")
    info = mo.database_info(missing)
    assert info["source"] == "builtin" and info["path"] is None and info["count"] >= 9
    assert mo.scan_motifs(INS, db_path=missing)[0]["id"] == "PS00262"


def test_yerlesik_set_tum_kaliplari_cevrilebilir(tmp_path):
    info = mo.database_info(str(tmp_path / "yok.dat"))
    assert info["skipped"] == 0


def test_bos_dosya_hata_verir(tmp_path):
    p = tmp_path / "bos.dat"
    p.write_text("bu bir PROSITE dosyasi degil\n")
    with pytest.raises(mo.MotifDatabaseError, match="kalıp bulunamadı"):
        mo.scan_motifs(INS, db_path=str(p))


def test_dosya_degisince_yeniden_okunur(tmp_path):
    p = tmp_path / "x.dat"
    p.write_text("ID   AA; PATTERN.\nAC   PS00001;\nDE   t.\nPA   A-A.\n//\n")
    assert mo.database_info(str(p))["count"] == 1
    p.write_text("ID   AA; PATTERN.\nAC   PS00001;\nDE   t.\nPA   A-A.\n//\n"
                 "ID   BB; PATTERN.\nAC   PS00002;\nDE   t.\nPA   B-B.\n//\n")
    assert mo.database_info(str(p))["count"] == 2


def test_gercek_boyutta_dosya_hizli(tmp_path):
    # ~1500 kalıplık sentetik dosya makul sürede okunmalı ve taranmalı
    import time
    rec = "ID   P{i}; PATTERN.\nAC   PS{i:05d};\nDE   t.\nPA   C-x(2)-[ST]-{{P}}-G.\n//\n"
    p = tmp_path / "big.dat"
    p.write_text("".join(rec.format(i=i) for i in range(1500)))
    t0 = time.time()
    mo.scan_motifs(INS * 20, include_frequent=True, db_path=str(p))
    assert time.time() - t0 < 5
