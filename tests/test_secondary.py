"""secondary (Chou-Fasman) testleri. Çalıştırma: python -m pytest tests -q"""
import pathlib
import random
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from biomarker.analysis import physchem as pc  # noqa: E402
from biomarker.analysis import secondary as sec  # noqa: E402

INS = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    "RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)
STAGES = {"core", "extended", "resolved"}


def seg(start, end, type_, stage, length=None):
    return {"start": start, "end": end, "type": type_, "stage": stage,
            "length": length or end - start + 1}


# ---------- basit ve bilinen durumlar ----------
def test_poli_alanin_tek_heliks():
    assert sec.predict_secondary_structure("A" * 20) == [seg(1, 20, "h", "extended")]


def test_poli_valin_tek_yaprak():
    assert sec.predict_secondary_structure("V" * 15) == [seg(1, 15, "s", "extended")]


@pytest.mark.parametrize("s", ["G" * 30, "P" * 30, "N" * 12])
def test_yapi_olusturmayan_diziler_bos(s):
    assert sec.predict_secondary_structure(s) == []


def test_tam_pencere_uzunlugunda_cekirdek():
    assert sec.predict_secondary_structure("AAAAAA") == [seg(1, 6, "h", "core")]
    assert sec.predict_secondary_structure("VVVVV") == [seg(1, 5, "s", "core")]


def test_pencereden_kisa_dizi_bos():
    assert sec.predict_secondary_structure("AAAAA") == []
    assert sec.predict_secondary_structure("VVVV") == []


def test_prolin_heliksi_boler():
    got = sec.predict_secondary_structure("A" * 10 + "P" + "A" * 10)
    assert got == [seg(1, 10, "h", "extended"), seg(12, 21, "h", "extended")]


def test_zayif_uclar_genislemez():
    # Glisin bloğu heliksin dışında kalmalı
    got = sec.predict_secondary_structure("G" * 10 + "A" * 10 + "G" * 10)
    assert got == [seg(11, 20, "h", "extended")]


# ---------- çakışma çözümü (iç fonksiyon) ----------
def _regions(h, s):
    return ([{"s": a, "e": b, "resolved": False} for a, b in h],
            [{"s": a, "e": b, "resolved": False} for a, b in s])


def test_cakisma_heliks_kazanir_yaprak_kesilir():
    n = 15
    pa = [1.2] * n
    pb = [0.8] * 5 + [0.8] * 5 + [1.3] * 5
    H, S = _regions([(0, 9)], [(5, 14)])
    H, S = sec._resolve(H, S, pa, pb)
    assert [(r["s"], r["e"]) for r in H] == [(0, 9)] and H[0]["resolved"]
    assert [(r["s"], r["e"]) for r in S] == [(10, 14)] and S[0]["resolved"]


def test_cakisma_yaprak_kazanir_kisa_heliks_silinir():
    n = 15
    pa = [1.0] * n
    pb = [1.5] * n
    H, S = _regions([(0, 9)], [(5, 14)])
    H, S = sec._resolve(H, S, pa, pb)
    assert H == []                      # kalan (0-4) 5 kalıntı < 6, silinir
    assert [(r["s"], r["e"]) for r in S] == [(5, 14)] and S[0]["resolved"]


def test_cakisma_ortada_ise_iki_parcaya_boler():
    n = 30
    pa = [0.5] * 10 + [2.0] * 10 + [0.5] * 10      # ortada heliks eğilimi yüksek
    pb = [1.5] * 30
    H, S = _regions([(10, 19)], [(0, 29)])
    H, S = sec._resolve(H, S, pa, pb)
    assert [(r["s"], r["e"]) for r in H] == [(10, 19)]
    assert sorted((r["s"], r["e"]) for r in S) == [(0, 9), (20, 29)]


def test_esitlikte_heliks_kazanir():
    pa = pb = [1.2] * 12
    H, S = _regions([(0, 11)], [(0, 11)])
    H, S = sec._resolve(H, S, pa, pb)
    assert len(H) == 1 and S == []


# ---------- genel özellikler (rastgele diziler) ----------
@pytest.mark.parametrize("seed", range(40))
def test_rastgele_dizilerde_tutarlilik(seed):
    rnd = random.Random(seed)
    n = rnd.randint(1, 250)
    s = "".join(rnd.choice(pc.STANDARD) for _ in range(n))
    segs = sec.predict_secondary_structure(s)
    last_end = 0
    for g in segs:
        assert 1 <= g["start"] <= g["end"] <= n
        assert g["type"] in ("h", "s") and g["stage"] in STAGES
        assert g["length"] == g["end"] - g["start"] + 1
        assert g["length"] >= (6 if g["type"] == "h" else 5)
        assert g["start"] > last_end          # sıralı ve çakışmasız
        last_end = g["end"]
    assert len(sec.secondary_structure_string(s)) == n


def test_proteinde_prolin_heliksin_icinde_yok():
    for g in sec.predict_secondary_structure(INS):
        if g["type"] == "h":
            assert "P" not in INS[g["start"] - 1:g["end"]]


# ---------- insülin (gerçek protein) ----------
def test_insulin_segmentleri_tutarli_ve_deterministik():
    a = sec.predict_secondary_structure(INS)
    assert a == sec.predict_secondary_structure(INS)
    assert a and all(1 <= g["start"] and g["end"] <= 110 for g in a)


def test_insulin_baslangic_heliksi():
    first = sec.predict_secondary_structure(INS)[0]
    assert (first["start"], first["type"]) == (1, "h")


# ---------- yardımcı fonksiyonlar ----------
def test_yapi_dizgesi_ve_ozet():
    s = "A" * 6 + "G" * 6 + "V" * 5
    string = sec.secondary_structure_string(s)
    assert string == "HHHHHH" + "CCCCCC" + "EEEEE"
    summ = sec.structure_summary(s)
    assert summ["helix"] == pytest.approx(6 * 100 / 17)
    assert sum(summ.values()) == pytest.approx(100.0)


def test_ozet_insulin_toplam_yuz():
    assert sum(sec.structure_summary(INS).values()) == pytest.approx(100.0)


def test_belirsiz_harf_konumlari_kaydirmaz():
    s = "G" * 6 + "X" + "A" * 6
    string = sec.secondary_structure_string(s)
    assert string == "CCCCCCC" + "HHHHHH"      # X atılmaz, heliks 8. konumdan başlar


def test_zayif_uclar_kirpilir():
    # Çekirdek penceresine giren zayıf uç kalıntılar segmente dahil olmaz
    assert sec.predict_secondary_structure("G" * 6 + "A" * 6 + "G" * 6) == [
        seg(7, 12, "h", "core")]


# ---------- hata ve etiketler ----------
@pytest.mark.parametrize("bad", ["", "AC1D", None])
def test_gecersiz_dizi_hata_verir(bad):
    with pytest.raises(pc.SequenceError):
        sec.predict_secondary_structure(bad)


def test_fasta_basligi_ve_kucuk_harf_kabul():
    assert sec.predict_secondary_structure(">x\naaaaaa") == [seg(1, 6, "h", "core")]


def test_etiketler():
    assert set(sec.STAGE_LABELS) == STAGES
    assert sec.TYPE_LABELS == {"h": "Heliks", "s": "Yaprak"}


def test_parametre_tablosu_yirmi_aminoasit():
    assert set(sec.CHOU_FASMAN) == set(pc.STANDARD)
