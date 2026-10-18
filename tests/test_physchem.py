"""physchem testleri. Çalıştırma: pytest -q (proje kökünden)."""
import math

import pytest
from Bio.SeqUtils.ProtParam import ProteinAnalysis

import pathlib
import sys

# Proje kökündeki src klasörünü otomatik ekler (pytest.ini gerekmez).
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from biomarker.analysis import physchem as pc  # noqa: E402

# İnsan preproinsülin (UniProt P01308, 110 aa)
INS = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    "RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)


# ---------- temizleme ve doğrulama ----------
def test_temizleme_fasta_bosluk_kucuk_harf_yildiz():
    raw = ">sp|P01308|INS_HUMAN\nmalw mrll\n  ppll*\n"
    clean, excluded = pc.clean_sequence(raw)
    assert clean == "MALWMRLLPPLL"
    assert excluded == 0


def test_belirsiz_harfler_haric_tutulur_ve_sayilir():
    clean, excluded = pc.clean_sequence("ACBXZDUOJ")
    assert clean == "ACD"
    assert excluded == 6


@pytest.mark.parametrize("bad", ["ACD1EF", "AC-DE", "AC DE#", "AC@"])
def test_gecersiz_karakter_hata_verir(bad):
    with pytest.raises(pc.SequenceError):
        pc.clean_sequence(bad)


def test_hata_mesaji_karakter_ve_konumu_icerir():
    with pytest.raises(pc.SequenceError) as e:
        pc.clean_sequence("ACD1EF")
    assert "'1'" in str(e.value) and "konum 4" in str(e.value)


@pytest.mark.parametrize("empty", ["", "   \n", ">sadece baslik", "*"])
def test_bos_dizi_hata_verir(empty):
    with pytest.raises(pc.SequenceError):
        pc.clean_sequence(empty)


def test_sadece_belirsiz_harf_hata_verir():
    with pytest.raises(pc.SequenceError):
        pc.clean_sequence("XXXX")


def test_metin_olmayan_girdi_hata_verir():
    with pytest.raises(pc.SequenceError):
        pc.clean_sequence(None)


# ---------- moleküler ağırlık, GRAVY, stabilite ----------
def test_molekuler_agirlik_insulin():
    assert pc.molecular_weight_kda(INS) == pytest.approx(11.98, abs=0.01)


def test_molekuler_agirlik_biopython_ile_ayni():
    assert pc.molecular_weight_kda(INS) == pytest.approx(
        ProteinAnalysis(INS).molecular_weight() / 1000, rel=1e-9)


def test_belirsiz_harfler_hesabi_degistirmez():
    assert pc.molecular_weight_kda(INS + "XBZ") == pytest.approx(
        pc.molecular_weight_kda(INS))


def test_stabilite_indeksi_insulin():
    # Biopython ile hesaplanan beklenen değer: 40,3 (kararsız, eşik 40)
    assert pc.instability_index(INS) == pytest.approx(40.3, abs=0.1)


def test_gravy_insulin():
    assert pc.gravy(INS) == pytest.approx(0.19, abs=0.01)


# ---------- aminoasit bileşimi ----------
def test_bilesim_yirmi_anahtar_ve_toplam_yuz():
    comp = pc.aa_composition(INS)
    assert list(comp) == list(pc.STANDARD)
    assert sum(comp.values()) == pytest.approx(100.0)


def test_bilesim_basit_dizi():
    comp = pc.aa_composition("AAAC")
    assert comp["A"] == pytest.approx(75.0)
    assert comp["C"] == pytest.approx(25.0)
    assert comp["W"] == 0.0


def test_bilesim_losin_sayisi_insulin():
    assert pc.aa_composition(INS)["L"] == pytest.approx(INS.count("L") * 100 / 110)


# ---------- net yük, pI, titrasyon ----------
def test_pi_asidik_ve_bazik_uc_durumlar():
    assert pc.isoelectric_point("D" * 20) < 4.0
    assert pc.isoelectric_point("K" * 20) > 10.0


def test_net_yuk_pi_noktasinda_sifir():
    pi = pc.isoelectric_point(INS)
    assert pc.net_charge(INS, pi) == pytest.approx(0.0, abs=1e-6)


def test_net_yuk_ph_arttikca_azalir():
    charges = [pc.net_charge(INS, ph / 2) for ph in range(0, 29)]
    assert all(a > b for a, b in zip(charges, charges[1:]))


def test_pi_biopython_ile_makul_aralikta():
    # Farklı pKa setleri kullanıldığı için birebir aynı olması beklenmez.
    bio = ProteinAnalysis(INS).isoelectric_point()
    assert abs(pc.isoelectric_point(INS) - bio) < 0.7


def test_pi_ozel_pka_seti_kabul_eder():
    shifted = dict(pc.PKA_EMBOSS, D=2.0, E=2.0)
    assert pc.isoelectric_point(INS, shifted) < pc.isoelectric_point(INS)


def test_titrasyon_egrisi_sekli():
    t = pc.titration_curve(INS)
    assert t["ph"][0] == 0.0 and t["ph"][-1] == pytest.approx(14.0)
    assert len(t["ph"]) == len(t["charge"]) == 141
    assert t["charge"][0] > 0 > t["charge"][-1]
    i = min(range(len(t["ph"])), key=lambda k: abs(t["ph"][k] - t["pi"]))
    assert abs(t["charge"][i]) < 0.3  # pI civarında yük sıfıra yakın


def test_titrasyon_gecersiz_adim():
    with pytest.raises(ValueError):
        pc.titration_curve(INS, step=0)


# ---------- analyze ----------
def test_analyze_anahtarlar_ve_degerler():
    r = pc.analyze(INS)
    for key in ("length", "mw_kda", "pi", "gravy", "instability", "window", "profile"):
        assert key in r
    assert r["length"] == 110
    assert r["window"] == 9
    assert len(r["profile"]) == 110 - 9 + 1
    assert r["excluded"] == 0
    assert r["mw_kda"] == pytest.approx(11.98, abs=0.01)


def test_analyze_haric_tutulan_ve_uzunluk():
    r = pc.analyze(INS + "XX")
    assert r["length"] == 112
    assert r["excluded"] == 2


def test_analyze_kisa_dizide_profil_bos():
    assert pc.analyze("ACDEF", window=9)["profile"] == []


@pytest.mark.parametrize("w", [0, 4, -3, 9.0])
def test_analyze_gecersiz_pencere(w):
    with pytest.raises(ValueError):
        pc.analyze(INS, window=w)


def test_analyze_gecersiz_dizi_hata_verir():
    with pytest.raises(pc.SequenceError):
        pc.analyze("ACD1")


# ---------- hidropati profili ----------
from Bio.SeqUtils.ProtParamData import hw, kd  # noqa: E402


def test_kyte_doolittle_biopython_ile_ayni():
    got = pc.hydropathy_profile(INS, "kyte-doolittle", 9)
    ref = ProteinAnalysis(INS).protein_scale(kd, 9)
    assert len(got) == len(ref) == 110 - 9 + 1
    assert got == pytest.approx(ref)


def test_hopp_woods_biopython_ile_ayni():
    got = pc.hydropathy_profile(INS, "hopp-woods", 9)
    assert got == pytest.approx(ProteinAnalysis(INS).protein_scale(hw, 9))


@pytest.mark.parametrize("name", ["Kyte-Doolittle", "kd", "KD", "kyte_doolittle", "kyte doolittle"])
def test_olcek_adi_varyantlari_kd(name):
    assert pc.hydropathy_profile("ACDEFGHIK", name, 3) == pc.hydropathy_profile(
        "ACDEFGHIK", "kyte-doolittle", 3)


@pytest.mark.parametrize("name", ["Hopp-Woods", "hw", "hopp_woods"])
def test_olcek_adi_varyantlari_hw(name):
    assert pc.hydropathy_profile("ACDEFGHIK", name, 3) == pc.hydropathy_profile(
        "ACDEFGHIK", "hopp-woods", 3)


def test_bilinmeyen_olcek_hata_verir():
    with pytest.raises(ValueError, match="Seçenekler"):
        pc.hydropathy_profile(INS, "yok-boyle-olcek")


def test_iki_olcegin_isareti_ters():
    assert pc.hydropathy_profile("K" * 9, "kd", 9) == [pytest.approx(-3.9)]
    assert pc.hydropathy_profile("K" * 9, "hw", 9) == [pytest.approx(3.0)]


def test_olcek_bilgisi():
    assert pc.scale_info("kd") == {"label": "Kyte-Doolittle", "hydrophobic_is_positive": True}
    assert pc.scale_info("hopp-woods")["hydrophobic_is_positive"] is False
    assert set(pc.SCALES) == {"kyte-doolittle", "hopp-woods"}


def test_pencere_kisa_dizide_bos_liste():
    assert pc.hydropathy_profile("ACDEF", "kd", 9) == []


@pytest.mark.parametrize("w", [0, 2, -1, 9.0, True])
def test_hidropati_gecersiz_pencere(w):
    with pytest.raises(ValueError):
        pc.hydropathy_profile(INS, "kd", w)


def test_hidropati_gecersiz_karakter_hata_verir():
    with pytest.raises(pc.SequenceError):
        pc.hydropathy_profile("ACD1EF", "kd", 3)


def test_belirsiz_harf_konumlari_bozmaz():
    p = pc.hydropathy_profile("ACDXEF", "kd", 3)
    assert len(p) == 4                       # X atılmaz, uzunluk korunur
    assert p[1] == pytest.approx((kd["C"] + kd["D"]) / 2)   # "CDX": X yok sayılır
    assert p[2] == pytest.approx((kd["D"] + kd["E"]) / 2)   # "DXE"


def test_pencerede_hic_bilinen_yoksa_none():
    p = pc.hydropathy_profile("AXXXA", "kd", 3)
    assert p[0] == pytest.approx(kd["A"]) and p[1] is None and p[2] == pytest.approx(kd["A"])


def test_profil_konumlari():
    pos = pc.profile_positions(110, 9)
    assert pos[0] == 5 and pos[-1] == 106 and len(pos) == 102
    assert pc.profile_positions(5, 9) == []
    assert len(pc.profile_positions(110, 9)) == len(pc.hydropathy_profile(INS, "kd", 9))


def test_analyze_profili_hidropati_fonksiyonuyla_ayni():
    assert pc.analyze(INS, window=7)["profile"] == pc.hydropathy_profile(INS, "kd", 7)


def test_analyze_haric_harfte_profil_konumlari_hizali():
    r = pc.analyze(INS + "X")
    assert r["length"] == 111 and len(r["profile"]) == 111 - 9 + 1
