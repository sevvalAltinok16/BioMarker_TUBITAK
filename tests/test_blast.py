"""blast testleri. Ağa çıkmaz: sahte qblast ve elle kurulmuş BLAST XML'i kullanır.
Çalıştırma: python -m pytest tests -q"""
import io
import pathlib
import sys
import threading
import time
import urllib.error

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from biomarker.analysis import blast as bl  # noqa: E402
from biomarker.analysis import physchem as pc  # noqa: E402

INS = (
    "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKT"
    "RREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN"
)
QSEQ = "FVNQHLCGSHLVEALYLVCGERGFFYTPKT"
HSEQ = "FVNQHLCGSHLVEALYLVCGERGFFYTPKA"
MID = "FVNQHLCGSHLVEALYLVCGERGFFYTPK "


def hsp_xml(n, h):
    ident = sum(1 for a, b in zip(h["qseq"], h["hseq"]) if a == b)
    pos = sum(1 for c in h["mid"] if c != " ")
    gaps = h["qseq"].count("-") + h["hseq"].count("-")
    return f"""<Hsp><Hsp_num>{n}</Hsp_num><Hsp_bit-score>{h['bits']}</Hsp_bit-score>
<Hsp_score>{h['score']}</Hsp_score><Hsp_evalue>{h['evalue']}</Hsp_evalue>
<Hsp_query-from>{h['qfrom']}</Hsp_query-from><Hsp_query-to>{h['qto']}</Hsp_query-to>
<Hsp_hit-from>{h['hfrom']}</Hsp_hit-from><Hsp_hit-to>{h['hto']}</Hsp_hit-to>
<Hsp_query-frame>0</Hsp_query-frame><Hsp_hit-frame>0</Hsp_hit-frame>
<Hsp_identity>{ident}</Hsp_identity><Hsp_positive>{pos}</Hsp_positive><Hsp_gaps>{gaps}</Hsp_gaps>
<Hsp_align-len>{len(h['qseq'])}</Hsp_align-len><Hsp_qseq>{h['qseq']}</Hsp_qseq>
<Hsp_hseq>{h['hseq']}</Hsp_hseq><Hsp_midline>{h['mid']}</Hsp_midline></Hsp>"""


def make_xml(hits, qlen=110, db="swissprot"):
    out = []
    for i, ht in enumerate(hits, 1):
        hsps = "".join(hsp_xml(j, h) for j, h in enumerate(ht["hsps"], 1))
        out.append(f"""<Hit><Hit_num>{i}</Hit_num><Hit_id>{ht['id']}</Hit_id>
<Hit_def>{ht['def']}</Hit_def><Hit_accession>{ht['acc']}</Hit_accession>
<Hit_len>{ht['len']}</Hit_len><Hit_hsps>{hsps}</Hit_hsps></Hit>""")
    return f"""<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput><BlastOutput_program>blastp</BlastOutput_program>
<BlastOutput_version>BLASTP 2.15.0+</BlastOutput_version>
<BlastOutput_reference>ref</BlastOutput_reference><BlastOutput_db>{db}</BlastOutput_db>
<BlastOutput_query-ID>Query_1</BlastOutput_query-ID><BlastOutput_query-def>q</BlastOutput_query-def>
<BlastOutput_query-len>{qlen}</BlastOutput_query-len>
<BlastOutput_param><Parameters><Parameters_matrix>BLOSUM62</Parameters_matrix>
<Parameters_expect>0.001</Parameters_expect><Parameters_gap-open>11</Parameters_gap-open>
<Parameters_gap-extend>1</Parameters_gap-extend><Parameters_filter>F</Parameters_filter></Parameters></BlastOutput_param>
<BlastOutput_iterations><Iteration><Iteration_iter-num>1</Iteration_iter-num>
<Iteration_query-ID>Query_1</Iteration_query-ID><Iteration_query-def>q</Iteration_query-def>
<Iteration_query-len>{qlen}</Iteration_query-len><Iteration_hits>{''.join(out)}</Iteration_hits>
<Iteration_stat><Statistics><Statistics_db-num>1</Statistics_db-num><Statistics_db-len>1</Statistics_db-len>
<Statistics_hsp-len>0</Statistics_hsp-len><Statistics_eff-space>0</Statistics_eff-space>
<Statistics_kappa>0.041</Statistics_kappa><Statistics_lambda>0.267</Statistics_lambda>
<Statistics_entropy>0.14</Statistics_entropy></Statistics></Iteration_stat>
</Iteration></BlastOutput_iterations></BlastOutput>"""


def hsp(bits, evalue, qfrom=25, qto=54, hfrom=25, hto=54, qseq=QSEQ, hseq=HSEQ, mid=MID):
    return dict(bits=bits, score=int(bits * 2.5), evalue=evalue, qfrom=qfrom, qto=qto,
                hfrom=hfrom, hto=hto, qseq=qseq, hseq=hseq, mid=mid)


PIG = dict(id="sp|P01315.2|INS_PIG", acc="P01315", len=108, hsps=[hsp(61.2, "2e-75")],
           **{"def": "RecName: Full=Insulin; Contains: RecName: Full=Insulin B chain; "
                     "Flags: Precursor [Sus scrofa]"})
HUMAN_NP = dict(id="ref|NP_000198.1|", acc="NP_000198", len=110, hsps=[hsp(200.0, "1e-60")],
                **{"def": "insulin preproprotein [Homo sapiens] >gb|AAA59172.1| insulin [Homo sapiens]"})
WEAK = dict(id="sp|Q99999.1|X", acc="Q99999", len=90, hsps=[hsp(20.0, "0.5")],
            **{"def": "RecName: Full=Zayif; [Mus musculus]"})
TWO = dict(id="sp|P11111.1|Y", acc="P11111", len=100, **{"def": "RecName: Full=Iki parca [Bos taurus]"},
           hsps=[hsp(50, "1e-9", 1, 40, 1, 40, "A" * 40, "A" * 40, "A" * 40),
                 hsp(60, "1e-12", 31, 70, 31, 70, "A" * 40, "A" * 40, "A" * 40)])


def fake_qblast(xml, delay=0.0, calls=None):
    def f(program, database, sequence, **kw):
        if calls is not None:
            calls.append(dict(program=program, database=database, sequence=sequence, **kw))
        if delay:
            time.sleep(delay)
        return io.StringIO(xml)
    return f


def raising(exc):
    def f(*a, **k):
        raise exc
    return f


# ---------- XML ayrıştırma ----------
def test_ayristirma_alanlar_ve_siralama():
    r = bl.parse_blast_xml(make_xml([PIG, HUMAN_NP]))
    assert r["program"] == "blastp" and r["database"] == "swissprot"
    assert r["query_length"] == 110 and r["total_hits"] == 2
    # bit skoruna göre azalan: 200 önce
    assert [h["accession"] for h in r["hits"]] == ["NP_000198", "P01315"]


def test_swissprot_hit_alanlari():
    h = bl.parse_blast_xml(make_xml([PIG]))["hits"][0]
    assert (h["accession"], h["name"], h["organism"]) == ("P01315", "Insulin", "Sus scrofa")
    assert h["bit_score"] == pytest.approx(61.2) and h["evalue"] == pytest.approx(2e-75)
    assert h["length"] == 108
    assert h["identity_percent"] == pytest.approx(29 * 100 / 30)
    assert h["coverage_percent"] == pytest.approx(30 * 100 / 110)
    assert h["url"] == "https://www.uniprot.org/uniprotkb/P01315"


def test_nr_hit_birlesik_baslik_ve_ncbi_baglantisi():
    h = bl.parse_blast_xml(make_xml([HUMAN_NP], db="nr"))["hits"][0]
    assert (h["name"], h["organism"]) == ("insulin preproprotein", "Homo sapiens")
    assert h["description"] == "insulin preproprotein [Homo sapiens]"
    assert h["url"] == "https://www.ncbi.nlm.nih.gov/protein/NP_000198"


def test_hizalama_alani():
    a = bl.parse_blast_xml(make_xml([PIG]))["hits"][0]["alignment"]
    assert (a["query_start"], a["query_end"], a["subject_start"], a["subject_end"]) == (25, 54, 25, 54)
    assert a["query_seq"] == QSEQ and a["subject_seq"] == HSEQ and a["midline"] == MID
    assert a["identities"] == 29 and a["positives"] == 29 and a["gaps"] == 0 and a["length"] == 30


def test_birden_fazla_hsp_kapsam_birlesimi_ve_en_iyi_hsp():
    h = bl.parse_blast_xml(make_xml([TWO], qlen=100))["hits"][0]
    assert h["coverage_percent"] == pytest.approx(70.0)       # 1-40 ve 31-70 birleşimi
    assert h["bit_score"] == 60 and h["evalue"] == pytest.approx(1e-12)
    assert h["alignment"]["query_start"] == 31                 # en iyi HSP


def test_evalue_filtresi_ve_max_sonuc():
    xml = make_xml([PIG, HUMAN_NP, WEAK])
    assert [h["accession"] for h in bl.parse_blast_xml(xml, evalue=0.001)["hits"]] == ["NP_000198", "P01315"]
    r = bl.parse_blast_xml(xml, evalue=0.001, max_results=1)
    assert len(r["hits"]) == 1 and r["total_hits"] == 2


def test_hit_yok_bos_liste():
    r = bl.parse_blast_xml(make_xml([]))
    assert r["hits"] == [] and r["total_hits"] == 0


@pytest.mark.parametrize("bad", ["", "<html>Error</html>", "<BlastOutput><kirik>"])
def test_bozuk_xml_hata_verir(bad):
    with pytest.raises(bl.BlastError, match="beklenmeyen"):
        bl.parse_blast_xml(bad)


# ---------- hizalama biçimlendirme ----------
def test_format_alignment_bloklar_ve_numaralar():
    a = bl.parse_blast_xml(make_xml([PIG]))["hits"][0]["alignment"]
    text = bl.format_alignment(a, width=20)
    lines = text.splitlines()
    assert lines[0] == "Sorgu  25  FVNQHLCGSHLVEALYLVCG  44"
    assert lines[2] == "Hedef  25  FVNQHLCGSHLVEALYLVCG  44"
    assert "Sorgu  45  ERGFFYTPKT  54" in lines and "Hedef  45  ERGFFYTPKA  54" in lines


def test_format_alignment_bosluk_konumu_saymaz():
    a = dict(query_seq="AC-DE", midline="|| ||", subject_seq="ACGDE", query_start=1,
             query_end=4, subject_start=10, subject_end=14)
    lines = bl.format_alignment(a, width=60).splitlines()
    assert lines[0].endswith("AC-DE  4") and lines[2].endswith("ACGDE  14")


# ---------- run_blast ----------
def test_run_blast_basarili_ve_cagri_argumanlari():
    calls = []
    r = bl.run_blast(INS, database="swissprot", evalue=0.01, max_results=25,
                     qblast=fake_qblast(make_xml([PIG, HUMAN_NP]), calls=calls))
    c = calls[0]
    assert (c["program"], c["database"], c["sequence"]) == ("blastp", "swissprot", INS)
    assert c["expect"] == 0.01 and c["hitlist_size"] == 25 and c["format_type"] == "XML"
    assert len(r["hits"]) == 2 and r["evalue_threshold"] == 0.01 and r["elapsed_seconds"] >= 0


def test_run_blast_diziyi_temizler_belirsiz_harf_kalir():
    calls = []
    bl.run_blast(">x\n" + INS[:30].lower() + "X*", qblast=fake_qblast(make_xml([]), calls=calls))
    assert calls[0]["sequence"] == INS[:30] + "X"


def test_zaman_asimi():
    t0 = time.monotonic()
    with pytest.raises(bl.BlastTimeout, match="yanıt vermedi"):
        bl.run_blast(INS, timeout=0.3, qblast=fake_qblast(make_xml([]), delay=2))
    assert time.monotonic() - t0 < 1.5


def test_zaman_asimi_sunucu_hatasi_ayni_siniftan():
    assert issubclass(bl.BlastTimeout, bl.BlastServerError)


def test_iptal():
    ev = threading.Event()
    threading.Timer(0.2, ev.set).start()
    t0 = time.monotonic()
    with pytest.raises(bl.BlastCancelled, match="iptal"):
        bl.run_blast(INS, cancel_event=ev, qblast=fake_qblast(make_xml([]), delay=2))
    assert time.monotonic() - t0 < 1.5


def test_onceden_iptal_edilmis_hic_baslamaz():
    ev = threading.Event()
    ev.set()
    calls = []
    with pytest.raises(bl.BlastCancelled):
        bl.run_blast(INS, cancel_event=ev, qblast=fake_qblast(make_xml([]), calls=calls))
    assert calls == []


def test_on_tick_gecen_sureyi_bildirir():
    ticks = []
    bl.run_blast(INS, on_tick=ticks.append, tick_interval=0.1,
                 qblast=fake_qblast(make_xml([]), delay=0.6))
    assert len(ticks) >= 3 and ticks == sorted(ticks) and ticks[-1] <= 1.5


def test_internet_yok():
    with pytest.raises(bl.BlastConnectionError, match="İnternet"):
        bl.run_blast(INS, qblast=raising(urllib.error.URLError("dns")))
    with pytest.raises(bl.BlastConnectionError):
        bl.run_blast(INS, qblast=raising(OSError("x")))


def test_http_hatasi_sunucu_hatasi():
    err = urllib.error.HTTPError("u", 503, "x", {}, None)
    with pytest.raises(bl.BlastServerError, match="yanıt vermedi"):
        bl.run_blast(INS, qblast=raising(err))


def test_ncbi_hata_mesaji_genel_hata():
    with pytest.raises(bl.BlastError, match="başarısız") as e:
        bl.run_blast(INS, qblast=raising(ValueError("Error message from NCBI: bad")))
    assert type(e.value) is bl.BlastError


def test_html_cevabi_beklenmeyen_cevap():
    with pytest.raises(bl.BlastError, match="beklenmeyen"):
        bl.run_blast(INS, qblast=lambda *a, **k: io.StringIO("<html>503</html>"))


def test_bytes_cevap_kabul_edilir():
    r = bl.run_blast(INS, qblast=lambda *a, **k: io.BytesIO(make_xml([PIG]).encode()))
    assert r["hits"][0]["accession"] == "P01315"


# ---------- doğrulama ----------
@pytest.mark.parametrize("kw", [
    dict(database="yok"), dict(evalue=0), dict(evalue=-1), dict(evalue=True),
    dict(max_results=0), dict(max_results=2.5), dict(timeout=0), dict(email="mailyok"),
])
def test_gecersiz_parametre(kw):
    with pytest.raises(ValueError):
        bl.run_blast(INS, qblast=fake_qblast(make_xml([])), **kw)


def test_kisa_ve_gecersiz_dizi():
    with pytest.raises(pc.SequenceError, match="çok kısa"):
        bl.run_blast("ACDEF", qblast=fake_qblast(make_xml([])))
    with pytest.raises(pc.SequenceError):
        bl.run_blast("ACDEFGHIKL1", qblast=fake_qblast(make_xml([])))
    with pytest.raises(pc.SequenceError):
        bl.run_blast("", qblast=fake_qblast(make_xml([])))


def test_hata_siniflari_valueerror_alt_sinifi():
    for cls in (bl.BlastError, bl.BlastConnectionError, bl.BlastServerError,
                bl.BlastTimeout, bl.BlastCancelled):
        assert issubclass(cls, ValueError)
