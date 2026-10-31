"""homology (EMBL-EBI BLAST) testleri. Ağa çıkmaz: sahte HTTP açıcı ve elle kurulmuş BLAST XML'i kullanır.
Çalıştırma: python -m pytest tests -q"""
import io
import pathlib
import sys
import threading
import time
import socket
import urllib.error
import urllib.parse

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from biomarker.analysis import homology as bl  # noqa: E402
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


EBI_HUMAN = dict(id="sp|P01315|INS_PIG", acc="P01315", len=108, hsps=[hsp(61.2, "2e-75")],
                 **{"def": "Insulin OS=Sus scrofa OX=9823 GN=INS PE=1 SV=2"})
EBI_NOACC = dict(id="sp|P01317|INS_BOVIN", acc="", len=105, hsps=[hsp(55.0, "1e-60")],
                 **{"def": "Insulin OS=Bos taurus OX=9913 GN=INS PE=1 SV=1"})


# ---------- başlık ve erişim no biçimleri ----------
@pytest.mark.parametrize("title,expected", [
    ("Insulin OS=Sus scrofa OX=9823 GN=INS PE=1 SV=2", ("Insulin", "Sus scrofa")),
    ("Insulin-like growth factor I OS=Homo sapiens OX=9606 PE=1 SV=1",
     ("Insulin-like growth factor I", "Homo sapiens")),
    ("insulin preproprotein [Homo sapiens] >gb|AAA59172.1| insulin [Homo sapiens]",
     ("insulin preproprotein", "Homo sapiens")),
    ("RecName: Full=Insulin; Flags: Precursor [Sus scrofa]", ("Insulin", "Sus scrofa")),
    ("sadece ad", ("sadece ad", "")),
    ("", ("", "")),
])
def test_baslik_cozumleme(title, expected):
    assert bl._clean_title(title) == expected


def test_ebi_bicimli_hit_ve_erisim_no_kimlikten():
    r = bl.parse_blast_xml(make_xml([EBI_HUMAN, EBI_NOACC]))
    by = {h["accession"]: h for h in r["hits"]}
    assert set(by) == {"P01315", "P01317"}                 # ikincisi sp|P01317|INS_BOVIN'den
    assert by["P01315"]["organism"] == "Sus scrofa" and by["P01317"]["name"] == "Insulin"
    assert by["P01317"]["url"] == "https://www.uniprot.org/uniprotkb/P01317"


# ---------- sahte EBI sunucusu ----------
JOB = "ncbiblast-R20261031-120000-0001-12345678-p1m"


class FakeResp:
    def __init__(self, text):
        self._data = text.encode("utf-8")

    def read(self):
        return self._data

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


class FakeEBI:
    """run / status / result uç noktalarını taklit eder. Çağrıları kaydeder."""

    def __init__(self, xml="", statuses=("FINISHED",), run_reply=JOB, errors=None):
        self.xml, self.statuses, self.run_reply = xml, list(statuses), run_reply
        self.errors = errors or {}          # {"run"|"status"|"xml": Exception}
        self.calls = []

    def __call__(self, req, timeout=None):
        url = req.full_url
        data = urllib.parse.parse_qs(req.data.decode()) if req.data else None
        self.calls.append((req.get_method(), url, data))
        kind = "run" if url.endswith("/run") else "status" if "/status/" in url else "xml"
        if kind in self.errors:
            raise self.errors[kind]
        if kind == "run":
            return FakeResp(self.run_reply)
        if kind == "status":
            s = self.statuses.pop(0) if len(self.statuses) > 1 else self.statuses[0]
            return FakeResp(s)
        return FakeResp(self.xml)

    def count(self, kind):
        key = {"run": "/run", "status": "/status/", "xml": "/xml"}[kind]
        return sum(1 for _, u, _ in self.calls if key in u)


def run(fake, **kw):
    kw.setdefault("email", "ben@ornek.com")
    kw.setdefault("poll_interval", 0.02)
    return bl.run_blast(INS, opener=fake, **kw)


# ---------- run_blast ----------
def test_run_blast_basarili_ve_istek_icerigi():
    fake = FakeEBI(make_xml([EBI_HUMAN, HUMAN_NP]), statuses=["QUEUED", "RUNNING", "FINISHED"])
    r = run(fake)
    method, url, data = fake.calls[0]
    assert (method, url) == ("POST", "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast/run")
    assert set(data) == {"email", "program", "stype", "database", "sequence"}   # yalnız zorunlular
    assert data["program"] == ["blastp"] and data["stype"] == ["protein"]
    assert data["database"] == ["uniprotkb_swissprot"] and data["email"] == ["ben@ornek.com"]
    assert data["sequence"] == [">BioMarker\n" + INS]
    assert fake.count("status") == 3 and fake.count("xml") == 1
    assert fake.calls[-1][1].endswith(f"/result/{JOB}/xml")
    assert r["job_id"] == JOB and r["provider"] == "EMBL-EBI" and r["database_label"] == "Swiss-Prot"
    assert r["evalue_threshold"] == 0.001 and r["elapsed_seconds"] >= 0
    assert {h["accession"] for h in r["hits"]} == {"P01315", "NP_000198"}


def test_veritabani_eslemesi():
    fake = FakeEBI(make_xml([]))
    run(fake, database="uniprotkb")
    assert fake.calls[0][2]["database"] == ["uniprotkb"]


def test_evalue_ve_max_sonuc_istemci_tarafinda_uygulanir():
    fake = FakeEBI(make_xml([EBI_HUMAN, HUMAN_NP, WEAK]))
    r = run(fake, evalue=0.001, max_results=1)
    assert len(r["hits"]) == 1 and r["total_hits"] == 2          # WEAK (0.5) elendi
    assert r["hits"][0]["accession"] == "NP_000198"              # bit skoru en yüksek


def test_diziyi_temizler_belirsiz_harf_kalir():
    fake = FakeEBI(make_xml([]))
    bl.run_blast(">x\n" + INS[:30].lower() + "X*", email="a@b.com", opener=fake, poll_interval=0.02)
    assert fake.calls[0][2]["sequence"] == [">BioMarker\n" + INS[:30] + "X"]


# ---------- durumlar ve hatalar ----------
def test_is_basarisiz():
    with pytest.raises(bl.BlastServerError, match="işleyemedi") as e:
        run(FakeEBI(statuses=["RUNNING", "FAILURE"]))
    assert "FAILURE" in e.value.detail


def test_is_bulunamadi():
    with pytest.raises(bl.BlastServerError, match="bulunamadı"):
        run(FakeEBI(statuses=["NOT_FOUND"]))


def test_beklenmeyen_durum():
    with pytest.raises(bl.BlastError, match="beklenmeyen") as e:
        run(FakeEBI(statuses=["<html>bakim</html>"]))
    assert "BAKIM" in e.value.detail.upper()


def test_gecersiz_is_kimligi():
    with pytest.raises(bl.BlastError, match="beklenmeyen") as e:
        run(FakeEBI(run_reply="Invalid email address"))
    assert "Invalid email" in e.value.detail


def test_xml_degil_beklenmeyen_cevap():
    with pytest.raises(bl.BlastError, match="beklenmeyen") as e:
        run(FakeEBI("<EBIApplicationResult>baska bicim</EBIApplicationResult>"))
    assert "EBIApplicationResult" in e.value.detail


def http_error(code):
    return urllib.error.HTTPError("u", code, "x", {}, io.BytesIO(b"sunucu mesaji"))


def test_http_400_istek_reddedildi_ayrinti_korunur():
    with pytest.raises(bl.BlastError, match="reddetti") as e:
        run(FakeEBI(errors={"run": http_error(400)}))
    assert "HTTP 400" in e.value.detail and "sunucu mesaji" in e.value.detail
    assert type(e.value) is bl.BlastError


def test_http_429_cok_istek():
    with pytest.raises(bl.BlastError, match="Çok fazla"):
        run(FakeEBI(errors={"run": http_error(429)}))


def test_http_503_sunucu_hatasi():
    with pytest.raises(bl.BlastServerError, match="yanıt vermedi"):
        run(FakeEBI(statuses=["RUNNING"], errors={"status": http_error(503)}))


def test_internet_yok_ve_soket_zaman_asimi():
    with pytest.raises(bl.BlastConnectionError, match="İnternet"):
        run(FakeEBI(errors={"run": urllib.error.URLError("dns")}))
    with pytest.raises(bl.BlastConnectionError) as e:
        run(FakeEBI(errors={"run": socket.timeout("okuma")}))
    assert "timeout" in e.value.detail.lower()


# ---------- zaman aşımı, iptal, tick ----------
def test_zaman_asimi_gercekten_durur():
    fake = FakeEBI(statuses=["RUNNING"])
    t0 = time.monotonic()
    with pytest.raises(bl.BlastTimeout, match="yanıt vermedi") as e:
        run(fake, timeout=0.3)
    assert time.monotonic() - t0 < 1.0 and "0.3" in e.value.detail
    n = len(fake.calls)
    time.sleep(0.2)
    assert len(fake.calls) == n                # zaman aşımından sonra hiç istek gitmiyor


def test_zaman_asimi_sunucu_hatasi_ayni_siniftan():
    assert issubclass(bl.BlastTimeout, bl.BlastServerError)


def test_iptal_gercekten_durur():
    fake = FakeEBI(statuses=["RUNNING"])
    ev = threading.Event()
    threading.Timer(0.2, ev.set).start()
    t0 = time.monotonic()
    with pytest.raises(bl.BlastCancelled, match="iptal"):
        run(fake, cancel_event=ev)
    assert time.monotonic() - t0 < 1.0
    n = len(fake.calls)
    time.sleep(0.2)
    assert len(fake.calls) == n and fake.count("xml") == 0


def test_onceden_iptal_edilmis_hic_istek_gonderilmez():
    ev = threading.Event()
    ev.set()
    fake = FakeEBI(make_xml([]))
    with pytest.raises(bl.BlastCancelled):
        run(fake, cancel_event=ev)
    assert fake.calls == []


def test_on_tick_gecen_sureyi_bildirir():
    ticks = []
    fake = FakeEBI(make_xml([]), statuses=["RUNNING"] * 12 + ["FINISHED"])
    run(fake, on_tick=ticks.append, tick_interval=0.1, poll_interval=0.05)
    assert len(ticks) >= 3 and ticks == sorted(ticks)


def test_on_tick_hatasi_sorguyu_bozmaz():
    def bad(_):
        raise RuntimeError("arayuz hatasi")
    fake = FakeEBI(make_xml([EBI_HUMAN]), statuses=["RUNNING", "FINISHED"])
    r = run(fake, on_tick=bad, tick_interval=0.0)
    assert r["hits"][0]["accession"] == "P01315"


# ---------- doğrulama ----------
@pytest.mark.parametrize("kw", [
    dict(database="nr"), dict(database="yok"), dict(evalue=0), dict(evalue=-1), dict(evalue=True),
    dict(max_results=0), dict(max_results=2.5), dict(timeout=0), dict(poll_interval=0),
])
def test_gecersiz_parametre(kw):
    with pytest.raises(ValueError):
        run(FakeEBI(make_xml([])), **kw)


@pytest.mark.parametrize("email", [None, "", "mailyok"])
def test_eposta_zorunlu(email):
    fake = FakeEBI(make_xml([]))
    with pytest.raises(bl.BlastError, match="e-posta"):
        bl.run_blast(INS, email=email, opener=fake)
    assert fake.calls == []


def test_kisa_ve_gecersiz_dizi_istek_gonderilmez():
    fake = FakeEBI(make_xml([]))
    for bad in ("ACDEF", "ACDEFGHIKL1", ""):
        with pytest.raises(pc.SequenceError):
            bl.run_blast(bad, email="a@b.com", opener=fake)
    assert fake.calls == []


def test_hata_siniflari_valueerror_alt_sinifi_ve_detay_alani():
    for cls in (bl.BlastError, bl.BlastConnectionError, bl.BlastServerError,
                bl.BlastTimeout, bl.BlastCancelled):
        assert issubclass(cls, ValueError)
    assert bl.BlastError("x").detail is None and bl.BlastError("x", "d").detail == "d"
