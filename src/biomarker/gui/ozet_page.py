"""Özet ekranı: kartlar, hidropati profili, ikincil yapı bandı, motif listesi."""
import pyqtgraph as pg
from PyQt5.QtCore import Qt, QRectF, pyqtSignal
from PyQt5.QtGui import QColor, QPainter, QPainterPath
from PyQt5.QtWidgets import QFrame, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

import theme

SS_COLORS = {"h": theme.COMMON["helix"], "s": theme.COMMON["sheet"], "t": theme.COMMON["turn"]}

# ÖRNEK VERİ: motif taraması yazılınca gerçek sonuçla değişecek
SAMPLE_MOTIFS = [
    ("PS00262", "İnsülin ailesi imzası", "95–109"),
    ("PS00005", "PKC fosforilasyon bölgesi", "51–53, 54–56"),
    ("PS00006", "CK2 fosforilasyon bölgesi", "54–57"),
]


def tr_num(x: float, decimals: int, sign: bool = False) -> str:
    s = f"{x:+.{decimals}f}" if sign else f"{x:.{decimals}f}"
    return s.replace(".", ",")


class SSBand(QWidget):
    """Tek bantlı ikincil yapı haritası."""

    def __init__(self):
        super().__init__()
        self.setFixedHeight(22)
        self.segments = []

    def set_segments(self, segs):
        self.segments = segs
        self.update()

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        r = QRectF(self.rect())
        path = QPainterPath()
        path.addRoundedRect(r, 4, 4)
        p.setClipPath(path)
        p.fillRect(r, QColor(SS_COLORS["t"]))  # kıvrım: segment dışı bölgeler
        for a, b, k in self.segments:
            p.fillRect(
                QRectF(r.x() + a * r.width(), r.y(), (b - a) * r.width() + 0.5, r.height()),
                QColor(SS_COLORS[k]),
            )


def _legend_item(color, text):
    w = QWidget()
    lay = QHBoxLayout(w)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(5)
    box = QFrame()
    box.setFixedSize(10, 10)
    box.setStyleSheet(f"background:{color}; border:none; border-radius:2px;")
    lab = QLabel(text)
    lab.setObjectName("muted")
    lay.addWidget(box)
    lay.addWidget(lab)
    return w


class OzetPage(QWidget):
    navigate = pyqtSignal(str)  # sidebar anahtarı

    def __init__(self):
        super().__init__()
        self.mode = "dark"
        root = QVBoxLayout(self)
        root.setContentsMargins(20, 16, 20, 16)
        root.setSpacing(12)

        # başlık satırı
        head = QHBoxLayout()
        title = QLabel("Özet")
        title.setObjectName("pageTitle")
        self.badge = QLabel("")
        self.badge.setObjectName("badge")
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(self.badge)
        root.addLayout(head)

        # özet kartları
        row = QHBoxLayout()
        row.setSpacing(10)
        self.vals, self.notes = {}, {}
        for key, label in (
            ("protein", "Protein"), ("uzunluk", "Uzunluk"), ("mw", "Mol. ağırlık"),
            ("pi", "pI"), ("gravy", "GRAVY"), ("stab", "Stabilite indeksi"),
        ):
            row.addWidget(self._metric(key, label), 1)
        root.addLayout(row)

        # alt bölüm
        body = QHBoxLayout()
        body.setSpacing(10)

        left = QVBoxLayout()
        left.setSpacing(10)
        left.addWidget(self._build_hydro_card(), 1)
        left.addWidget(self._build_ss_card())
        body.addLayout(left, 3)
        body.addWidget(self._build_motif_card(), 2)
        root.addLayout(body, 1)

    # ---------- bileşenler ----------
    def _metric(self, key, label):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(14, 10, 14, 10)
        lay.setSpacing(2)
        lab = QLabel(label)
        lab.setObjectName("muted")
        val = QLabel("–")
        val.setObjectName("metric")
        note = QLabel("")
        note.setMinimumHeight(16)
        note.setObjectName("muted")
        lay.addWidget(lab)
        lay.addWidget(val)
        lay.addWidget(note)
        self.vals[key], self.notes[key] = val, note
        return card

    def _card(self, title, link_key=None, tag=None):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(8)
        head = QHBoxLayout()
        t = QLabel(title)
        t.setObjectName("cardTitle")
        head.addWidget(t)
        if tag:
            tl = QLabel(tag)
            tl.setObjectName("muted")
            head.addWidget(tl)
        head.addStretch(1)
        if link_key:
            b = QPushButton("Ayrıntı →")
            b.setObjectName("link")
            b.setCursor(Qt.PointingHandCursor)
            b.clicked.connect(lambda _=False, k=link_key: self.navigate.emit(k))
            head.addWidget(b)
        lay.addLayout(head)
        return card, lay

    def _build_hydro_card(self):
        card, lay = self._card("Hidropati profili", "fizikokimya")
        top = QLabel("hidrofobik")
        top.setObjectName("muted")
        bottom = QLabel("hidrofilik")
        bottom.setObjectName("muted")
        self.plot = pg.PlotWidget()
        self.plot.setMenuEnabled(False)
        self.plot.setMouseEnabled(False, False)
        self.plot.hideButtons()
        self.plot.hideAxis("left")
        self.plot.hideAxis("bottom")
        self.plot.setFrameShape(QFrame.NoFrame)
        self.plot.setMinimumHeight(120)
        self.plot.setYRange(-4.5, 4.5, padding=0)
        self.zero = pg.InfiniteLine(pos=0, angle=0)
        self.plot.addItem(self.zero)
        self.curve = self.plot.plot([], [])
        lay.addWidget(top)
        lay.addWidget(self.plot, 1)
        lay.addWidget(bottom)
        return card

    def _build_ss_card(self):
        card, lay = self._card("İkincil yapı", "ikincil")
        self.band = SSBand()
        lay.addWidget(self.band)
        leg = QHBoxLayout()
        leg.setSpacing(14)
        leg.addWidget(_legend_item(SS_COLORS["h"], "Heliks"))
        leg.addWidget(_legend_item(SS_COLORS["s"], "Yaprak"))
        leg.addWidget(_legend_item(SS_COLORS["t"], "Kıvrım"))
        leg.addStretch(1)
        lay.addLayout(leg)
        return card

    def _build_motif_card(self):
        card, lay = self._card("Motifler", "motif")
        for pid, desc, rng in SAMPLE_MOTIFS:
            sep = QFrame()
            sep.setObjectName("sep")
            sep.setFixedHeight(1)
            lay.addWidget(sep)
            row = QHBoxLayout()
            row.setContentsMargins(0, 4, 0, 4)
            col = QVBoxLayout()
            col.setSpacing(0)
            a = QLabel(pid)
            a.setObjectName("cardTitle")
            b = QLabel(desc)
            b.setObjectName("muted")
            col.addWidget(a)
            col.addWidget(b)
            r = QLabel(rng)
            row.addLayout(col)
            row.addStretch(1)
            row.addWidget(r)
            lay.addLayout(row)
        lay.addStretch(1)
        return card

    # ---------- veri ve tema ----------
    def set_result(self, info: dict):
        a = info["analysis"]
        self.badge.setText(f"{info['label']} · {a['length']} aa")
        self.vals["protein"].setText(info["accession"])
        self.vals["uzunluk"].setText(f"{a['length']} aa")
        self.vals["mw"].setText(f"{tr_num(a['mw_kda'], 2)} kDa")
        self.vals["pi"].setText(tr_num(a["pi"], 2))
        self.vals["gravy"].setText(tr_num(a["gravy"], 2, sign=True))
        self.vals["stab"].setText(tr_num(a["instability"], 1))
        stable = a["instability"] < 40
        t = theme.COMMON
        self.notes["stab"].setText("kararlı" if stable else "kararsız")
        self.notes["stab"].setStyleSheet(f"color:{t['success'] if stable else t['warning']};")

        n = max(a["length"], 1)
        self.band.set_segments(
            [((g["start"] - 1) / n, g["end"] / n, g["type"]) for g in info.get("ss", [])]
        )

        w = a["window"]
        xs = [i + (w + 1) / 2 for i in range(len(a["profile"]))]
        self.curve.setData(xs, a["profile"])
        self.plot.setXRange(1, max(a["length"], 2), padding=0)

    def apply_theme(self, mode: str):
        self.mode = mode
        t = theme.tokens(mode)
        self.plot.setBackground(t["card"])
        self.curve.setPen(pg.mkPen(t["accent"], width=2))
        self.zero.setPen(pg.mkPen(t["muted"], style=Qt.DashLine))