"""Özet ekranı: kartlar, hidropati profili, ikincil yapı bandı, motif listesi."""
import pyqtgraph as pg
from PyQt5.QtCore import Qt, QRectF, pyqtSignal
from PyQt5.QtGui import QColor, QPainter, QPainterPath
from PyQt5.QtWidgets import QFrame, QHBoxLayout, QLabel, QPushButton, QVBoxLayout, QWidget

import theme

SS_COLORS = {"h": theme.COMMON["helix"], "s": theme.COMMON["sheet"], "t": theme.COMMON["turn"]}



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
        self.motif_rows = QVBoxLayout()
        self.motif_rows.setSpacing(0)
        lay.addLayout(self.motif_rows)
        lay.addStretch(1)
        return card

    def _fill_motifs(self, motifs):
        while self.motif_rows.count():
            w = self.motif_rows.takeAt(0).widget()
            if w:
                w.setParent(None)
                w.deleteLater()
        groups = {}
        for m in motifs:
            g = groups.setdefault(m["id"], {"name": m["name"], "ranges": []})
            g["ranges"].append(f"{m['start']}–{m['end']}")
        if not groups:
            none = QLabel("Motif bulunamadı.")
            none.setObjectName("muted")
            self.motif_rows.addWidget(none)
            return
        for mid, g in groups.items():
            sep = QFrame()
            sep.setObjectName("sep")
            sep.setFixedHeight(1)
            self.motif_rows.addWidget(sep)
            row = QWidget()
            rl = QHBoxLayout(row)
            rl.setContentsMargins(0, 6, 0, 6)
            col = QVBoxLayout()
            col.setSpacing(0)
            a = QLabel(mid)
            a.setObjectName("cardTitle")
            b = QLabel(g["name"])
            b.setObjectName("muted")
            col.addWidget(a)
            col.addWidget(b)
            rl.addLayout(col)
            rl.addStretch(1)
            rl.addWidget(QLabel(", ".join(g["ranges"])))
            self.motif_rows.addWidget(row)

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

        self._fill_motifs(info.get("motifs", []))
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