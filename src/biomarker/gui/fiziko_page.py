"""Fizikokimyasal analiz sayfası: kartlar, hidropati profili, titrasyon eğrisi, aminoasit bileşimi."""
import pyqtgraph as pg
from pyqtgraph.exporters import ImageExporter
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QFont
from PyQt5.QtWidgets import (
    QComboBox, QFileDialog, QFrame, QHBoxLayout, QLabel, QPushButton,
    QScrollArea, QSpinBox, QVBoxLayout, QWidget,
)

import analysis
import theme
from ozet_page import tr_num

SCALES = ["Kyte-Doolittle", "Hopp-Woods"]
YRANGE = {"Kyte-Doolittle": (-4.5, 4.5), "Hopp-Woods": (-3.5, 3.5)}
AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"


class FizikoPage(QWidget):
    def __init__(self):
        super().__init__()
        self.info = None
        self.seq = ""
        self.mode = "dark"
        self._profile = []
        self._window = 9
        self.bars = None

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.NoFrame)
        outer.addWidget(scroll)
        inner = QWidget()
        inner.setObjectName("workspace")
        scroll.setWidget(inner)

        root = QVBoxLayout(inner)
        root.setContentsMargins(20, 16, 20, 16)
        root.setSpacing(12)

        head = QHBoxLayout()
        title = QLabel("Fizikokimyasal analiz")
        title.setObjectName("pageTitle")
        self.badge = QLabel("")
        self.badge.setObjectName("badge")
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(self.badge)
        root.addLayout(head)

        # kartlar
        row = QHBoxLayout()
        row.setSpacing(10)
        self.vals, self.notes = {}, {}
        for key, label in (
            ("uzunluk", "Uzunluk"), ("mw", "Mol. ağırlık"), ("pi", "pI"),
            ("gravy", "GRAVY"), ("stab", "Stabilite indeksi"),
        ):
            row.addWidget(self._metric(key, label), 1)
        root.addLayout(row)

        root.addWidget(self._build_hydro_card())

        bottom = QHBoxLayout()
        bottom.setSpacing(10)
        bottom.addWidget(self._build_titration_card(), 1)
        bottom.addWidget(self._build_comp_card(), 1)
        root.addLayout(bottom)
        root.addStretch(1)

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
        note.setObjectName("muted")
        note.setMinimumHeight(16)
        lay.addWidget(lab)
        lay.addWidget(val)
        lay.addWidget(note)
        self.vals[key], self.notes[key] = val, note
        return card

    def _plain_plot(self, height):
        p = pg.PlotWidget()
        p.setMenuEnabled(False)
        p.hideButtons()
        p.setFrameShape(QFrame.NoFrame)
        p.setMinimumHeight(height)
        f = QFont()
        f.setPointSize(9)
        for ax in ("left", "bottom"):
            p.getAxis(ax).setTickFont(f)
        return p

    def _build_hydro_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(6)

        head = QHBoxLayout()
        head.setSpacing(8)
        t = QLabel("Hidropati profili")
        t.setObjectName("cardTitle")
        self.scale_box = QComboBox()
        self.scale_box.addItems(SCALES)
        pl = QLabel("Pencere")
        pl.setObjectName("muted")
        self.win_spin = QSpinBox()
        self.win_spin.setRange(3, 31)
        self.win_spin.setSingleStep(2)
        self.win_spin.setValue(9)
        png = QPushButton("PNG kaydet")
        png.setObjectName("secondary")
        png.setCursor(Qt.PointingHandCursor)
        png.clicked.connect(self._save_png)
        head.addWidget(t)
        head.addStretch(1)
        head.addWidget(self.scale_box)
        head.addWidget(pl)
        head.addWidget(self.win_spin)
        head.addWidget(png)
        lay.addLayout(head)

        self.hover = QLabel("İmleci grafiğin üzerine getirin")
        self.hover.setObjectName("muted")
        top = QLabel("hidrofobik")
        top.setObjectName("muted")
        bottom = QLabel("hidrofilik")
        bottom.setObjectName("muted")

        self.h_plot = self._plain_plot(170)
        self.h_plot.setMouseEnabled(x=True, y=False)
        self.h_plot.setLabel("bottom", "Amino asit konumu", **{"font-size": "9pt"})
        self.h_plot.setLabel("left", "Skor", **{"font-size": "9pt"})
        self.h_zero = pg.InfiniteLine(pos=0, angle=0)
        self.h_plot.addItem(self.h_zero)
        self.h_curve = self.h_plot.plot([], [])
        self.h_plot.scene().sigMouseMoved.connect(self._on_hover)
        self.h_plot.scene().sigMouseClicked.connect(self._on_click)

        lay.addWidget(self.hover)
        lay.addWidget(top)
        lay.addWidget(self.h_plot, 1)
        lay.addWidget(bottom)

        self.scale_box.currentIndexChanged.connect(self._reload_profile)
        self.win_spin.valueChanged.connect(self._on_window)
        return card

    def _build_titration_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(6)
        t = QLabel("Titrasyon eğrisi")
        t.setObjectName("cardTitle")
        lay.addWidget(t)

        self.t_plot = self._plain_plot(150)
        self.t_plot.setMouseEnabled(x=False, y=False)
        self.t_plot.hideAxis("left")
        self.t_plot.setLabel("bottom", "pH", **{"font-size": "9pt"})
        self.t_plot.getAxis("bottom").setTicks([[(0, "0"), (7, "7"), (14, "14")]])
        self.t_plot.setXRange(0, 14, padding=0.06)
        self.t_zero = pg.InfiniteLine(pos=0, angle=0)
        self.t_plot.addItem(self.t_zero)
        self.t_up = self.t_plot.plot([], [])
        self.t_down = self.t_plot.plot([], [])
        self.t_dot = pg.ScatterPlotItem(size=9)
        self.t_plot.addItem(self.t_dot)
        self.t_pi = pg.TextItem("", anchor=(0, 1))
        self.t_plus = pg.TextItem("+", anchor=(0, 0.5))
        self.t_minus = pg.TextItem("−", anchor=(0, 0.5))
        for it in (self.t_pi, self.t_plus, self.t_minus):
            self.t_plot.addItem(it)
        lay.addWidget(self.t_plot, 1)
        return card

    def _build_comp_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(6)
        head = QHBoxLayout()
        t = QLabel("Aminoasit bileşimi")
        t.setObjectName("cardTitle")
        sub = QLabel("% · alfabetik")
        sub.setObjectName("muted")
        head.addWidget(t)
        head.addStretch(1)
        head.addWidget(sub)
        lay.addLayout(head)

        self.c_plot = self._plain_plot(150)
        self.c_plot.setMouseEnabled(x=False, y=False)
        self.c_plot.hideAxis("left")
        self.c_plot.getAxis("bottom").setTicks([list(enumerate(AA_ORDER))])
        self.c_plot.setXRange(-0.7, 19.7, padding=0)
        lay.addWidget(self.c_plot, 1)
        return card

    # ---------- veri ----------
    def set_result(self, info: dict):
        self.info = info
        self.seq = info.get("seq", "")
        a = info["analysis"]
        self.badge.setText(f"{info['label']} · {a['length']} aa")
        self.vals["uzunluk"].setText(f"{a['length']} aa")
        self.vals["mw"].setText(f"{tr_num(a['mw_kda'], 2)} kDa")
        self.vals["pi"].setText(tr_num(a["pi"], 2))
        self.vals["gravy"].setText(tr_num(a["gravy"], 2, sign=True))
        self.vals["stab"].setText(tr_num(a["instability"], 1))
        stable = a["instability"] < 40
        c = theme.COMMON
        self.notes["stab"].setText("kararlı" if stable else "kararsız")
        self.notes["stab"].setStyleSheet(f"color:{c['success'] if stable else c['warning']};")

        self._reload_profile()
        self._load_titration()
        self._load_composition()
        self._restyle()

    def _on_window(self, v):
        if v % 2 == 0:  # yalnızca tek sayı
            self.win_spin.blockSignals(True)
            self.win_spin.setValue(v + 1)
            self.win_spin.blockSignals(False)
        self._reload_profile()

    def _reload_profile(self, *_):
        if not self.info:
            return
        scale = self.scale_box.currentText()
        w = max(1, min(self.win_spin.value(), len(self.seq) or 1))
        self._window = w
        self._profile = analysis.hydropathy_profile(self.seq, scale, w)
        xs = [i + (w + 1) // 2 for i in range(len(self._profile))]
        self.h_curve.setData(xs, self._profile)
        n = max(len(self.seq), 2)
        lo, hi = YRANGE[scale]
        self.h_plot.setLimits(xMin=1, xMax=n)
        self.h_plot.setYRange(lo, hi, padding=0)
        self.h_plot.setXRange(1, n, padding=0)

    def _load_titration(self):
        d = analysis.titration_curve(self.seq)
        pi, ph, ch = d["pi"], d["ph"], d["charge"]
        up = [(p, c) for p, c in zip(ph, ch) if p < pi] + [(pi, 0)]
        down = [(pi, 0)] + [(p, c) for p, c in zip(ph, ch) if p > pi]
        self.t_up.setData([p for p, _ in up], [c for _, c in up])
        self.t_down.setData([p for p, _ in down], [c for _, c in down])
        self.t_dot.setData([pi], [0])
        m = max(abs(c) for c in ch) or 1
        self.t_plot.setYRange(-m * 1.15, m * 1.15, padding=0)
        self.t_pi.setText(f"pI = {tr_num(pi, 2)}")
        self.t_pi.setPos(pi + 0.3, 0)
        self.t_plus.setPos(0.1, m)
        self.t_minus.setPos(0.1, -m)

    def _load_composition(self):
        d = analysis.aa_composition(self.seq)
        vals = [d[a] for a in AA_ORDER]
        if self.bars is not None:
            self.c_plot.removeItem(self.bars)
        self.bars = pg.BarGraphItem(x=list(range(20)), height=vals, width=0.7,
                                    brush=theme.tokens(self.mode)["accent"], pen=None)
        self.c_plot.addItem(self.bars)
        self.c_plot.setYRange(0, max(vals) * 1.15, padding=0)

    # ---------- etkileşim ----------
    def _on_hover(self, pos):
        if not self._profile or not self.h_plot.sceneBoundingRect().contains(pos):
            return
        vb = self.h_plot.getPlotItem().vb
        res = int(round(vb.mapSceneToView(pos).x()))
        idx = res - (self._window + 1) // 2
        if 1 <= res <= len(self.seq) and 0 <= idx < len(self._profile):
            self.hover.setText(
                f"Konum {res} · {self.seq[res - 1]} · skor {tr_num(self._profile[idx], 2, sign=True)}"
            )

    def _on_click(self, ev):
        if ev.double() and self.seq:  # çift tık: yakınlaştırmayı sıfırla
            self.h_plot.setXRange(1, max(len(self.seq), 2), padding=0)

    def _save_png(self):
        path, _ = QFileDialog.getSaveFileName(self, "Grafiği kaydet", "hidropati.png", "PNG (*.png)")
        if not path:
            return
        try:
            ImageExporter(self.h_plot.plotItem).export(path)
        except Exception:
            self.hover.setText("Grafik kaydedilemedi. Dosya yolunu kontrol edin.")

    # ---------- tema ----------
    def apply_theme(self, mode: str):
        self.mode = mode
        self._restyle()

    def _restyle(self):
        t = theme.tokens(mode := self.mode)
        muted = pg.mkColor(t["muted"])
        for plot in (self.h_plot, self.t_plot, self.c_plot):
            plot.setBackground(t["card"])
            for ax in ("left", "bottom"):
                a = plot.getAxis(ax)
                a.setPen(muted)
                a.setTextPen(muted)
        dash = pg.mkPen(t["muted"], style=Qt.DashLine)
        self.h_zero.setPen(dash)
        self.t_zero.setPen(dash)
        self.h_curve.setPen(pg.mkPen(t["accent"], width=2))
        self.t_up.setPen(pg.mkPen(t["charge_pos"], width=2))
        self.t_down.setPen(pg.mkPen(t["charge_neg"], width=2))
        self.t_dot.setBrush(pg.mkBrush(t["charge_pos"]))
        self.t_dot.setPen(pg.mkPen(None))
        for it in (self.t_pi, self.t_plus, self.t_minus):
            it.setColor(t["text"])
        if self.bars is not None:
            self.bars.setOpts(brush=pg.mkBrush(t["accent"]))
