"""İkincil yapı sayfası: Chou-Fasman bandı, işaretli dizi, segment tablosu."""
import math

from PyQt5.QtCore import QRectF, Qt, pyqtSignal
from PyQt5.QtGui import QColor, QFont, QPainter, QPainterPath, QPen
from PyQt5.QtWidgets import (
    QAbstractItemView, QCheckBox, QFrame, QHBoxLayout, QHeaderView, QLabel,
    QScrollArea, QTableWidget, QTableWidgetItem, QTextEdit, QVBoxLayout, QWidget,
)

import theme
from ozet_page import SS_COLORS, _legend_item

TYPE_NAMES = {"h": "Heliks", "s": "Yaprak"}
LINE = 40  # dizi satırında karakter sayısı


def blend(fg: str, bg: str, a: float) -> str:
    c1, c2 = QColor(fg), QColor(bg)
    return QColor(
        int(c1.red() * a + c2.red() * (1 - a)),
        int(c1.green() * a + c2.green() * (1 - a)),
        int(c1.blue() * a + c2.blue() * (1 - a)),
    ).name()


class SSBar(QWidget):
    """Tam genişlikte bant. Tıklanan segment seçilir."""
    clicked = pyqtSignal(int)

    def __init__(self):
        super().__init__()
        self.setFixedHeight(38)
        self.setCursor(Qt.PointingHandCursor)
        self.n, self.segs, self.sel = 1, [], -1
        self.outline = "#FFFFFF"

    def set_data(self, n, segs, sel, outline):
        self.n, self.segs, self.sel, self.outline = max(n, 1), segs, sel, outline
        self.update()

    def _rect_of(self, seg):
        w = self.width()
        x = (seg["start"] - 1) / self.n * w
        return QRectF(x, 0, (seg["end"] - seg["start"] + 1) / self.n * w, self.height())

    def paintEvent(self, _):
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        r = QRectF(self.rect())
        path = QPainterPath()
        path.addRoundedRect(r, 6, 6)
        p.setClipPath(path)
        p.fillRect(r, QColor(SS_COLORS["t"]))
        for seg in self.segs:
            p.fillRect(self._rect_of(seg), QColor(SS_COLORS[seg["type"]]))
        p.setClipping(False)
        if 0 <= self.sel < len(self.segs):
            p.setPen(QPen(QColor(self.outline), 2))
            p.setBrush(Qt.NoBrush)
            p.drawRect(self._rect_of(self.segs[self.sel]).adjusted(1, 1, -1, -1))

    def mousePressEvent(self, e):
        for i, seg in enumerate(self.segs):
            if self._rect_of(seg).contains(e.pos().x(), 1):
                self.clicked.emit(i)
                return


class IkincilPage(QWidget):
    def __init__(self):
        super().__init__()
        self.mode = "dark"
        self.seq = ""
        self.all_segs = []
        self.segs = []
        self.sel = -1

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
        title = QLabel("İkincil yapı tahmini")
        title.setObjectName("pageTitle")
        self.badge = QLabel("")
        self.badge.setObjectName("badge")
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(self.badge)
        root.addLayout(head)

        root.addWidget(self._build_band_card())
        root.addWidget(self._build_seq_card())
        root.addWidget(self._build_table_card())
        root.addStretch(1)

    # ---------- bileşenler ----------
    def _build_band_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(8)

        head = QHBoxLayout()
        head.setSpacing(14)
        t = QLabel("Chou-Fasman bandı")
        t.setObjectName("cardTitle")
        self.cb_core = QCheckBox("Çekirdek")
        self.cb_ext = QCheckBox("Genişletilmiş")
        self.cb_core.setChecked(True)
        self.cb_ext.setChecked(True)
        self.cb_core.toggled.connect(self._rebuild)
        self.cb_ext.toggled.connect(self._rebuild)
        head.addWidget(t)
        head.addStretch(1)
        head.addWidget(self.cb_core)
        head.addWidget(self.cb_ext)
        lay.addLayout(head)

        self.bar = SSBar()
        self.bar.clicked.connect(self._select)
        lay.addWidget(self.bar)

        ticks = QHBoxLayout()
        self.tick_l, self.tick_m, self.tick_r = QLabel("1"), QLabel(""), QLabel("")
        for lab, al in ((self.tick_l, Qt.AlignLeft), (self.tick_m, Qt.AlignCenter), (self.tick_r, Qt.AlignRight)):
            lab.setObjectName("muted")
            lab.setAlignment(al)
            ticks.addWidget(lab, 1)
        lay.addLayout(ticks)

        leg = QHBoxLayout()
        leg.setSpacing(14)
        leg.addWidget(_legend_item(SS_COLORS["h"], "Heliks"))
        leg.addWidget(_legend_item(SS_COLORS["s"], "Yaprak"))
        leg.addWidget(_legend_item(SS_COLORS["t"], "Kıvrım"))
        leg.addStretch(1)
        lay.addLayout(leg)

        sep = QFrame()
        sep.setObjectName("sep")
        sep.setFixedHeight(1)
        lay.addWidget(sep)
        self.summary = QLabel("")
        lay.addWidget(self.summary)
        return card

    def _build_seq_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(8)
        t = QLabel("Dizi")
        t.setObjectName("cardTitle")
        lay.addWidget(t)
        self.seq_view = QTextEdit()
        self.seq_view.setReadOnly(True)
        self.seq_view.setFrameShape(QFrame.NoFrame)
        self.seq_view.setStyleSheet("background: transparent;")
        lay.addWidget(self.seq_view)
        return card

    def _build_table_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(12, 8, 12, 8)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Tür", "Başlangıç", "Bitiş", "Uzunluk", "Evre"])
        self.table.verticalHeader().hide()
        self.table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.table.horizontalHeader().setHighlightSections(False)
        self.table.horizontalHeader().setDefaultAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.table.horizontalHeader().setFixedHeight(34)
        self.table.verticalHeader().setDefaultSectionSize(34)
        self.table.setSelectionBehavior(QAbstractItemView.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table.setShowGrid(False)
        self.table.setFocusPolicy(Qt.NoFocus)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.table.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.table.itemSelectionChanged.connect(self._on_table_select)
        lay.addWidget(self.table)
        return card

    # ---------- veri ----------
    def set_result(self, info: dict):
        self.seq = info.get("seq", "")
        a = info["analysis"]
        self.badge.setText(f"{info['label']} · {a['length']} aa")
        self.all_segs = list(info.get("ss", []))
        n = len(self.seq)
        self.tick_l.setText("1")
        self.tick_m.setText(str(n // 2))
        self.tick_r.setText(str(n))
        self.sel = 0
        self._rebuild()

    def _visible(self):
        ok = {"Çakışma çözüldü"}
        if self.cb_core.isChecked():
            ok.add("Çekirdek")
        if self.cb_ext.isChecked():
            ok.add("Genişletilmiş")
        return [s for s in self.all_segs if s["stage"] in ok]

    def _rebuild(self, *_):
        if not self.seq:
            return
        self.segs = self._visible()
        self.sel = 0 if self.segs else -1
        self._fill_table()
        self._update_summary()
        self._refresh_selection()

    def _fill_table(self):
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.segs))
        for r, s in enumerate(self.segs):
            vals = [TYPE_NAMES[s["type"]], s["start"], s["end"], s["end"] - s["start"] + 1, s["stage"]]
            for c, v in enumerate(vals):
                self.table.setItem(r, c, QTableWidgetItem(str(v)))
        self.table.setFixedHeight(34 + 34 * len(self.segs) + 4)
        self.table.blockSignals(False)

    def _update_summary(self):
        n = max(len(self.seq), 1)
        h = sum(s["end"] - s["start"] + 1 for s in self.segs if s["type"] == "h")
        sh = sum(s["end"] - s["start"] + 1 for s in self.segs if s["type"] == "s")
        hp, sp = round(h / n * 100), round(sh / n * 100)
        self.summary.setText(f"Heliks %{hp} · Yaprak %{sp} · Kıvrım %{100 - hp - sp}")

    # ---------- seçim ----------
    def _select(self, idx):
        self.sel = idx
        self._refresh_selection()

    def _on_table_select(self):
        rows = self.table.selectionModel().selectedRows()
        if rows:
            self.sel = rows[0].row()
            self._refresh_selection(from_table=True)

    def _refresh_selection(self, from_table=False):
        t = theme.tokens(self.mode)
        self.bar.set_data(len(self.seq), self.segs, self.sel, t["text"])
        if not from_table and 0 <= self.sel < self.table.rowCount():
            self.table.blockSignals(True)
            self.table.selectRow(self.sel)
            self.table.blockSignals(False)
        self._render_sequence()

    def _render_sequence(self):
        t = theme.tokens(self.mode)
        n = len(self.seq)
        owner = [-1] * n  # her kalıntının hangi segmente ait olduğu
        for i, s in enumerate(self.segs):
            for p in range(s["start"] - 1, s["end"]):
                owner[p] = i
        parts = []
        for start in range(0, n, LINE):
            end = min(start + LINE, n)
            parts.append(f'<span style="color:{t["muted"]}">{start + 1:>4}</span>  ')
            p = start
            while p < end:
                q = p
                while q < end and owner[q] == owner[p]:
                    q += 1
                o = owner[p]
                if o == -1:
                    bg = blend(SS_COLORS["t"], t["card"], 0.22)
                    style = f"background-color:{bg};"
                else:
                    kind = self.segs[o]["type"]
                    alpha = 0.6 if o == self.sel else 0.34
                    style = f"background-color:{blend(SS_COLORS[kind], t['card'], alpha)};"
                    if o == self.sel:
                        style += "text-decoration:underline;"
                parts.append(f'<span style="{style}">{self.seq[p:q]}</span>')
                p = q
            parts.append("\n")
        html = (
            f'<pre style="margin:0; color:{t["text"]}; '
            f"font-family:'Consolas','DejaVu Sans Mono',monospace; font-size:13px;\">"
            + "".join(parts).rstrip("\n") + "</pre>"
        )
        self.seq_view.setHtml(html)
        lines = math.ceil(n / LINE)
        self.seq_view.setFixedHeight(min(max(lines * 22 + 16, 60), 300))

    # ---------- tema ----------
    def apply_theme(self, mode: str):
        self.mode = mode
        if self.seq:
            self._refresh_selection()
