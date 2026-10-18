"""Motif ve alan taraması sayfası: konum şeritleri, motif tablosu, ayrıntı paneli."""
from PyQt5.QtCore import QRectF, Qt, QUrl, pyqtSignal
from PyQt5.QtGui import QColor, QDesktopServices, QFont, QPainter, QPen
from PyQt5.QtWidgets import (
    QAbstractItemView, QFrame, QHBoxLayout, QHeaderView, QLabel, QPushButton,
    QScrollArea, QTableWidget, QTableWidgetItem, QVBoxLayout, QWidget,
)

import theme

ROW_H = 28
LABEL_W = 84
AXIS_H = 24


class MotifTracks(QWidget):
    """Her motif için ayrı şerit; bulunan konumlar blok olarak çizilir."""
    clicked = pyqtSignal(int)  # eşleşme indeksi

    def __init__(self):
        super().__init__()
        self.setCursor(Qt.PointingHandCursor)
        self.n = 1
        self.ids = []
        self.matches = []
        self.sel = -1
        self.colors = theme.tokens("dark")

    def set_data(self, n, matches, sel, colors):
        self.n = max(n, 1)
        self.matches = matches
        self.ids = list(dict.fromkeys(m["id"] for m in matches))
        self.sel, self.colors = sel, colors
        self.setFixedHeight(max(len(self.ids), 1) * ROW_H + AXIS_H)
        self.update()

    def _track(self):
        return LABEL_W, max(self.width() - LABEL_W, 1)

    def _block(self, m):
        x0, w = self._track()
        row = self.ids.index(m["id"])
        x = x0 + (m["start"] - 1) / self.n * w
        bw = max((m["end"] - m["start"] + 1) / self.n * w, 4)
        return QRectF(x, row * ROW_H + 4, bw, ROW_H - 8)

    def paintEvent(self, _):
        c = self.colors
        p = QPainter(self)
        p.setRenderHint(QPainter.Antialiasing)
        x0, w = self._track()
        f = QFont(self.font())
        f.setPointSize(9)
        p.setFont(f)
        for i, mid in enumerate(self.ids):
            y = i * ROW_H
            p.setPen(QColor(c["muted"]))
            p.drawText(QRectF(0, y, LABEL_W - 8, ROW_H), Qt.AlignVCenter | Qt.AlignLeft, mid)
            p.setPen(Qt.NoPen)
            p.setBrush(QColor(c["bg"]))
            p.drawRoundedRect(QRectF(x0, y + 4, w, ROW_H - 8), 4, 4)
        for idx, m in enumerate(self.matches):
            p.setPen(Qt.NoPen)
            p.setBrush(QColor(c["accent"]))
            r = self._block(m)
            p.drawRoundedRect(r, 3, 3)
            if idx == self.sel:
                p.setPen(QPen(QColor(c["text"]), 2))
                p.setBrush(Qt.NoBrush)
                p.drawRoundedRect(r.adjusted(-1, -1, 1, 1), 3, 3)
        # eksen: 1, orta, son
        p.setPen(QColor(c["muted"]))
        y = len(self.ids) * ROW_H
        p.drawText(QRectF(x0, y, w, AXIS_H), Qt.AlignVCenter | Qt.AlignLeft, "1")
        p.drawText(QRectF(x0, y, w, AXIS_H), Qt.AlignVCenter | Qt.AlignHCenter, str(self.n // 2))
        p.drawText(QRectF(x0, y, w, AXIS_H), Qt.AlignVCenter | Qt.AlignRight, str(self.n))

    def mousePressEvent(self, e):
        for idx, m in enumerate(self.matches):
            if self._block(m).adjusted(-2, -2, 2, 2).contains(e.pos().x(), e.pos().y()):
                self.clicked.emit(idx)
                return


class MotifPage(QWidget):
    show_in_3d = pyqtSignal(int, int, str)  # başlangıç, bitiş, motif kimliği

    def __init__(self):
        super().__init__()
        self.mode = "dark"
        self.seq = ""
        self.matches = []
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
        title = QLabel("Motif ve alan taraması")
        title.setObjectName("pageTitle")
        self.badge = QLabel("")
        self.badge.setObjectName("badge")
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(self.badge)
        root.addLayout(head)

        root.addWidget(self._build_tracks_card())
        root.addWidget(self._build_table_card())
        root.addWidget(self._build_detail_card())
        root.addStretch(1)

    # ---------- bileşenler ----------
    def _build_tracks_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 10)
        lay.setSpacing(8)
        t = QLabel("Dizi üzerinde konumlar")
        t.setObjectName("cardTitle")
        lay.addWidget(t)
        self.tracks = MotifTracks()
        self.tracks.clicked.connect(self._select)
        lay.addWidget(self.tracks)
        return card

    def _build_table_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(12, 8, 12, 8)
        self.table = QTableWidget(0, 5)
        self.table.setHorizontalHeaderLabels(["Kimlik", "Ad", "Kalıp", "Konum", "Eşleşen dizi"])
        self.table.verticalHeader().hide()
        hh = self.table.horizontalHeader()
        hh.setHighlightSections(False)
        hh.setFixedHeight(34)
        hh.setDefaultAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        hh.setSectionResizeMode(QHeaderView.Stretch)
        hh.setSectionResizeMode(0, QHeaderView.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.ResizeToContents)
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

    def _build_detail_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(6)
        self.d_title = QLabel("")
        self.d_title.setObjectName("cardTitle")
        self.d_text = QLabel("")
        self.d_text.setObjectName("muted")
        self.d_text.setWordWrap(True)
        btns = QHBoxLayout()
        btns.setSpacing(8)
        self.btn_prosite = QPushButton("PROSITE kaydını aç")
        self.btn_3d = QPushButton("3B yapıda göster")
        for b in (self.btn_prosite, self.btn_3d):
            b.setObjectName("secondary")
            b.setCursor(Qt.PointingHandCursor)
            btns.addWidget(b)
        btns.addStretch(1)
        self.btn_prosite.clicked.connect(self._open_prosite)
        self.btn_3d.clicked.connect(self._emit_3d)
        lay.addWidget(self.d_title)
        lay.addWidget(self.d_text)
        lay.addLayout(btns)
        self.detail_card = card
        return card

    # ---------- veri ----------
    def set_result(self, info: dict):
        self.seq = info.get("seq", "")
        a = info["analysis"]
        self.matches = list(info.get("motifs", []))
        n_motifs = len({m["id"] for m in self.matches})
        self.badge.setText(f"{n_motifs} motif · {len(self.matches)} eşleşme")
        self._fill_table()
        self.sel = 0 if self.matches else -1
        self._refresh()

    def _fill_table(self):
        mono = QFont("Consolas")
        mono.setStyleHint(QFont.Monospace)
        self.table.blockSignals(True)
        self.table.setRowCount(len(self.matches))
        for r, m in enumerate(self.matches):
            vals = [m["id"], m["name"], m["pattern"], f"{m['start']}–{m['end']}", m["match"]]
            for c, v in enumerate(vals):
                it = QTableWidgetItem(v)
                it.setToolTip(v)
                if c in (2, 4):
                    it.setFont(mono)
                self.table.setItem(r, c, it)
        self.table.setFixedHeight(34 + 34 * max(len(self.matches), 1) + 4)
        self.table.blockSignals(False)

    # ---------- seçim ----------
    def _select(self, idx):
        self.sel = idx
        self._refresh()

    def _on_table_select(self):
        rows = self.table.selectionModel().selectedRows()
        if rows:
            self.sel = rows[0].row()
            self._refresh(from_table=True)

    def _refresh(self, from_table=False):
        c = theme.tokens(self.mode)
        self.tracks.set_data(len(self.seq), self.matches, self.sel, c)
        if not from_table and 0 <= self.sel < self.table.rowCount():
            self.table.blockSignals(True)
            self.table.selectRow(self.sel)
            self.table.blockSignals(False)
        has = 0 <= self.sel < len(self.matches)
        if has:
            m = self.matches[self.sel]
            self.d_title.setText(f"{m['id']} · {m['name']}")
            self.d_text.setText(m.get("description", ""))
        else:
            self.d_title.setText("Eşleşme bulunamadı")
            self.d_text.setText("Bu dizide bilinen bir motif ya da alan bulunamadı.")
        self.btn_prosite.setEnabled(has)
        self.btn_3d.setEnabled(has)

    # ---------- eylemler ----------
    def _open_prosite(self):
        if 0 <= self.sel < len(self.matches):
            QDesktopServices.openUrl(QUrl(f"https://prosite.expasy.org/{self.matches[self.sel]['id']}"))

    def _emit_3d(self):
        if 0 <= self.sel < len(self.matches):
            m = self.matches[self.sel]
            self.show_in_3d.emit(m["start"], m["end"], m["id"])

    # ---------- tema ----------
    def apply_theme(self, mode: str):
        self.mode = mode
        if self.seq:
            self._refresh()
