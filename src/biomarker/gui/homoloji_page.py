"""Homoloji taraması sayfası: parametreler, çalışma durumu, sonuç tablosu, hizalama."""
import threading
import time

from PyQt5.QtCore import QObject, Qt, QTimer, QUrl, pyqtSignal
from PyQt5.QtGui import QDesktopServices, QFont
from PyQt5.QtWidgets import (
    QAbstractItemView, QComboBox, QFrame, QHBoxLayout, QHeaderView, QLabel,
    QLineEdit, QProgressBar, QPushButton, QScrollArea, QSpinBox, QTableWidget,
    QTableWidgetItem, QVBoxLayout, QWidget,
)

import analysis
import theme
from giris_page import _repolish

PAGE = 10
DATABASES = ["Swiss-Prot", "nr", "PDB"]
PROGRAMS = ["blastp"]
LINE_W = 60  # hizalama satırında karakter sayısı


def mmss(sec: float) -> str:
    s = int(sec)
    return f"{s // 60:02d}:{s % 60:02d}"


def fmt_ev(x: float) -> str:
    return f"{x:.0e}" if x else "0"


class _Bridge(QObject):
    """Arka plan iş parçacığından arayüze sonuç taşır."""
    done = pyqtSignal(int, list)
    failed = pyqtSignal(int, str)


class HomolojiPage(QWidget):
    state_changed = pyqtSignal(str)  # "running" | "done" | "error" | "idle"

    def __init__(self):
        super().__init__()
        self.info = None
        self.shown = PAGE
        self.sel = -1
        self._rid = 0
        self._runs = {}
        self.bridge = _Bridge()
        self.bridge.done.connect(self._on_done)
        self.bridge.failed.connect(self._on_failed)
        self.timer = QTimer(self)
        self.timer.setInterval(1000)
        self.timer.timeout.connect(self._tick)

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
        title = QLabel("Homoloji taraması")
        title.setObjectName("pageTitle")
        self.badge = QLabel("")
        self.badge.setObjectName("badge")
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(self.badge)
        root.addLayout(head)

        root.addWidget(self._build_params_card())
        root.addWidget(self._build_warn_box())
        root.addLayout(self._build_status_row())
        root.addWidget(self._build_error_box())
        root.addWidget(self._build_empty_box())
        root.addWidget(self._build_results_card())
        root.addWidget(self._build_align_card())
        root.addStretch(1)
        self._render()

    # ---------- bileşenler ----------
    def _build_params_card(self):
        card = QFrame()
        card.setObjectName("card")
        row = QHBoxLayout(card)
        row.setContentsMargins(16, 12, 16, 12)
        row.setSpacing(10)

        def field(label, widget):
            col = QVBoxLayout()
            col.setSpacing(4)
            lab = QLabel(label)
            lab.setObjectName("muted")
            col.addWidget(lab)
            col.addWidget(widget)
            row.addLayout(col, 1)

        self.db_box = QComboBox()
        self.db_box.addItems(DATABASES)
        self.prog_box = QComboBox()
        self.prog_box.addItems(PROGRAMS)
        self.ev_edit = QLineEdit("0,001")
        self.ev_edit.setObjectName("paramEdit")
        self.hits_spin = QSpinBox()
        self.hits_spin.setRange(1, 500)
        self.hits_spin.setValue(50)
        field("Veritabanı", self.db_box)
        field("Program", self.prog_box)
        field("E-değeri eşiği", self.ev_edit)
        field("En çok sonuç", self.hits_spin)

        col = QVBoxLayout()
        col.addStretch(1)
        self.run_btn = QPushButton("BLAST çalıştır")
        self.run_btn.setObjectName("primary")
        self.run_btn.setCursor(Qt.PointingHandCursor)
        self.run_btn.setFixedWidth(150)
        self.run_btn.clicked.connect(self._on_run_clicked)
        col.addWidget(self.run_btn)
        row.addLayout(col)
        return card

    def _build_warn_box(self):
        self.warn = QFrame()
        self.warn.setObjectName("warnbox")
        lay = QHBoxLayout(self.warn)
        lay.setContentsMargins(14, 10, 14, 10)
        lab = QLabel("İnternet bağlantısı gerekiyor. Diğer modüller çevrimdışı çalışır.")
        lab.setObjectName("warnText")
        lay.addWidget(lab)
        return self.warn

    def _build_status_row(self):
        row = QHBoxLayout()
        row.setSpacing(10)
        self.st_icon = QLabel("✓")
        self.st_icon.setObjectName("okText")
        self.st_text = QLabel("")
        self.bar = QProgressBar()
        self.bar.setRange(0, 0)
        self.bar.setTextVisible(False)
        self.bar.setFixedWidth(120)
        self.st_time = QLabel("")
        self.st_time.setObjectName("muted")
        row.addWidget(self.st_icon)
        row.addWidget(self.st_text)
        row.addWidget(self.bar)
        row.addStretch(1)
        row.addWidget(self.st_time)
        self.status_widgets = [self.st_icon, self.st_text, self.bar, self.st_time]
        return row

    def _build_error_box(self):
        self.err = QFrame()
        self.err.setObjectName("errorbox")
        lay = QHBoxLayout(self.err)
        lay.setContentsMargins(14, 10, 14, 10)
        self.err_text = QLabel("")
        self.err_text.setObjectName("errorText")
        self.err_text.setWordWrap(True)
        self.err_retry = QPushButton("Yeniden dene")
        self.err_retry.setObjectName("secondary")
        self.err_retry.setCursor(Qt.PointingHandCursor)
        self.err_retry.clicked.connect(self._start)
        lay.addWidget(self.err_text, 1)
        lay.addWidget(self.err_retry)
        return self.err

    def _build_empty_box(self):
        self.empty = QFrame()
        self.empty.setObjectName("emptybox")
        self.empty.setMinimumHeight(90)
        lay = QVBoxLayout(self.empty)
        self.empty_text = QLabel("Sonuçlar tamamlanınca burada görünecek.")
        self.empty_text.setObjectName("muted")
        self.empty_text.setAlignment(Qt.AlignCenter)
        lay.addWidget(self.empty_text)
        return self.empty

    def _build_results_card(self):
        self.results = QFrame()
        self.results.setObjectName("card")
        lay = QVBoxLayout(self.results)
        lay.setContentsMargins(12, 8, 12, 10)
        lay.setSpacing(6)
        self.table = QTableWidget(0, 6)
        self.table.setHorizontalHeaderLabels(["Erişim no", "Tür", "Skor ↓", "E-değeri", "Özd. %", "Kapsam %"])
        self.table.verticalHeader().hide()
        hh = self.table.horizontalHeader()
        hh.setHighlightSections(False)
        hh.setFixedHeight(34)
        hh.setDefaultAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        hh.setSectionResizeMode(QHeaderView.Stretch)
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
        foot = QHBoxLayout()
        self.note = QLabel("")
        self.note.setObjectName("muted")
        self.more_btn = QPushButton("Daha fazla göster")
        self.more_btn.setObjectName("link")
        self.more_btn.setCursor(Qt.PointingHandCursor)
        self.more_btn.clicked.connect(self._show_more)
        foot.addWidget(self.note)
        foot.addStretch(1)
        foot.addWidget(self.more_btn)
        lay.addLayout(foot)
        return self.results

    def _build_align_card(self):
        self.align = QFrame()
        self.align.setObjectName("card")
        lay = QVBoxLayout(self.align)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(8)
        head = QHBoxLayout()
        self.a_title = QLabel("")
        self.a_title.setObjectName("cardTitle")
        self.a_range = QLabel("")
        self.a_range.setObjectName("muted")
        head.addWidget(self.a_title)
        head.addStretch(1)
        head.addWidget(self.a_range)
        lay.addLayout(head)
        self.a_text = QLabel("")
        mono = QFont("Consolas")
        mono.setStyleHint(QFont.Monospace)
        self.a_text.setFont(mono)
        self.a_text.setTextInteractionFlags(Qt.TextSelectableByMouse)
        lay.addWidget(self.a_text)
        btns = QHBoxLayout()
        btns.setSpacing(8)
        self.btn_uniprot = QPushButton("UniProt'ta aç")
        self.btn_report = QPushButton("Rapora ekle")
        for b in (self.btn_uniprot, self.btn_report):
            b.setObjectName("secondary")
            b.setCursor(Qt.PointingHandCursor)
            b.setMinimumWidth(140)
            btns.addWidget(b)
        btns.addStretch(1)
        self.btn_uniprot.clicked.connect(self._open_uniprot)
        self.btn_report.clicked.connect(self._add_report)
        lay.addLayout(btns)
        return self.align

    # ---------- veri ----------
    def set_result(self, info: dict):
        self.info = info
        self.shown = PAGE
        self.sel = 0
        self.badge.setText(f"{info['label']} · {info['analysis']['length']} aa")
        b = info.get("blast")
        if b and b["state"] == "running":
            self.state_changed.emit("running")
        self._render()

    def _blast(self):
        return self.info.get("blast") if self.info else None

    # ---------- çalıştırma ----------
    def _on_run_clicked(self):
        b = self._blast()
        if b and b["state"] == "running":
            self._cancel()
        else:
            self._start()

    def _start(self):
        if not self.info:
            return
        try:
            evalue = float(self.ev_edit.text().strip().replace(",", "."))
            if evalue <= 0:
                raise ValueError
        except ValueError:
            self._show_input_error("E-değeri pozitif bir sayı olmalı. Örnek: 0,001")
            return
        if not analysis.is_online():
            self._render()
            return
        params = {
            "database": self.db_box.currentText(), "program": self.prog_box.currentText(),
            "evalue": evalue, "max_hits": self.hits_spin.value(),
        }
        self._rid += 1
        rid = self._rid
        b = {"state": "running", "started": time.time(), "elapsed": 0.0,
             "hits": [], "error": "", "rid": rid}
        self.info["blast"] = b
        self._runs[rid] = b
        self.sel, self.shown = 0, PAGE
        seq = self.info["seq"]

        def work():
            try:
                self.bridge.done.emit(rid, analysis.run_blast(seq, **params))
            except TimeoutError:
                self.bridge.failed.emit(rid, "BLAST sunucusu yanıt vermedi. Yeniden deneyin.")
            except analysis.AnalysisError as e:
                self.bridge.failed.emit(rid, str(e))
            except Exception:
                self.bridge.failed.emit(rid, "BLAST isteği başarısız oldu. Bağlantıyı kontrol edip yeniden deneyin.")

        threading.Thread(target=work, daemon=True).start()
        self.state_changed.emit("running")
        self._render()

    def _cancel(self):
        b = self._blast()
        if b:
            self._runs.pop(b["rid"], None)
        self.info["blast"] = None
        self.state_changed.emit("idle")
        self._render()

    def _on_done(self, rid, hits):
        b = self._runs.pop(rid, None)
        if b is None:
            return  # iptal edilmiş
        b["state"] = "done"
        b["elapsed"] = time.time() - b["started"]
        b["hits"] = sorted(hits, key=lambda h: -h["score"])
        self.state_changed.emit("done")
        if self.info and self.info.get("blast") is b:
            self.sel = 0
            self._render()

    def _on_failed(self, rid, message):
        b = self._runs.pop(rid, None)
        if b is None:
            return
        b["state"], b["error"] = "error", message
        b["elapsed"] = time.time() - b["started"]
        self.state_changed.emit("error")
        if self.info and self.info.get("blast") is b:
            self._render()

    def _show_input_error(self, text):
        self.info["blast"] = {"state": "error", "started": 0, "elapsed": 0.0, "hits": [],
                              "error": text, "rid": -1, "retry": False}
        self.state_changed.emit("error")
        self._render()

    def _tick(self):
        b = self._blast()
        if b and b["state"] == "running":
            self.st_time.setText(f"Geçen süre {mmss(time.time() - b['started'])}")

    # ---------- ekranı çiz ----------
    def _render(self):
        b = self._blast()
        state = b["state"] if b else "idle"
        online = analysis.is_online()
        running = state == "running"

        self.warn.setVisible(not online)
        self.run_btn.setText("İptal" if running else "BLAST çalıştır")
        self.run_btn.setObjectName("secondary" if running else "primary")
        _repolish(self.run_btn)
        self.run_btn.setEnabled(running or (online and self.info is not None))

        for w in self.status_widgets:
            w.setVisible(False)
        if running:
            self.st_text.setText("İşlem devam ediyor…")
            self.st_time.setText(f"Geçen süre {mmss(time.time() - b['started'])}")
            for w in (self.st_text, self.bar, self.st_time):
                w.setVisible(True)
            self.timer.start()
        else:
            self.timer.stop()
        if state == "done":
            hits = b["hits"]
            self.st_text.setText(f"Tamamlandı · {len(hits)} sonuç · {mmss(b['elapsed'])}")
            self.st_icon.setVisible(True)
            self.st_text.setVisible(True)

        self.err.setVisible(state == "error")
        if state == "error":
            self.err_text.setText(b["error"])
            self.err_retry.setVisible(b.get("retry", True))

        has_hits = state == "done" and bool(b["hits"])
        self.results.setVisible(has_hits)
        self.align.setVisible(has_hits)
        self.empty.setVisible(not has_hits)
        self.empty_text.setText(
            "Bu dizi için sonuç bulunamadı." if state == "done" else "Sonuçlar tamamlanınca burada görünecek."
        )
        if has_hits:
            self._fill_table(b["hits"])

    def _fill_table(self, hits):
        n = min(self.shown, len(hits))
        self.table.blockSignals(True)
        self.table.setRowCount(n)
        for r in range(n):
            h = hits[r]
            vals = [h["accession"], h["species"], h["score"], fmt_ev(h["evalue"]), h["identity"], h["coverage"]]
            for c, v in enumerate(vals):
                self.table.setItem(r, c, QTableWidgetItem(str(v)))
        self.table.setFixedHeight(34 + 34 * n + 4)
        self.sel = min(max(self.sel, 0), n - 1)
        self.table.selectRow(self.sel)
        self.table.blockSignals(False)
        self.note.setText(f"{n} / {len(hits)} sonuç gösteriliyor")
        self.more_btn.setVisible(n < len(hits))
        self._fill_alignment(hits[self.sel])

    def _show_more(self):
        self.shown += PAGE
        self._render()

    def _on_table_select(self):
        rows = self.table.selectionModel().selectedRows()
        b = self._blast()
        if rows and b and b["state"] == "done":
            self.sel = rows[0].row()
            self._fill_alignment(b["hits"][self.sel])

    def _fill_alignment(self, h):
        a = h["alignment"]
        self.a_title.setText(f"Hizalama · {h['accession']}")
        self.a_range.setText(f"Sorgu {a['query_start']}–{a['query_end']}")
        q, m, s = a["query"], a["match"], a["subject"]
        qp, sp = a["query_start"], a["subject_start"]
        lines = []
        for i in range(0, len(q), LINE_W):
            qc, mc, sc = q[i:i + LINE_W], m[i:i + LINE_W], s[i:i + LINE_W]
            qn, sn = sum(ch != "-" for ch in qc), sum(ch != "-" for ch in sc)
            lines.append(f"Sorgu {qp:>5}  {qc}  {qp + qn - 1}")
            lines.append(f"      {'':>5}  {mc}")
            lines.append(f"Hedef {sp:>5}  {sc}  {sp + sn - 1}")
            lines.append("")
            qp, sp = qp + qn, sp + sn
        self.a_text.setText("\n".join(lines).rstrip())
        done = h["accession"] in (self.info.get("report_hits") or [])
        self.btn_report.setText("Rapora eklendi" if done else "Rapora ekle")
        self.btn_report.setEnabled(not done)

    # ---------- eylemler ----------
    def _current_hit(self):
        b = self._blast()
        if b and b["state"] == "done" and 0 <= self.sel < len(b["hits"]):
            return b["hits"][self.sel]
        return None

    def _open_uniprot(self):
        h = self._current_hit()
        if h:
            QDesktopServices.openUrl(QUrl(f"https://www.uniprot.org/uniprotkb/{h['accession']}"))

    def _add_report(self):
        h = self._current_hit()
        if h:
            self.info.setdefault("report_hits", []).append(h["accession"])
            self._fill_alignment(h)

    def refresh_connection(self):
        """Sayfa her gösterildiğinde bağlantı durumunu yeniden kontrol et."""
        self._render()

    def apply_theme(self, mode: str):
        pass  # renkler QSS'ten geliyor
