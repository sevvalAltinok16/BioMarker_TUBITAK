"""Raporlar sayfası: bölüm seçimi, rapor bilgileri, PDF rapor al, son raporlar."""
import os
from datetime import datetime

from PyQt5.QtCore import Qt, QUrl
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import (
    QCheckBox, QFileDialog, QFrame, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit,
    QPushButton, QScrollArea, QVBoxLayout, QWidget,
)

import analysis
from giris_page import _repolish

SECTIONS = [
    ("ozet", "Özet", "Protein bilgileri ve özet kartları"),
    ("fizikokimya", "Fizikokimyasal analiz", "Mol. ağırlık, pI, GRAVY, hidropati profili, titrasyon eğrisi"),
    ("ikincil", "İkincil yapı", "Chou-Fasman bandı, işaretli dizi, segment tablosu"),
    ("motif", "Motif ve alan taraması", "Eşleşen motifler ve konumları"),
    ("homoloji", "Homoloji taraması", "BLAST sonuçları ve hizalamalar"),
]
MAX_HISTORY = 5


class RaporlarPage(QWidget):
    def __init__(self):
        super().__init__()
        self.info = None
        self.history = []  # oturum boyunca alınan raporlar
        self._loading = False

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
        title = QLabel("Raporlar")
        title.setObjectName("pageTitle")
        self.badge = QLabel("")
        self.badge.setObjectName("badge")
        head.addWidget(title)
        head.addStretch(1)
        head.addWidget(self.badge)
        root.addLayout(head)

        self.empty = QLabel("Rapor almak için önce bir dizi analiz edin.")
        self.empty.setObjectName("muted")
        root.addWidget(self.empty)

        mid = QHBoxLayout()
        mid.setSpacing(10)
        mid.addWidget(self._build_sections_card(), 3)
        mid.addWidget(self._build_info_card(), 2)
        root.addLayout(mid)

        act = QHBoxLayout()
        self.msg = QLabel("")
        self.msg.setWordWrap(True)
        self.msg.hide()
        self.open_btn = QPushButton("Aç")
        self.open_btn.setObjectName("link")
        self.open_btn.setCursor(Qt.PointingHandCursor)
        self.open_btn.hide()
        self.open_btn.clicked.connect(self._open_last)
        self.run_btn = QPushButton("PDF rapor al (Ctrl+E)")
        self.run_btn.setObjectName("primary")
        self.run_btn.setCursor(Qt.PointingHandCursor)
        self.run_btn.clicked.connect(self.export_pdf)
        act.addWidget(self.msg, 1)
        act.addWidget(self.open_btn)
        act.addWidget(self.run_btn)
        root.addLayout(act)

        root.addWidget(self._build_history_card())
        root.addStretch(1)
        self._apply_state()

    # ---------- bileşenler ----------
    def _build_sections_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(10)
        t = QLabel("Rapor içeriği")
        t.setObjectName("cardTitle")
        lay.addWidget(t)
        self.checks, self.descs = {}, {}
        for key, label, desc in SECTIONS:
            col = QVBoxLayout()
            col.setSpacing(0)
            cb = QCheckBox(label)
            cb.setChecked(True)
            cb.toggled.connect(self._save_opts)
            d = QLabel(desc)
            d.setObjectName("muted")
            d.setContentsMargins(24, 0, 0, 0)
            col.addWidget(cb)
            col.addWidget(d)
            lay.addLayout(col)
            self.checks[key], self.descs[key] = cb, d
        return card

    def _build_info_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 14)
        lay.setSpacing(6)
        t = QLabel("Rapor bilgileri")
        t.setObjectName("cardTitle")
        l1 = QLabel("Başlık")
        l1.setObjectName("muted")
        self.title_edit = QLineEdit()
        self.title_edit.setObjectName("paramEdit")
        self.title_edit.textChanged.connect(self._save_opts)
        l2 = QLabel("Not (isteğe bağlı)")
        l2.setObjectName("muted")
        self.note_edit = QPlainTextEdit()
        self.note_edit.setObjectName("noteEdit")
        self.note_edit.setPlaceholderText("Rapora eklenecek kısa bir not yazın")
        self.note_edit.setMinimumHeight(90)
        self.note_edit.textChanged.connect(self._save_opts)
        for w in (t, l1, self.title_edit, l2, self.note_edit):
            lay.addWidget(w)
        lay.addStretch(1)
        return card

    def _build_history_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 12)
        lay.setSpacing(0)
        t = QLabel("Son raporlar")
        t.setObjectName("cardTitle")
        lay.addWidget(t)
        self.hist_box = QVBoxLayout()
        self.hist_box.setSpacing(0)
        lay.addLayout(self.hist_box)
        self._fill_history()
        return card

    # ---------- durum ----------
    def set_result(self, info):
        """info None ise bu sekmede analiz yok."""
        self.info = info
        self._clear_message()
        if info:
            self.badge.setText(f"{info['label']} · {info['analysis']['length']} aa")
            opts = info.setdefault("report_opts", {
                "sections": {k: True for k, _, _ in SECTIONS},
                "title": f"BioMarker raporu · {info['label']}",
                "note": "",
            })
            self._loading = True
            self.title_edit.setText(opts["title"])
            self.note_edit.setPlainText(opts["note"])
            self._loading = False
        else:
            self.badge.setText("")
        self._apply_state()

    def refresh(self):
        """Sayfa her gösterildiğinde bölümlerin uygunluğunu yeniden değerlendir."""
        self._apply_state()

    def _blast_hits(self):
        b = self.info.get("blast") if self.info else None
        return b["hits"] if b and b["state"] == "done" and b["hits"] else None

    def _apply_state(self):
        has = self.info is not None
        self.empty.setVisible(not has)
        self._loading = True
        opts = self.info.get("report_opts") if has else None
        hits = self._blast_hits() if has else None
        for key, _, desc in SECTIONS:
            cb = self.checks[key]
            avail = has and (key != "homoloji" or hits is not None)
            cb.setEnabled(avail)
            want = opts["sections"].get(key, True) if opts else False
            cb.setChecked(bool(want and avail))
            if key == "homoloji":
                if not has or hits is None:
                    self.descs[key].setText("Tarama yapılmadı. Homoloji sayfasından BLAST çalıştırın.")
                else:
                    added = len(self.info.get("report_hits") or [])
                    extra = f" · {added} sonuç rapora eklendi" if added else ""
                    self.descs[key].setText(f"BLAST sonuçları ve hizalamalar · {len(hits)} sonuç{extra}")
            else:
                self.descs[key].setText(desc)
        self.title_edit.setEnabled(has)
        self.note_edit.setEnabled(has)
        self.run_btn.setEnabled(has)
        self._loading = False

    def _save_opts(self, *_):
        if self._loading or not self.info:
            return
        opts = self.info.get("report_opts")
        if opts is None:
            return
        for key, cb in self.checks.items():
            if cb.isEnabled():  # uygun olmayan bölümün kullanıcı seçimini bozma
                opts["sections"][key] = cb.isChecked()
        opts["title"] = self.title_edit.text()
        opts["note"] = self.note_edit.toPlainText()

    # ---------- rapor alma ----------
    def selected_sections(self):
        return [k for k, _, _ in SECTIONS if self.checks[k].isEnabled() and self.checks[k].isChecked()]

    def export_pdf(self):
        """Düğme ve Ctrl+E bunu çağırır."""
        if not self.info:
            self._show_message("Rapor almak için önce bir dizi analiz edin.", ok=False)
            return
        if not self.selected_sections():
            self._show_message("Rapora en az bir bölüm ekleyin.", ok=False)
            return
        default = f"{self.info['accession']}_rapor.pdf".replace(" ", "_")
        path, _ = QFileDialog.getSaveFileName(self, "Raporu kaydet", default, "PDF (*.pdf)")
        if path:
            self.write_report(path)

    def write_report(self, path):
        if not path.lower().endswith(".pdf"):
            path += ".pdf"
        try:
            analysis.create_pdf_report(
                self.info, self.selected_sections(),
                self.title_edit.text().strip(), self.note_edit.toPlainText().strip(), path,
            )
        except analysis.AnalysisError as e:
            self._show_message(str(e), ok=False)
            return
        except Exception:
            self._show_message("Rapor oluşturulamadı. Dosya yolunu kontrol edip yeniden deneyin.", ok=False)
            return
        self._last_path = path
        self.history.insert(0, {"name": os.path.basename(path), "path": path,
                                "time": datetime.now().strftime("%d.%m.%Y %H:%M")})
        del self.history[MAX_HISTORY:]
        self._show_message(f"Rapor kaydedildi · {os.path.basename(path)}", ok=True)
        self._fill_history()

    # ---------- iletiler ve liste ----------
    def _show_message(self, text, ok):
        self.msg.setText(("✓ " if ok else "") + text)
        self.msg.setObjectName("okText" if ok else "errorText")
        _repolish(self.msg)
        self.msg.show()
        self.open_btn.setVisible(ok)

    def _clear_message(self):
        self.msg.hide()
        self.open_btn.hide()

    def _fill_history(self):
        while self.hist_box.count():
            wdg = self.hist_box.takeAt(0).widget()
            if wdg:
                wdg.setParent(None)
                wdg.deleteLater()
        if not self.history:
            none = QLabel("Henüz rapor oluşturulmadı.")
            none.setObjectName("muted")
            none.setContentsMargins(0, 8, 0, 4)
            self.hist_box.addWidget(none)
            return
        for h in self.history:
            sep = QFrame()
            sep.setObjectName("sep")
            sep.setFixedHeight(1)
            self.hist_box.addWidget(sep)
            row = QWidget()
            rl = QHBoxLayout(row)
            rl.setContentsMargins(0, 8, 0, 8)
            name = QLabel(h["name"])
            when = QLabel(h["time"])
            when.setObjectName("muted")
            b = QPushButton("Aç")
            b.setObjectName("link")
            b.setCursor(Qt.PointingHandCursor)
            b.clicked.connect(lambda _=False, p=h["path"]: QDesktopServices.openUrl(QUrl.fromLocalFile(p)))
            rl.addWidget(name)
            rl.addStretch(1)
            rl.addWidget(when)
            rl.addSpacing(10)
            rl.addWidget(b)
            self.hist_box.addWidget(row)

    def _open_last(self):
        p = getattr(self, "_last_path", None)
        if p:
            QDesktopServices.openUrl(QUrl.fromLocalFile(p))

    def apply_theme(self, mode: str):
        pass  # renkler QSS'ten geliyor
