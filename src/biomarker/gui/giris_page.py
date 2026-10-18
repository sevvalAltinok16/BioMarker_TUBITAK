"""Giriş ekranı: dizi yapıştır / dosya / kimlik girişi."""
import os
import re

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QFont, QTextCharFormat, QTextCursor
from PyQt5.QtWidgets import (
    QButtonGroup, QFileDialog, QFrame, QHBoxLayout, QLabel, QLineEdit,
    QPlainTextEdit, QPushButton, QStackedWidget, QTextEdit, QVBoxLayout, QWidget,
)

import theme

# Geçerli harfler (20 standart + B, X, Z, U, O). İleride karar değişirse sadece burası.
VALID = set("ACDEFGHIKLMNPQRSTVWYBXZUO")

UNIPROT_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
NCBI_RE = re.compile(r"^[A-Z]{2}_\d+(\.\d+)?$")

EXAMPLES = {
    "İnsülin": (
        ">sp|P01308|INS_HUMAN",
        "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
    ),
    "Hemoglobin": (
        ">sp|P69905|HBA_HUMAN",
        "MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR",
    ),
    "Miyoglobin": (
        ">sp|P02144|MYG_HUMAN",
        "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKTEAEMKASEDLKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHPGDFGADAQGAMNKALELFRKDIAAKYKELGYQG",
    ),
}

# Örnek veri: gerçek geçmiş sistemi sonra bağlanacak
RECENT = [
    ("P69905 · Hemoglobin α", "dün"),
    ("P01308 · İnsülin", "2 gün önce"),
    ("Özel dizi (86 aa)", "geçen hafta"),
]


def _repolish(w):
    w.style().unpolish(w)
    w.style().polish(w)


class GirisPage(QWidget):
    analyze_requested = pyqtSignal(str, str)  # tür ("dizi" | "kimlik"), değer

    def __init__(self):
        super().__init__()
        root = QVBoxLayout(self)
        root.setContentsMargins(20, 16, 20, 16)
        root.setSpacing(12)

        title = QLabel("Yeni analiz")
        title.setObjectName("pageTitle")
        sub = QLabel("Bir protein dizisi girin, UniProt/NCBI kimliği yazın ya da dosya bırakın.")
        sub.setObjectName("muted")
        root.addWidget(title)
        root.addWidget(sub)

        root.addWidget(self._build_input_card())
        root.addWidget(self._build_recent_card())
        root.addStretch(1)

    # ---------- girdi kartı ----------
    def _build_input_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(10)

        # sekmeler (hap düğmeler)
        pills = QHBoxLayout()
        pills.setSpacing(6)
        self.pill_group = QButtonGroup(self)
        self.pill_group.setExclusive(True)
        for i, text in enumerate(("Dizi yapıştır", "Dosya", "UniProt / NCBI kimliği")):
            b = QPushButton(text)
            b.setObjectName("pill")
            b.setCheckable(True)
            b.setCursor(Qt.PointingHandCursor)
            self.pill_group.addButton(b, i)
            pills.addWidget(b)
        pills.addStretch(1)
        self.pill_group.button(0).setChecked(True)
        self.pill_group.buttonClicked[int].connect(self._on_pill)
        lay.addLayout(pills)

        # sayfalar
        self.stack = QStackedWidget()
        self.stack.setObjectName("inputStack")

        self.seq_edit = QPlainTextEdit()
        self.seq_edit.setObjectName("seqEdit")
        self.seq_edit.setPlaceholderText("FASTA dizisini buraya yapıştırın")
        self.seq_edit.setMinimumHeight(120)
        self.seq_edit.setAcceptDrops(False)  # dosya bırakma pencereye gitsin
        mono = QFont("Consolas")
        mono.setStyleHint(QFont.Monospace)
        self.seq_edit.setFont(mono)
        self.seq_edit.textChanged.connect(self._check_sequence)
        self.stack.addWidget(self.seq_edit)

        self.dropzone = QFrame()
        self.dropzone.setObjectName("dropzone")
        self.dropzone.setMinimumHeight(120)
        dz = QVBoxLayout(self.dropzone)
        dz.setAlignment(Qt.AlignCenter)
        self.drop_label = QLabel(".fasta, .txt veya .pdb dosyasını buraya bırakın")
        self.drop_label.setAlignment(Qt.AlignCenter)
        pick = QPushButton("Dosya seç")
        pick.setObjectName("secondary")
        pick.setCursor(Qt.PointingHandCursor)
        pick.clicked.connect(self.open_file_dialog)
        dz.addWidget(self.drop_label)
        dz.addSpacing(6)
        dz.addWidget(pick, 0, Qt.AlignCenter)
        self.stack.addWidget(self.dropzone)

        id_page = QWidget()
        il = QVBoxLayout(id_page)
        il.setContentsMargins(0, 0, 0, 0)
        self.id_edit = QLineEdit()
        self.id_edit.setObjectName("idEdit")
        self.id_edit.setPlaceholderText("P01308 veya NP_000198")
        self.id_edit.setAcceptDrops(False)
        self.id_edit.returnPressed.connect(self.submit)
        il.addWidget(self.id_edit)
        il.addStretch(1)
        self.stack.addWidget(id_page)

        lay.addWidget(self.stack)

        self.msg = QLabel("")
        self.msg.setObjectName("error")
        self.msg.setWordWrap(True)
        self.msg.hide()
        lay.addWidget(self.msg)

        # alt satır: örnekler + Analiz et
        bottom = QHBoxLayout()
        bottom.setSpacing(6)
        ex = QLabel("Örnekler:")
        ex.setObjectName("muted")
        bottom.addWidget(ex)
        for name in EXAMPLES:
            c = QPushButton(name)
            c.setObjectName("chip")
            c.setCursor(Qt.PointingHandCursor)
            c.clicked.connect(lambda _=False, n=name: self.load_example(n))
            bottom.addWidget(c)
        bottom.addStretch(1)
        self.run_btn = QPushButton("Analiz et (F5)")
        self.run_btn.setObjectName("primary")
        self.run_btn.setCursor(Qt.PointingHandCursor)
        self.run_btn.clicked.connect(self.submit)
        bottom.addWidget(self.run_btn)
        lay.addLayout(bottom)
        return card

    # ---------- son analizler ----------
    def _build_recent_card(self):
        card = QFrame()
        card.setObjectName("card")
        lay = QVBoxLayout(card)
        lay.setContentsMargins(16, 12, 16, 8)
        lay.setSpacing(0)
        head = QLabel("Son analizler")
        head.setObjectName("cardTitle")
        lay.addWidget(head)
        lay.addSpacing(6)
        for i, (name, when) in enumerate(RECENT):
            sep = QFrame()
            sep.setObjectName("sep")
            sep.setFixedHeight(1)
            lay.addWidget(sep)
            row = QHBoxLayout()
            row.setContentsMargins(0, 8, 0, 8)
            n = QLabel(name)
            w = QLabel(when)
            w.setObjectName("muted")
            row.addWidget(n)
            row.addStretch(1)
            row.addWidget(w)
            lay.addLayout(row)
        return card

    # ---------- davranış ----------
    def _on_pill(self, index):
        self.stack.setCurrentIndex(index)
        self._clear_message()

    def load_example(self, name):
        header, seq = EXAMPLES[name]
        self.pill_group.button(0).setChecked(True)
        self.stack.setCurrentIndex(0)
        self.seq_edit.setPlainText(f"{header}\n{seq}")

    def open_file_dialog(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Dosya seç", "", "Dizi dosyaları (*.fasta *.fa *.txt *.pdb)"
        )
        if path:
            self.load_file(path)

    def load_file(self, path):
        ext = os.path.splitext(path)[1].lower()
        if ext not in (".fasta", ".fa", ".txt", ".pdb"):
            self._show_message("Desteklenmeyen dosya türü. .fasta, .txt veya .pdb seçin.")
            return
        if ext == ".pdb":
            self.pill_group.button(1).setChecked(True)
            self.stack.setCurrentIndex(1)
            self.drop_label.setText(os.path.basename(path))
            self._show_message("PDB dosyasından dizi okuma henüz eklenmedi.")
            return
        try:
            with open(path, "r", encoding="utf-8", errors="replace") as f:
                text = f.read()
        except OSError:
            self._show_message("Dosya okunamadı. Dosyayı kontrol edip tekrar deneyin.")
            return
        self.pill_group.button(0).setChecked(True)
        self.stack.setCurrentIndex(0)
        self.seq_edit.setPlainText(text.strip())

    def show_error(self, text):
        self._show_message(text)

    def clear(self):
        self.seq_edit.clear()
        self.id_edit.clear()
        self.pill_group.button(0).setChecked(True)
        self.stack.setCurrentIndex(0)
        self.drop_label.setText(".fasta, .txt veya .pdb dosyasını buraya bırakın")
        self._clear_message()

    def set_drag(self, on: bool):
        self.dropzone.setProperty("drag", on)
        _repolish(self.dropzone)

    def _clear_message(self):
        self.msg.hide()

    def _show_message(self, text):
        self.msg.setText(text)
        self.msg.show()

    def _check_sequence(self):
        """IUPAC dışı karakterlere dalgalı kırmızı alt çizgi."""
        doc = self.seq_edit.document()
        fmt = QTextCharFormat()
        fmt.setUnderlineStyle(QTextCharFormat.WaveUnderline)
        fmt.setUnderlineColor(QColor(theme.COMMON["error"]))
        sels, invalid, pos = [], 0, 0
        for line in self.seq_edit.toPlainText().split("\n"):
            if not line.startswith(">"):
                for i, ch in enumerate(line):
                    if ch.isspace() or ch.upper() in VALID:
                        continue
                    invalid += 1
                    if len(sels) < 500:
                        cur = QTextCursor(doc)
                        cur.setPosition(pos + i)
                        cur.setPosition(pos + i + 1, QTextCursor.KeepAnchor)
                        s = QTextEdit.ExtraSelection()
                        s.cursor = cur
                        s.format = fmt
                        sels.append(s)
            pos += len(line) + 1
        self.seq_edit.setExtraSelections(sels)
        if invalid:
            self._show_message("Dizide geçersiz karakterler var. Vurgulanan konumları düzeltin.")
        else:
            self._clear_message()
        return invalid == 0

    def submit(self):
        mode = self.stack.currentIndex()
        if mode == 0:
            text = self.seq_edit.toPlainText().strip()
            if not text:
                self._show_message("Önce bir dizi girin ya da örnek seçin.")
                return
            if not self._check_sequence():
                return
            self.analyze_requested.emit("dizi", text)
        elif mode == 1:
            self._show_message("Önce bir dosya seçin ya da buraya bırakın.")
        else:
            value = self.id_edit.text().strip().upper()
            if not (UNIPROT_RE.match(value) or NCBI_RE.match(value)):
                self._show_message("Kimlik biçimi tanınmadı. Örnek: P01308 veya NP_000198.")
                return
            self._clear_message()
            self.analyze_requested.emit("kimlik", value)