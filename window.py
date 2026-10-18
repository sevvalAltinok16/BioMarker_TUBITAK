from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QKeySequence
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QTabBar, QToolButton,
    QFrame, QPushButton, QButtonGroup, QStackedWidget, QLabel, QAction,
    QActionGroup, QApplication,
)

import theme
from giris_page import GirisPage

# (anahtar, etiket, analiz gerektirir mi)
NAV_ITEMS = [
    ("giris", "Giriş", False),
    ("gecmis", "Geçmiş", False),
    ("fizikokimya", "Fizikokimyasal", True),
    ("ikincil", "İkincil yapı", True),
    ("motif", "Motif", True),
    ("homoloji", "Homoloji", True),
    ("yapi3b", "3B yapı", False),
    ("raporlar", "Raporlar", False),
]


class PlaceholderPage(QWidget):
    """Sayfalar yazılana kadar yer tutucu."""

    def __init__(self, title: str):
        super().__init__()
        lay = QVBoxLayout(self)
        lay.setContentsMargins(20, 16, 20, 16)
        head = QLabel(title)
        head.setObjectName("pageTitle")
        note = QLabel("Bu sayfa henüz yapılmadı.")
        note.setObjectName("muted")
        lay.addWidget(head)
        lay.addWidget(note)
        lay.addStretch(1)


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("BioMarker")
        self.setMinimumSize(1100, 700)
        self.resize(1280, 800)
        self.setAcceptDrops(True)
        self.theme_mode = "dark"
        self.theme_choice = "system"

        self._build_menu()
        self._build_body()
        self._build_status()
        self.set_theme("system")
        self.new_tab()

    # ---------- menü ----------
    def _build_menu(self):
        mb = self.menuBar()

        m_file = mb.addMenu("Dosya")
        self.act_new = self._action(m_file, "Yeni analiz", "Ctrl+N", self.new_tab)
        self._action(m_file, "FASTA içe aktar", "Ctrl+O", lambda: self.giris.open_file_dialog())
        self._action(m_file, "Projeyi kaydet", "Ctrl+S")
        self._action(m_file, "Raporu dışa aktar", "Ctrl+E")

        m_edit = mb.addMenu("Düzen")
        self._action(m_edit, "Diziyi kopyala")

        m_view = mb.addMenu("Görünüm")
        grp = QActionGroup(self)
        for key, label in (("light", "Aydınlık"), ("dark", "Koyu"), ("system", "Sistem")):
            a = QAction(label, self, checkable=True)
            a.setChecked(key == "system")
            a.triggered.connect(lambda _=False, k=key: self.set_theme(k))
            grp.addAction(a)
            m_view.addAction(a)

        m_tools = mb.addMenu("Araçlar")
        self._action(m_tools, "Dil seçimi")

        m_help = mb.addMenu("Yardım")
        self._action(m_help, "Belgeler")
        self._action(m_help, "Hakkında")

        # F5: Analiz et (şimdilik bağlı değil)
        self.act_run = QAction("Analiz et", self)
        self.act_run.setShortcut(QKeySequence("F5"))
        self.act_run.triggered.connect(lambda _=False: self.giris.submit())
        self.addAction(self.act_run)

    def _action(self, menu, text, shortcut=None, slot=None):
        a = QAction(text, self)
        if shortcut:
            a.setShortcut(QKeySequence(shortcut))
        if slot:
            a.triggered.connect(lambda _=False: slot())
        menu.addAction(a)
        return a

    # ---------- gövde ----------
    def _build_body(self):
        root = QWidget()
        root.setObjectName("workspace")
        col = QVBoxLayout(root)
        col.setContentsMargins(0, 0, 0, 0)
        col.setSpacing(0)

        # sekme satırı
        row = QFrame()
        row.setObjectName("tabRow")
        rl = QHBoxLayout(row)
        rl.setContentsMargins(8, 6, 8, 0)
        rl.setSpacing(0)
        self.tabbar = QTabBar()
        self.tabbar.setTabsClosable(True)
        self.tabbar.setExpanding(False)
        self.tabbar.setDrawBase(False)
        self.tabbar.tabCloseRequested.connect(self.close_tab)
        plus = QToolButton()
        plus.setObjectName("newTab")
        plus.setText("+")
        plus.clicked.connect(self.new_tab)
        rl.addWidget(self.tabbar)
        rl.addWidget(plus)
        rl.addStretch(1)
        col.addWidget(row)

        # sidebar + sayfalar
        body = QHBoxLayout()
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(0)
        body.addWidget(self._build_sidebar())
        self.stack = QStackedWidget()
        for key, label, _ in NAV_ITEMS:
            if key == "giris":
                self.giris = GirisPage()
                self.giris.analyze_requested.connect(self.on_analyze)
                self.stack.addWidget(self.giris)
            else:
                self.stack.addWidget(PlaceholderPage(label))
        self.stack.addWidget(PlaceholderPage("Ayarlar"))
        body.addWidget(self.stack, 1)
        col.addLayout(body, 1)

        self.setCentralWidget(root)

    def _build_sidebar(self):
        side = QFrame()
        side.setObjectName("sidebar")
        side.setFixedWidth(200)
        lay = QVBoxLayout(side)
        lay.setContentsMargins(8, 0, 8, 10)
        lay.setSpacing(2)

        brand = QLabel("BioMarker")
        brand.setObjectName("brand")
        lay.addWidget(brand)

        self.nav_group = QButtonGroup(self)
        self.nav_group.setExclusive(True)
        self.nav_buttons = {}
        for i, (key, label, needs_analysis) in enumerate(NAV_ITEMS):
            lay.addWidget(self._nav_button(i, key, label, needs_analysis))
        lay.addStretch(1)
        lay.addWidget(self._nav_button(len(NAV_ITEMS), "ayarlar", "Ayarlar", False))

        self.nav_buttons["giris"].setChecked(True)
        return side

    def _nav_button(self, index, key, label, needs_analysis):
        b = QPushButton(label)
        b.setObjectName("navBtn")
        b.setCheckable(True)
        b.setCursor(Qt.PointingHandCursor)
        b.setEnabled(not needs_analysis)  # analiz yapılana kadar soluk
        self.nav_group.addButton(b, index)
        self.nav_buttons[key] = b
        b.clicked.connect(lambda _=False, i=index: self.stack.setCurrentIndex(i))
        return b

    # ---------- durum çubuğu ----------
    def _build_status(self):
        sb = self.statusBar()
        sb.setSizeGripEnabled(False)
        self.status_left = QLabel("Durum: Hazır")
        self.status_right = QLabel("Çevrimdışı analiz: açık")
        sb.addWidget(self.status_left, 1)
        sb.addPermanentWidget(self.status_right)

    # ---------- sekmeler ----------
    def new_tab(self):
        i = self.tabbar.addTab("Yeni analiz")
        self.tabbar.setCurrentIndex(i)

    def close_tab(self, index):
        if self.tabbar.count() > 1:
            self.tabbar.removeTab(index)

    # ---------- tema ----------
    def set_theme(self, choice: str):
        self.theme_choice = choice
        if choice == "system":
            pal = QApplication.instance().palette()
            mode = "dark" if pal.window().color().lightness() < 128 else "light"
        else:
            mode = choice
        self.theme_mode = mode
        QApplication.instance().setStyleSheet(theme.build_qss(mode))
    
        # ---------- analiz (şimdilik sahte, API sonra bağlanacak) ----------
    def on_analyze(self, kind, value):
        self.status_left.setText("Durum: Analiz ediliyor…")
        QTimer.singleShot(1200, lambda: self.status_left.setText("Durum: Hazır"))

    # ---------- sürükle-bırak: pencerenin her yerine ----------
    @staticmethod
    def _first_local_file(event):
        for url in event.mimeData().urls():
            if url.isLocalFile():
                return url.toLocalFile()
        return None

    def dragEnterEvent(self, event):
        if self._first_local_file(event):
            event.acceptProposedAction()
            self.giris.set_drag(True)

    def dragLeaveEvent(self, event):
        self.giris.set_drag(False)

    def dropEvent(self, event):
        self.giris.set_drag(False)
        path = self._first_local_file(event)
        if path:
            self.stack.setCurrentIndex(0)
            self.nav_buttons["giris"].setChecked(True)
            self.giris.load_file(path)
            event.acceptProposedAction()