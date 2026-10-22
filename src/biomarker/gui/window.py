from PyQt5.QtCore import Qt, QTimer
from PyQt5.QtGui import QKeySequence
from PyQt5.QtWidgets import (
    QMainWindow, QWidget, QHBoxLayout, QVBoxLayout, QTabBar, QToolButton,
    QFrame, QPushButton, QButtonGroup, QStackedWidget, QLabel, QAction,
    QActionGroup, QApplication,
)

import analysis
import theme
from fiziko_page import FizikoPage
from giris_page import GirisPage
from homoloji_page import HomolojiPage
from ikincil_page import IkincilPage
from motif_page import MotifPage
from ozet_page import OzetPage
from raporlar_page import RaporlarPage

# (anahtar, etiket, analiz gerektirir mi)
NAV_ITEMS = [
    ("giris", "Giriş", False),
    ("gecmis", "Geçmiş", False),
    ("fizikokimya", "Fizikokim.", True),
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
        self.tab_infos = {}   # sekme kimliği -> analiz sonucu
        self._next_tab_id = 1

        self._build_menu()
        self._build_body()
        self._build_status()
        self.set_theme("system")
        self.new_tab()

    # ---------- menü ----------
    def _build_menu(self):
        mb = self.menuBar()

        m_file = mb.addMenu("Dosya")
        self._action(m_file, "Yeni analiz", "Ctrl+N", self.new_tab)
        self._action(m_file, "FASTA içe aktar", "Ctrl+O", self.import_fasta)
        self._action(m_file, "Projeyi kaydet", "Ctrl+S")
        self._action(m_file, "Raporu dışa aktar", "Ctrl+E", self.export_report)

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

        # F5: yalnızca Giriş ekranındayken analiz başlatır
        self.act_run = QAction("Analiz et", self)
        self.act_run.setShortcut(QKeySequence("F5"))
        self.act_run.triggered.connect(self.run_shortcut)
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
        self.tabbar.currentChanged.connect(self.refresh_for_tab)
        plus = QToolButton()
        plus.setObjectName("newTab")
        plus.setText("+")
        plus.clicked.connect(self.new_tab)
        rl.addWidget(self.tabbar)
        rl.addWidget(plus)
        rl.addStretch(1)
        col.addWidget(row)

        # sayfalar: önce NAV_ITEMS sırasıyla, sonra Ayarlar, en sonda Özet
        self.stack = QStackedWidget()
        self.page_index = {}
        for key, label, _ in NAV_ITEMS:
            if key == "giris":
                self.giris = GirisPage()
                self.giris.analyze_requested.connect(self.on_analyze)
                page = self.giris
            elif key == "raporlar":
                self.raporlar = RaporlarPage()
                page = self.raporlar
            elif key == "homoloji":
                self.homoloji = HomolojiPage()
                self.homoloji.state_changed.connect(self.on_blast_state)
                page = self.homoloji
            elif key == "motif":
                self.motif = MotifPage()
                self.motif.show_in_3d.connect(self.on_show_in_3d)
                page = self.motif
            elif key == "ikincil":
                self.ikincil = IkincilPage()
                page = self.ikincil
            elif key == "fizikokimya":
                self.fiziko = FizikoPage()
                page = self.fiziko
            else:
                page = PlaceholderPage(label)
            self.page_index[key] = self.stack.addWidget(page)
        self.page_index["ayarlar"] = self.stack.addWidget(PlaceholderPage("Ayarlar"))
        self.ozet = OzetPage()
        self.ozet.navigate.connect(self.go)
        self.page_index["ozet"] = self.stack.addWidget(self.ozet)

        self.stack.currentChanged.connect(self._update_status_right)

        body = QHBoxLayout()
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(0)
        body.addWidget(self._build_sidebar())
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
        self.needs_analysis = []
        for i, (key, label, needs) in enumerate(NAV_ITEMS):
            lay.addWidget(self._nav_button(i, key, label, needs))
            if needs:
                self.needs_analysis.append(key)
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
        if key == "giris":
            b.clicked.connect(lambda _=False: self.show_home())
        else:
            b.clicked.connect(lambda _=False, k=key: self.go(k))
        return b

    # ---------- durum çubuğu ----------
    def _build_status(self):
        sb = self.statusBar()
        sb.setSizeGripEnabled(False)
        self.status_left = QLabel("Durum: Hazır")
        self.status_right = QLabel("Çevrimdışı analiz: açık")
        sb.addWidget(self.status_left, 1)
        sb.addPermanentWidget(self.status_right)

    # ---------- gezinme ----------
    def go(self, key):
        self.stack.setCurrentIndex(self.page_index[key])
        self.nav_buttons[key].setChecked(True)

    def current_result(self):
        i = self.tabbar.currentIndex()
        if i < 0:
            return None
        return self.tab_infos.get(self.tabbar.tabData(i))

    def show_home(self):
        """Analiz yoksa Giriş, varsa Özet."""
        if self.current_result() is not None:
            self.stack.setCurrentIndex(self.page_index["ozet"])
        else:
            self.stack.setCurrentIndex(self.page_index["giris"])
        self.nav_buttons["giris"].setChecked(True)

    def refresh_for_tab(self, *_):
        """Sekme değişince sidebar ve sayfayı o sekmenin durumuna göre ayarla."""
        res = self.current_result()
        has = res is not None
        self.nav_buttons["giris"].setText("Özet" if has else "Giriş")
        for key in self.needs_analysis:
            self.nav_buttons[key].setEnabled(has)
        self.raporlar.set_result(res)
        if has:
            self.ozet.set_result(res)
            self.fiziko.set_result(res)
            self.ikincil.set_result(res)
            self.motif.set_result(res)
            self.homoloji.set_result(res)
        self.show_home()

    # ---------- sekmeler ----------
    def new_tab(self):
        i = self.tabbar.addTab("Yeni analiz")
        self.tabbar.setTabData(i, self._next_tab_id)
        self._next_tab_id += 1
        self.giris.clear()
        self.tabbar.setCurrentIndex(i)
        self.refresh_for_tab()

    def close_tab(self, index):
        if self.tabbar.count() > 1:
            self.tab_infos.pop(self.tabbar.tabData(index), None)
            self.tabbar.removeTab(index)

    # ---------- durum çubuğu metinleri ----------
    def _update_status_right(self, index):
        if index == self.page_index["raporlar"]:
            self.raporlar.refresh()
        if index == self.page_index["homoloji"]:
            self.homoloji.refresh_connection()
            ok = analysis.is_online()
            self.status_right.setText("NCBI BLAST bağlantısı: " + ("açık" if ok else "kapalı"))
        else:
            self.status_right.setText("Çevrimdışı analiz: açık")

    def on_blast_state(self, state):
        self.status_left.setText({"running": "Durum: Analiz ediliyor…", "error": "Durum: Hata"}.get(state, "Durum: Hazır"))

    def export_report(self):
        """Ctrl+E: Raporlar sayfasına geç ve PDF raporu al."""
        self.go("raporlar")
        self.raporlar.export_pdf()

    # ---------- analiz ----------
    def run_shortcut(self):
        if self.stack.currentIndex() == self.page_index["giris"]:
            self.giris.submit()

    def import_fasta(self):
        if self.current_result() is not None:
            self.new_tab()
        self.giris.open_file_dialog()

    def on_analyze(self, kind, value):
        if kind == "kimlik":
            self.giris.show_error("Kimlikle sorgu henüz bağlı değil. Şimdilik diziyi yapıştırın.")
            return
        header, seq = analysis.parse_fasta(value)
        try:
            result = analysis.analyze(seq)
        except analysis.AnalysisError as e:
            self.giris.show_error(str(e))
            return
        accession, entry = analysis.describe(header)
        label = accession if accession == entry else f"{accession} · {entry}"
        info = {"accession": accession, "entry": entry, "label": label, "analysis": result, "seq": seq,
                "ss": analysis.predict_secondary_structure(seq),
                "motifs": analysis.scan_motifs(seq)}

        i = self.tabbar.currentIndex()
        self.tab_infos[self.tabbar.tabData(i)] = info
        self.tabbar.setTabText(i, accession)
        self.tabbar.updateGeometry()
        self.status_left.setText("Durum: Analiz tamamlandı")
        QTimer.singleShot(1500, lambda: self.status_left.setText("Durum: Hazır"))
        self.refresh_for_tab()

    def on_show_in_3d(self, start, end, motif_id):
        """Motif sayfasındaki "3B yapıda göster". 3B sayfası yazılınca aralık oraya iletilecek."""
        self.highlight_3d = (start, end, motif_id)
        self.go("yapi3b")

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
        self.ozet.apply_theme(mode)
        self.fiziko.apply_theme(mode)
        self.ikincil.apply_theme(mode)
        self.motif.apply_theme(mode)

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
            if self.current_result() is not None:
                self.new_tab()  # sonuçlu sekmenin üstüne yazma
            self.go("giris")
            self.giris.load_file(path)
            event.acceptProposedAction()