"""Tema tokenları ve QSS üretimi. Renkleri değiştirmek için sadece bu dosya."""

# Her iki temada da sabit olanlar
COMMON = {
    "sidebar_bg": "#1E152A",
    "sidebar_text": "#D9CFEA",
    "sidebar_hover": "#36224E",
    "primary": "#8A2BE2",
    "primary_text": "#FFFFFF",
    "accent": "#A855F7",      # grafik çizgisi, heliks
    "helix": "#A855F7",
    "sheet": "#EAB308",
    "turn": "#9CA3AF",
    "charge_pos": "#A855F7",
    "charge_neg": "#EC4899",
    "neutral": "#6B7280",
    "error": "#EF4444",
    "success": "#22C55E",
    "warning": "#F59E0B",
}

MODES = {
    "light": {
        "bg": "#F8F6FC",
        "card": "#FFFFFF",
        "text": "#1E152A",
        "muted": "#6B7280",
        "border": "#E4DDF0",
        "tabrow_bg": "#2A1B3D",
        "pill_on": "#E9DDF7",
        "link": "#2563EB",
        "pill_on_text": "#1E152A",
        "drag_bg": "#E9DDF7",
    },
    "dark": {
        "bg": "#161616",
        "card": "#1C1C1C",
        "text": "#F3EEFA",
        "muted": "#9A9A9A",
        "border": "#2E2E2E",
        "tabrow_bg": "#2A1B3D",
        "pill_on": "#EDE7FB",
        "link": "#60A5FA",
        "pill_on_text": "#1E152A",
        "drag_bg": "#2A1B3D",
    },
}


def tokens(mode: str) -> dict:
    return {**COMMON, **MODES[mode]}


_CHECK_PATH = None


def _check_icon_path() -> str:
    """Onay kutusu için ✓ simgesini bir kez çizip geçici klasöre kaydeder."""
    global _CHECK_PATH
    if _CHECK_PATH:
        return _CHECK_PATH
    import os
    import tempfile
    from PyQt5.QtCore import Qt, QPointF
    from PyQt5.QtGui import QColor, QPainter, QPen, QPixmap
    pm = QPixmap(16, 16)
    pm.fill(Qt.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.Antialiasing)
    pen = QPen(QColor("#FFFFFF"), 2)
    pen.setCapStyle(Qt.RoundCap)
    pen.setJoinStyle(Qt.RoundJoin)
    p.setPen(pen)
    p.drawPolyline(QPointF(3.8, 8.4), QPointF(6.8, 11.4), QPointF(12.2, 4.8))
    p.end()
    path = os.path.join(tempfile.gettempdir(), "biomarker_check.png")
    pm.save(path)
    _CHECK_PATH = path.replace("\\", "/")
    return _CHECK_PATH


def build_qss(mode: str) -> str:
    t = tokens(mode)
    check = _check_icon_path()
    from PyQt5.QtGui import QColor
    er, wa = QColor(t['error']), QColor(t['warning'])
    err_bg = f"rgba({er.red()}, {er.green()}, {er.blue()}, 40)"
    warn_bg = f"rgba({wa.red()}, {wa.green()}, {wa.blue()}, 40)"
    return f"""
* {{ font-size: 13px; }}
QMainWindow, QWidget#workspace, QStackedWidget {{
    background: {t['bg']}; color: {t['text']};
}}
QLabel {{ color: {t['text']}; background: transparent; }}
QLabel#pageTitle {{ font-size: 20px; font-weight: 500; }}
QLabel#muted {{ color: {t['muted']}; font-size: 12px; }}

/* Menü çubuğu */
QMenuBar {{ background: {t['sidebar_bg']}; color: {t['sidebar_text']}; }}
QMenuBar::item {{ padding: 5px 10px; background: transparent; }}
QMenuBar::item:selected {{ background: {t['sidebar_hover']}; }}
QMenu {{ background: {t['card']}; color: {t['text']}; border: 1px solid {t['border']}; }}
QMenu::item {{ padding: 6px 24px; }}
QMenu::item:selected {{ background: {t['primary']}; color: {t['primary_text']}; }}

/* Sekme çubuğu */
QFrame#tabRow {{ background: {t['tabrow_bg']}; border: none; }}
QTabBar {{ background: transparent; }}
QTabBar::tab {{
    background: transparent; color: {t['sidebar_text']};
    padding: 8px 14px; margin-right: 2px;
    border-top-left-radius: 6px; border-top-right-radius: 6px;
}}
QTabBar::tab:selected {{ background: {t['bg']}; color: {t['text']}; }}
QTabBar::tab:hover:!selected {{ background: {t['sidebar_hover']}; }}
QToolButton#newTab {{
    background: transparent; color: {t['sidebar_text']};
    border: none; font-size: 18px; padding: 2px 10px;
}}
QToolButton#newTab:hover {{ background: {t['sidebar_hover']}; border-radius: 6px; }}

/* Sidebar */
QFrame#sidebar {{ background: {t['sidebar_bg']}; border: none; }}
QLabel#brand {{ color: {t['accent']}; font-size: 16px; font-weight: 500; padding: 14px 14px 10px 14px; }}
QPushButton#navBtn {{
    background: transparent; color: {t['sidebar_text']};
    text-align: left; padding: 8px 12px; border: none; border-radius: 6px;
}}
QPushButton#navBtn:hover {{ background: {t['sidebar_hover']}; }}
QPushButton#navBtn:checked {{ background: {t['primary']}; color: {t['primary_text']}; }}
QPushButton#navBtn:disabled {{ color: {t['muted']}; background: transparent; }}

/* Kart ve düğmeler */
QFrame#card {{
    background: {t['card']}; border: 1px solid {t['border']}; border-radius: 8px;
}}
QPushButton#primary {{
    background: {t['primary']}; color: {t['primary_text']};
    border: none; border-radius: 6px; padding: 8px 16px; font-weight: 500;
}}
QPushButton#primary:hover {{ background: {t['accent']}; }}
QPushButton#secondary {{
    background: transparent; color: {t['text']};
    border: 1px solid {t['border']}; border-radius: 6px; padding: 7px 14px;
}}
QPushButton#secondary:hover {{ border-color: {t['primary']}; }}

/* Giriş ekranı */
QLabel#cardTitle {{ font-weight: 500; }}
QLabel#error {{ color: {t['error']}; font-size: 12px; }}
QFrame#sep {{ background: {t['border']}; border: none; }}
QPushButton#pill {{
    background: transparent; color: {t['muted']};
    border: none; border-radius: 6px; padding: 6px 14px;
}}
QPushButton#pill:hover {{ color: {t['text']}; }}
QPushButton#pill:checked {{ background: {t['pill_on']}; color: {t['pill_on_text']}; }}
QPushButton#chip {{
    background: transparent; color: {t['text']};
    border: 1px solid {t['border']}; border-radius: 10px; padding: 3px 12px;
}}
QPushButton#chip:hover {{ border-color: {t['primary']}; }}
QPlainTextEdit#seqEdit, QPlainTextEdit#noteEdit, QLineEdit#idEdit {{
    background: {t['bg']}; color: {t['text']};
    border: 1px solid {t['border']}; border-radius: 6px; padding: 8px;
    selection-background-color: {t['primary']}; selection-color: {t['primary_text']};
}}
QPlainTextEdit#seqEdit:focus, QPlainTextEdit#noteEdit:focus, QLineEdit#idEdit:focus {{ border-color: {t['primary']}; }}
QFrame#dropzone {{
    background: transparent; border: 1px dashed {t['primary']}; border-radius: 8px;
}}
QFrame#dropzone[drag="true"] {{ background: {t['drag_bg']}; }}
QFrame#dropzone QLabel {{ color: {t['text']}; }}

/* Özet ekranı */
QLabel#metric {{ font-size: 18px; font-weight: 500; }}
QLabel#badge {{
    background: {t['pill_on']}; color: {t['pill_on_text']};
    border-radius: 6px; padding: 4px 10px;
}}
QPushButton#link {{
    background: transparent; color: {t['link']};
    border: none; padding: 0; font-size: 12px;
}}
QPushButton#link:hover {{ color: {t['text']}; }}

/* Fizikokimyasal ekranı */
QScrollArea {{ background: transparent; border: none; }}
QComboBox, QSpinBox {{
    background: {t['bg']}; color: {t['text']};
    border: 1px solid {t['border']}; border-radius: 6px;
    padding: 4px 8px; min-height: 20px;
}}
QComboBox:hover, QSpinBox:hover {{ border-color: {t['primary']}; }}
QComboBox::drop-down {{ border: none; width: 20px; }}
QComboBox QAbstractItemView {{
    background: {t['card']}; color: {t['text']}; border: 1px solid {t['border']};
    selection-background-color: {t['primary']}; selection-color: {t['primary_text']};
}}
QScrollBar:vertical {{ background: transparent; width: 10px; }}
QScrollBar::handle:vertical {{ background: {t['border']}; border-radius: 5px; min-height: 24px; }}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{ height: 0; }}

/* İkincil yapı ekranı */
QCheckBox {{ color: {t['text']}; spacing: 6px; background: transparent; }}
QCheckBox::indicator {{
    width: 16px; height: 16px; border: 1px solid {t['border']};
    border-radius: 4px; background: {t['bg']};
}}
QCheckBox::indicator:checked {{
    background: {t['primary']}; border-color: {t['primary']}; image: url({check});
}}
QTableWidget {{
    background: transparent; border: none; color: {t['text']}; outline: 0;
    selection-background-color: {t['pill_on']}; selection-color: {t['pill_on_text']};
}}
QTableWidget::item {{ padding: 0 8px; border-bottom: 1px solid {t['border']}; }}
QTableWidget::item:selected,
QTableWidget::item:selected:active,
QTableWidget::item:selected:!active {{
    background: {t['pill_on']}; color: {t['pill_on_text']};
}}
QHeaderView {{ background: transparent; }}
QHeaderView::section {{
    background: transparent; color: {t['muted']}; border: none;
    border-bottom: 1px solid {t['border']}; padding: 0 8px; text-align: left;
}}

/* Homoloji ekranı */
QLineEdit#paramEdit {{
    background: {t['bg']}; color: {t['text']};
    border: 1px solid {t['border']}; border-radius: 6px; padding: 4px 8px; min-height: 20px;
}}
QLineEdit#paramEdit:focus {{ border-color: {t['primary']}; }}
QProgressBar {{ background: {t['border']}; border: none; border-radius: 2px; max-height: 4px; }}
QProgressBar::chunk {{ background: {t['primary']}; border-radius: 2px; }}
QFrame#emptybox {{ background: transparent; border: 1px dashed {t['border']}; border-radius: 8px; }}
QFrame#errorbox {{ background: {err_bg}; border: 1px solid {t['error']}; border-radius: 8px; }}
QFrame#warnbox {{ background: {warn_bg}; border: 1px solid {t['warning']}; border-radius: 8px; }}
QLabel#errorText {{ color: {t['error']}; }}
QLabel#warnText {{ color: {t['warning']}; }}
QLabel#okText {{ color: {t['success']}; font-weight: 500; }}

/* Durum çubuğu */
QStatusBar {{ background: {t['sidebar_bg']}; color: {t['sidebar_text']}; }}
QStatusBar QLabel {{ color: {t['sidebar_text']}; font-size: 12px; }}
QStatusBar::item {{ border: none; }}
"""