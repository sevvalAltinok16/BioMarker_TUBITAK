"""Tema tokenları ve QSS üretimi. Renkleri değiştirmek için sadece bu dosya."""

# Her iki temada da sabit olanlar
COMMON = {
    "sidebar_bg": "#1E152A",
    "sidebar_text": "#D9CFEA",
    "sidebar_hover": "#36224E",
    "primary": "#2A1B3D",
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
        "tabrow_bg": "#EDE7F6",
        "pill_on": "#E9DDF7",
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
        "pill_on_text": "#1E152A",
        "drag_bg": "#2A1B3D",
    },
}


def tokens(mode: str) -> dict:
    return {**COMMON, **MODES[mode]}


def build_qss(mode: str) -> str:
    t = tokens(mode)
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
QPlainTextEdit#seqEdit, QLineEdit#idEdit {{
    background: {t['bg']}; color: {t['text']};
    border: 1px solid {t['border']}; border-radius: 6px; padding: 8px;
    selection-background-color: {t['primary']}; selection-color: {t['primary_text']};
}}
QPlainTextEdit#seqEdit:focus, QLineEdit#idEdit:focus {{ border-color: {t['primary']}; }}
QFrame#dropzone {{
    background: transparent; border: 1px dashed {t['primary']}; border-radius: 8px;
}}
QFrame#dropzone[drag="true"] {{ background: {t['drag_bg']}; }}
QFrame#dropzone QLabel {{ color: {t['text']}; }}

/* Durum çubuğu */
QStatusBar {{ background: {t['sidebar_bg']}; color: {t['sidebar_text']}; }}
QStatusBar QLabel {{ color: {t['sidebar_text']}; font-size: 12px; }}
QStatusBar::item {{ border: none; }}
"""