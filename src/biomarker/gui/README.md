# BioMarker - Desktop Protein Analysis Interface
# GUI

**State:** In the development phase

BioMarker is a desktop-based bioinformatics analysis tool designed to calculate physicochemical properties, hydropathy profiles, secondary structure predictions, and motif analysis from protein sequences.


## Design & Theme Guidelines
* **Architecture:** PyQt5-based clean, modern, and scientific user interface
* **Theme Support:** System-adaptive **Light** and **Dark** themes


## ⌨️ Keyboard Shortcuts
| Shortcut | Action |
| :--- | :--- |
| **F5** | Start analysis[cite: 7] |
| **Ctrl+N** | New analysis tab[cite: 7] |
| **Ctrl+O** | Import FASTA file[cite: 4, 7] |
| **Ctrl+S** | Save project[cite: 4, 7] |
| **Ctrl+E** | Export PDF report[cite: 4, 7] |

## GUI Architecture
```text
gui/
│
├── BioMarker.sln         # Visual Studio solution file
├── GUI.pyproj            # Project configuration and build file
├── main.py               # Main application entry point and flow manager
├── giris_page.py         # User login and welcome page module
├── window.py             # Core window structure and layout definitions
├── theme.py              # Application color palette, styles, and theming
├── requirements.txt      # GUI-specific dependencies and library list
└── __init__.py           # Python package initializer

