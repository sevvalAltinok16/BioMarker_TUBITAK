import sys
from PyQt5.QtWidgets import QApplication
from window import MainWindow


def main():
    app = QApplication(sys.argv)
    app.setApplicationName("BioMarker")
    win = MainWindow()
    win.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
