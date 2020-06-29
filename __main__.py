from src.QtTest import myWidget
from PyQt5.QtWidgets import QApplication

if __name__ == '__main__':
    app = QApplication([])
    widget = myWidget()
    widget.show()
    app.exec()