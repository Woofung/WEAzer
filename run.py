import sys
from PyQt6.QtWidgets import QApplication
from main_app import EditingAnalyzer

def main():
    app = QApplication(sys.argv)
    analyzer = EditingAnalyzer()
    sys.exit(app.exec())

if __name__ == "__main__":
    main()