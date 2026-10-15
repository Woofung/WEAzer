from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox)
from PyQt6.QtCore import pyqtSignal

class LanguageSettings(QWidget):
    language_changed = pyqtSignal(str)

    def __init__(self):
        super().__init__()
        self.current_language = "en"
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()
        lang_layout = QHBoxLayout()
        self.lang_label = QLabel("Interface Language:")
        self.lang_combo = QComboBox()
        self.lang_combo.addItems([
            "Choose language", "English", "Traditional Chinese", "Simplified Chinese",
            "Japanese", "Korean", "French", "German", "Spanish", "Portuguese"
        ])
        self.lang_combo.setCurrentIndex(0)
        self.lang_combo.currentTextChanged.connect(self.change_language)
        lang_layout.addWidget(self.lang_label)
        lang_layout.addWidget(self.lang_combo)
        self.apply_btn = QPushButton("Apply")
        self.apply_btn.clicked.connect(self.apply_settings)
        layout.addLayout(lang_layout)
        layout.addWidget(self.apply_btn)
        layout.addStretch()
        self.setLayout(layout)

    def set_language(self, language):
        self.current_language = language
        translations = {
            "en": {"interface_language": "Interface Language:", "apply": "Apply"},
            "zh_TW": {"interface_language": "介面語言：", "apply": "應用"},
            "zh_CN": {"interface_language": "界面语言：", "apply": "应用"},
            "ja": {"interface_language": "インターフェース言語：", "apply": "適用"},
            "ko": {"interface_language": "인터페이스 언어:", "apply": "적용"},
            "fr": {"interface_language": "Langue de l'interface :", "apply": "Appliquer"},
            "de": {"interface_language": "Schnittstellensprache:", "apply": "Anwenden"},
            "es": {"interface_language": "Idioma de la interfaz:", "apply": "Aplicar"},
            "pt": {"interface_language": "Idioma da interface:", "apply": "Aplicar"}
        }
        t = translations.get(self.current_language, translations["en"])
        self.lang_label.setText(t["interface_language"])
        self.apply_btn.setText(t["apply"])

    def change_language(self, language):
        self.current_language = language

    def apply_settings(self):
        if self.current_language != "Choose language":
            language_map = {
                "English": "en", "Traditional Chinese": "zh_TW", "Simplified Chinese": "zh_CN",
                "Japanese": "ja", "Korean": "ko", "French": "fr", "German": "de",
                "Spanish": "es", "Portuguese": "pt"
            }
            self.language_changed.emit(language_map[self.current_language])