from PyQt6.QtWidgets import (QMainWindow, QTabWidget, QWidget, QVBoxLayout, QLabel)
from sequence_extractor import SequenceExtractor
from editing_analyzer import EditingAnalyzer
from language_settings import LanguageSettings
from project_manager import ProjectManager
from paired_end_merger import PairedEndMerger
from translations import get_translations

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Editing Analyzer")
        self.setGeometry(100, 100, 960, 600)
        self.current_language = "en"
        self.translations = get_translations()
        self.sequence_extractor = SequenceExtractor()
        self.editing_analyzer = EditingAnalyzer()
        self.paired_end_merger = PairedEndMerger()
        self.setup_ui()

    def setup_ui(self):
        central_widget = QWidget()
        layout = QVBoxLayout(central_widget)
        self.project_manager = ProjectManager(self.sequence_extractor, self.editing_analyzer)
        self.tabs = QTabWidget()
        self.language_settings = LanguageSettings()

        self.tabs.addTab(self.project_manager, "Project Management")
        self.tabs.addTab(self.sequence_extractor, "Sequence Extraction")
        self.tabs.addTab(self.editing_analyzer, "Editing Analysis")
        self.tabs.addTab(self.paired_end_merger, "Paired-end Sequence Merge")
        self.tabs.addTab(self.language_settings, "Settings")

        about_tab = QWidget()
        about_layout = QVBoxLayout(about_tab)
        self.about_label = QLabel("")
        self.about_label.setWordWrap(True)
        about_layout.addWidget(self.about_label)
        self.tabs.addTab(about_tab, "About")

        layout.addWidget(self.tabs)
        self.setCentralWidget(central_widget)
        self.setup_connections()
        self.update_ui_text()

    def setup_connections(self):
        self.sequence_extractor.extraction_finished.connect(self.editing_analyzer.set_files)
        self.project_manager.project_loaded.connect(self.editing_analyzer.set_project)
        self.project_manager.project_loaded.connect(self.sequence_extractor.set_project)
        self.paired_end_merger.merge_finished.connect(self.editing_analyzer.set_files)
        self.language_settings.language_changed.connect(self.change_language)

    def change_language(self, language):
        self.current_language = language
        self.update_ui_text()

    def update_ui_text(self):
        translations = self.translations[self.current_language]
        self.setWindowTitle(translations["window_title"])
        self.tabs.setTabText(0, translations["project_management"])
        self.tabs.setTabText(1, translations["sequence_extraction"])
        self.tabs.setTabText(2, translations["editing_analysis"])
        self.tabs.setTabText(3, translations["paired_end_merge"])
        self.tabs.setTabText(4, translations["settings"])
        self.tabs.setTabText(5, translations["about"])
        self.project_manager.set_language(self.current_language)
        self.sequence_extractor.set_language(self.current_language)
        self.editing_analyzer.set_language(self.current_language)
        self.paired_end_merger.set_language(self.current_language)
        self.language_settings.set_language(self.current_language)
        self.about_label.setText(
            translations["about_text"].format(
                version="1.0.5",
                developer="Woofung Tsai",
                copyright="© 2025 Woofung Tsai. All rights reserved.",
                email="woofung@gmail.com"
            )
        )

    def get_project_manager(self):
        return self.project_manager

    def get_sequence_extractor(self):
        return self.sequence_extractor

    def get_editing_analyzer(self):
        return self.editing_analyzer

    def get_paired_end_merger(self):
        return self.paired_end_merger