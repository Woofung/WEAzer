from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QTreeWidget, QTreeWidgetItem,
    QLabel, QPushButton, QTableWidget, QTableWidgetItem,
    QMessageBox, QProgressBar, QFileDialog, QRadioButton, QLineEdit,
    QSpinBox, QItemDelegate, QDialog, QComboBox, QScrollArea, QMainWindow,
    QTextEdit, QApplication
)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QPen, QColor, QBrush
from Bio import SeqIO
from Bio import AlignIO
from io import StringIO
import os
import shutil  # Added for potential PATH searching in non-bundled mode
import sys
from collections import Counter, defaultdict
import csv
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.patches import Rectangle
import math
import re
import subprocess
import tempfile
import logging
import numpy as np

# Set up logging
logging.basicConfig(level=logging.WARNING, format='%(asctime)s - %(levelname)s - %(message)s')

# Custom canvas to handle wheel events with tilt support
class ScrollableFigureCanvas(FigureCanvas):
    def __init__(self, figure, scroll_area=None):
        super().__init__(figure)
        self.scroll_area = scroll_area

    def wheelEvent(self, event):
        if self.scroll_area:
            v_scrollbar = self.scroll_area.verticalScrollBar()
            h_scrollbar = self.scroll_area.horizontalScrollBar()
            delta_y = event.angleDelta().y()  # Vertical scroll (up/down)
            delta_x = event.angleDelta().x()  # Horizontal scroll (tilt left/right)
            scroll_factor = 3  # Increase scrolling speed

            # Vertical scrolling
            if delta_y != 0 and v_scrollbar.isVisible():
                step = v_scrollbar.singleStep() * scroll_factor * (-1 if delta_y > 0 else 1)
                v_scrollbar.setValue(v_scrollbar.value() + step)
                event.accept()

            # Horizontal scrolling with tilt
            if delta_x != 0 and h_scrollbar.isVisible():
                step = h_scrollbar.singleStep() * scroll_factor * (-1 if delta_x > 0 else 1)
                h_scrollbar.setValue(h_scrollbar.value() + step)
                event.accept()

        else:
            super().wheelEvent(event)

class VariantSequenceDelegate(QItemDelegate):
    def paint(self, painter, option, index):
        sequence = index.data()
        font = QFont()
        font.setBold(True)
        painter.setFont(font)
        
        x, y = option.rect.x(), option.rect.y()
        for char in sequence:
            if char == 'A':
                painter.setPen(QPen(QColor("#FF5454")))  # Red
            elif char == 'C':
                painter.setPen(QPen(QColor("#26A8FF")))  # Blue
            elif char == 'G':
                painter.setPen(QPen(QColor("#F5F500")))  # Yellow
            elif char == 'T':
                painter.setPen(QPen(QColor("#00D118")))  # Green
            else:
                painter.setPen(QPen(QColor("black")))
            painter.drawText(x, y + option.rect.height() - 5, char)
            x += painter.fontMetrics().horizontalAdvance(char)

class LicenseDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("WEAzer License Agreement")
        self.setMinimumSize(600, 400)
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()
        license_text = """
Woofung Editing Analyzer (WEAzer) License Agreement

License Grant: Woofung Tsai ("Licensor") hereby grants you ("Licensee") a non-exclusive, non-transferable license to use the Woofung Editing Analyzer (WEAzer) software ("Software") for testing and personal use.

Ownership: The Software is owned by Woofung Tsai and is protected by copyright laws and international treaty provisions. This license does not transfer any ownership rights to the Licensee.

Permitted Use: This Software is free to use for personal and research purposes.

Restrictions: The Licensee shall not:
a) Copy, modify, or distribute the Software without explicit written permission from Woofung Tsai.
b) Reverse engineer, decompile, or disassemble the Software.
c) Use the Software for any commercial purposes without prior written agreement from the Licensor.
d) Redistribute the Software in any form without explicit written permission from Woofung Tsai.

No Warranty: The Software is provided "as is" without warranty of any kind, either express or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. Woofung Tsai does not guarantee the accuracy of analysis results produced by the Software.

Limitation of Liability: In no event shall the Licensor be liable for any damages whatsoever arising out of the use of or inability to use this Software. The Licensee is solely responsible for checking the integrity and accuracy of processed data and any results obtained through the use of the Software.

User Responsibility: The Licensee acknowledges and agrees that they are solely responsible for the accuracy and integrity of their own data, reports, documentation, and security when using the Software. The Licensee should independently verify all results obtained through the use of the Software.

Governing Law: This license shall be governed by and construed in accordance with the laws of the United States of America.
"""
        text_edit = QTextEdit()
        text_edit.setReadOnly(True)
        text_edit.setText(license_text)
        layout.addWidget(text_edit)

        button_layout = QHBoxLayout()
        agree_btn = QPushButton(" Agree")
        agree_btn.clicked.connect(self.accept)
        disagree_btn = QPushButton("Disagree")
        disagree_btn.clicked.connect(self.reject)
        button_layout.addStretch()
        button_layout.addWidget(agree_btn)
        button_layout.addWidget(disagree_btn)
        button_layout.addStretch()

        layout.addLayout(button_layout)
        self.setLayout(layout)

class EditingAnalyzer(QWidget):
    def __init__(self):
        super().__init__()
        self.files = []
        self.anchors = []
        self.anchor_type = "left"
        self.range = 0
        self.current_project = ""
        self.project_info = None
        self.barcodes = []
        self.current_language = "en"
        self.file_to_group = {}
        self.file_variants = {}
        self.reference_sequence = None
        self.summary_windows = []
        self.check_license_agreement()
        self.setup_ui()

    def check_license_agreement(self):
        if getattr(sys, 'frozen', False):
            base_path = os.path.join(os.path.dirname(sys.executable), "weazer_data")
            os.makedirs(base_path, exist_ok=True)
            flag_file = os.path.join(base_path, "license_accepted.txt")
        else:
            flag_file = os.path.expanduser("~/.weazer_license_accepted")
        
        if not os.path.exists(flag_file):
            app = QApplication.instance() or QApplication(sys.argv)
            dialog = LicenseDialog(self)
            if dialog.exec() == QDialog.DialogCode.Accepted:
                with open(flag_file, 'w') as f:
                    f.write("License accepted on {}".format(__import__('datetime').datetime.now()))
            else:
                QMessageBox.critical(self, "License Not Accepted", 
                                     "You must agree to the license to use WEAzer. The application will now close.")
                sys.exit(1)

    def setup_ui(self):
        layout = QVBoxLayout()

        file_source_layout = QHBoxLayout()
        file_source_label = QLabel("Step 1. Choose File Source:")
        file_source_label.setStyleSheet("color: #007BFF;")
        self.use_extractor_radio = QRadioButton("Use Files from Sequence Extractor")
        self.use_local_radio = QRadioButton("Choose Local Files")
        self.use_extractor_radio.setChecked(True)
        self.choose_files_btn = QPushButton("Choose Files")
        self.choose_files_btn.setFixedWidth(100)
        self.choose_files_btn.clicked.connect(self.choose_local_files)
        self.choose_files_btn.setEnabled(False)
        file_source_layout.addWidget(file_source_label)
        file_source_layout.addWidget(self.use_extractor_radio)
        file_source_layout.addWidget(self.use_local_radio)
        file_source_layout.addWidget(self.choose_files_btn)
        file_source_layout.addStretch()
        layout.addLayout(file_source_layout)

        self.use_extractor_radio.toggled.connect(self.toggle_file_source)
        self.use_local_radio.toggled.connect(self.toggle_file_source)

        self.project_label = QLabel("Current Project: ")
        layout.addWidget(self.project_label)

        self.folder_path_label = QLabel("Folder Path:")
        layout.addWidget(self.folder_path_label)

        self.tree_widget = QTreeWidget()
        self.tree_widget.setHeaderLabels(["File", "Group", "Replicate"])
        self.tree_widget.setColumnWidth(0, 400)
        self.tree_widget.setColumnWidth(1, 200)
        self.tree_widget.setColumnWidth(2, 100)
        layout.addWidget(self.tree_widget)

        controls_layout = QHBoxLayout()
        self.min_reads_label = QLabel("Minimal Reads Threshold:")
        self.min_reads_combo = QComboBox()
        self.min_reads_combo.addItems(["10", "20", "50", "100"])
        self.min_reads_combo.setCurrentText("10")
        self.min_reads_combo.setEditable(True)
        self.min_reads_combo.setFixedWidth(80)
        controls_layout.addWidget(self.min_reads_label)
        controls_layout.addWidget(self.min_reads_combo)

        self.display_limit_label = QLabel("Variants to Display:")
        self.display_limit_spinbox = QSpinBox()
        self.display_limit_spinbox.setMinimum(1)
        self.display_limit_spinbox.setMaximum(1000)
        self.display_limit_spinbox.setValue(10)
        self.display_limit_spinbox.setFixedWidth(80)
        controls_layout.addWidget(self.display_limit_label)
        controls_layout.addWidget(self.display_limit_spinbox)

        controls_layout.addStretch()
        self.analyze_button = QPushButton("Step 2: Start Analysis")
        self.analyze_button.setFixedWidth(300)
        self.analyze_button.setStyleSheet("QPushButton { background-color: #007BFF; color: white; } QPushButton:hover { color: #007BFF; }")
        self.analyze_button.clicked.connect(self.start_analysis)
        controls_layout.addWidget(self.analyze_button)
        layout.addLayout(controls_layout)

        self.progress_bar = QProgressBar()
        layout.addWidget(self.progress_bar)

        self.result_table = QTableWidget()
        self.result_table.setColumnCount(5)
        self.result_table.setHorizontalHeaderLabels(["Group", "Variant ID", "Variant Sequence", "Frequency", "(%)"])
        self.result_table.setColumnWidth(0, 100)
        self.result_table.setColumnWidth(1, 100)
        self.result_table.setColumnWidth(2, 300)
        self.result_table.setColumnWidth(3, 100)
        self.result_table.setColumnWidth(4, 100)
        self.result_table.setItemDelegateForColumn(2, VariantSequenceDelegate())
        layout.addWidget(self.result_table)

        self.summarize_button = QPushButton("Step 3: Summarize Results")
        self.summarize_button.setStyleSheet("""
            QPushButton { background-color: #28A745; color: white; }
            QPushButton:hover { color: #28A745; }
        """)
        self.summarize_button.clicked.connect(self.summarize)
        layout.addWidget(self.summarize_button)

        self.setLayout(layout)
        self.set_language(self.current_language)

    def set_language(self, language):
        self.current_language = language
        translations = {
            "en": {
                "file_source": "Step 1. Choose File Source:",
                "use_extractor": "Use Files from Sequence Extractor",
                "choose_local": "Choose Local Files",
                "choose_files": "Choose Files",
                "current_project": "Current Project: ",
                "folder_path": "Folder Path:",
                "file": "File",
                "group": "Group",
                "replicate": "Replicate",
                "min_reads": "Minimal Reads Threshold:",
                "display_limit": "Variants to Display:",
                "start_analysis": "Step 2: Start Analysis",
                "summarize": "Step 3: Summarize Results",
                "export_summary": "Export Summary to CSV",
                "plot_variants": "Visualize Result",
                "group_header": "Group",
                "variant_id_header": "Variant ID",
                "variant_sequence_header": "Variant Sequence",
                "frequency_header": "Frequency",
                "percentage_header": "(%)",
                "summary_title": "Summary",
                "sequence_vis": "Sequence visualization:"
            },
            "zh_TW": {
                "file_source": "步驟一：選擇檔案來源：",
                "use_extractor": "使用序列提取器的檔案",
                "choose_local": "選擇本地檔案",
                "choose_files": "選擇檔案",
                "current_project": "當前專案：",
                "folder_path": "檔案路徑：",
                "file": "檔案",
                "group": "組別",
                "replicate": "重複",
                "min_reads": "最小讀數閾值：",
                "display_limit": "顯示的變體數：",
                "start_analysis": "步驟二：開始分析",
                "summarize": "步驟三：總結結果",
                "export_summary": "匯出總結至CSV",
                "plot_variants": "視覺化結果",
                "group_header": "組別",
                "variant_id_header": "變體ID",
                "variant_sequence_header": "變體序列",
                "frequency_header": "頻率",
                "percentage_header": "(%)",
                "summary_title": "總結",
                "sequence_vis": "序列可視化："
            },
            "zh_CN": {
                "file_source": "步骤一：选择文件来源：",
                "use_extractor": "使用序列提取器的文件",
                "choose_local": "选择本地文件",
                "choose_files": "选择文件",
                "current_project": "当前项目：",
                "folder_path": "文件路径：",
                "file": "文件",
                "group": "组别",
                "replicate": "重复",
                "min_reads": "最小读数阈值：",
                "display_limit": "显示的变体数：",
                "start_analysis": "步骤二：开始分析",
                "summarize": "步骤三：总结结果",
                "export_summary": "导出总结至CSV",
                "plot_variants": "可视化结果",
                "group_header": "组别",
                "variant_id_header": "变体ID",
                "variant_sequence_header": "变体序列",
                "frequency_header": "频率",
                "percentage_header": "(%)",
                "summary_title": "总结",
                "sequence_vis": "序列可视化："
            },
            "ko": {
                "file_source": "1단계: 파일 소스 선택:",
                "use_extractor": "시퀀스 추출기에서 파일 사용",
                "choose_local": "로컬 파일 선택",
                "choose_files": "파일 선택",
                "current_project": "현재 프로젝트: ",
                "folder_path": "폴더 경로:",
                "file": "파일",
                "group": "그룹",
                "replicate": "반복",
                "min_reads": "최소 리드 임계값:",
                "display_limit": "표시할 변이 수:",
                "start_analysis": "2단계: 분석 시작",
                "summarize": "3단계: 결과 요약",
                "export_summary": "CSV로 요약 내보내기",
                "plot_variants": "결과 시각화",
                "group_header": "그룹",
                "variant_id_header": "변이 ID",
                "variant_sequence_header": "변이 시퀀스",
                "frequency_header": "빈도",
                "percentage_header": "(%)",
                "summary_title": "요약",
                "sequence_vis": "시퀀스 시각화:"
            }
        }
        t = translations.get(self.current_language, translations["en"])
        file_source_label = self.layout().itemAt(0).layout().itemAt(0).widget()
        file_source_label.setText(t["file_source"])
        self.use_extractor_radio.setText(t["use_extractor"])
        self.use_local_radio.setText(t["choose_local"])
        self.choose_files_btn.setText(t["choose_files"])
        self.project_label.setText(f"{t['current_project']}{self.current_project}")
        self.folder_path_label.setText(t["folder_path"])
        self.tree_widget.setHeaderLabels([t["file"], t["group"], t["replicate"]])
        self.min_reads_label.setText(t["min_reads"])
        self.display_limit_label.setText(t["display_limit"])
        self.analyze_button.setText(t["start_analysis"])
        self.summarize_button.setText(t["summarize"])
        self.result_table.setHorizontalHeaderLabels([
            t["group_header"],
            t["variant_id_header"],
            t["variant_sequence_header"],
            t["frequency_header"],
            t["percentage_header"]
        ])

    def toggle_file_source(self):
        self.choose_files_btn.setEnabled(self.use_local_radio.isChecked())
        self.tree_widget.clear()
        self.files = []
        if self.use_extractor_radio.isChecked() and self.project_info:
            self.set_files(self.files)

    def choose_local_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select FASTQ Files", "", "FASTQ Files (*.fastq)")
        if files:
            self.files = files
            self.set_files(files)

    def set_project(self, project_info):
        self.project_info = project_info
        self.current_project = project_info.get('name', '')
        self.anchors = project_info.get('anchors', [])
        self.anchor_type = project_info.get('anchor_type', 'left')
        self.range = int(project_info.get('range', '0')) if project_info.get('range', '') else 0
        self.barcodes = [seq.get('barcode', '') for seq in project_info.get('sequences', [])]
        self.reference_sequence = project_info.get('reference_sequence', '')
        self.project_label.setText(f"Current Project: {self.current_project}")
        if self.use_extractor_radio.isChecked():
            self.set_files(self.files)

    def parse_file_name(self, file_name):
        match = re.match(r"(.*)-\[Rep(\d+)\]-(\w+)\.fastq", file_name)
        if match:
            base_name, rep_number, barcode = match.groups()
            group = base_name
            replicate = rep_number
        else:
            base_name = os.path.splitext(file_name)[0]
            group = base_name
            replicate = "1"
        return group, replicate

    def set_files(self, files):
        self.files = files
        self.tree_widget.clear()
        if files:
            self.folder_path_label.setText(f"Folder Path: {os.path.dirname(files[0])}")
        for file in files:
            file_name = os.path.basename(file)
            group, replicate = self.parse_file_name(file_name)
            parent_item = QTreeWidgetItem(self.tree_widget, [file_name])
            child_item = QTreeWidgetItem(parent_item, [file_name, group, f"Replicate {replicate}"])
            group_edit = QLineEdit(group)
            self.tree_widget.setItemWidget(child_item, 1, group_edit)
        self.tree_widget.expandAll()

    def start_analysis(self):
        if not self.files:
            QMessageBox.warning(self, "Warning", "Please set the files to analyze first")
            return
        if not self.anchors:
            QMessageBox.warning(self, "Warning", "Please set the anchor sequences first")
            return
        if self.anchor_type != "flanking" and self.range <= 0:
            QMessageBox.warning(self, "Warning", "Please set a valid analyzing range")
            return

        self.file_to_group = {}
        for top_level_item in range(self.tree_widget.topLevelItemCount()):
            parent_item = self.tree_widget.topLevelItem(top_level_item)
            file_name = parent_item.text(0)
            child_item = parent_item.child(0)
            group_name = self.tree_widget.itemWidget(child_item, 1).text()
            self.file_to_group[file_name] = group_name

        self.result_table.setRowCount(0)
        self.progress_bar.setMaximum(len(self.files))
        self.progress_bar.setValue(0)

        min_reads = int(self.min_reads_combo.currentText()) if self.min_reads_combo.currentText().isdigit() else 10
        display_limit = self.display_limit_spinbox.value()

        self.file_variants = {}
        for file_idx, file in enumerate(self.files):
            variants = self.analyze_file(file)
            total_reads = sum(variants.values())
            sorted_variants = sorted(variants.items(), key=lambda x: x[1], reverse=True)
            filtered_variants = [(seq, count) for seq, count in sorted_variants if count >= min_reads][:display_limit]

            group_name = self.file_to_group.get(os.path.basename(file), "Unknown")
            for idx, (variant_seq, freq) in enumerate(filtered_variants, start=1):
                row = self.result_table.rowCount()
                self.result_table.insertRow(row)
                self.result_table.setItem(row, 0, QTableWidgetItem(group_name))
                self.result_table.setItem(row, 1, QTableWidgetItem(f"{idx:02d}"))
                self.result_table.setItem(row, 2, QTableWidgetItem(variant_seq))
                self.result_table.setItem(row, 3, QTableWidgetItem(str(freq)))
                percentage = (freq / total_reads * 100) if total_reads > 0 else 0
                self.result_table.setItem(row, 4, QTableWidgetItem(f"{percentage:.2f}%"))

            self.progress_bar.setValue(file_idx + 1)
            self.file_variants[file] = {variant_seq: freq for variant_seq, freq in filtered_variants}

        file_colors = {}
        for idx, file in enumerate(self.files):
            color = QBrush(QColor(173, 216, 230)) if idx % 2 == 0 else QBrush(QColor(255, 255, 255))
            file_colors[file] = color

        for row in range(self.result_table.rowCount()):
            group_name = self.result_table.item(row, 0).text()
            for file in self.files:
                if self.file_to_group.get(os.path.basename(file), "Unknown") == group_name:
                    color = file_colors.get(file, QBrush(QColor(255, 255, 255)))
                    break
            for col in range(self.result_table.columnCount()):
                item = self.result_table.item(row, col)
                item.setBackground(color)

    def analyze_file(self, file):
        variants = Counter()
        try:
            with open(file, 'r') as handle:
                for record in SeqIO.parse(handle, "fastq"):
                    try:
                        seq = str(record.seq)
                        if self.anchor_type == "left":
                            anchor = self.anchors[0]
                            if anchor in seq:
                                start = seq.index(anchor) + len(anchor)
                                variant = seq[start:start + self.range]
                                if len(variant) == self.range:
                                    variants[variant] += 1
                        elif self.anchor_type == "right":
                            anchor = self.anchors[0]
                            if anchor in seq:
                                end = seq.index(anchor)
                                variant = seq[max(0, end - self.range):end]
                                if len(variant) == self.range:
                                    variants[variant] += 1
                        elif self.anchor_type == "flanking":
                            left_anchor, right_anchor = self.anchors
                            if left_anchor in seq and right_anchor in seq:
                                start = seq.index(left_anchor) + len(left_anchor)
                                end = seq.index(right_anchor)
                                variant = seq[start:end]
                                variants[variant] += 1
                    except ValueError as e:
                        logging.warning(f"Skipping invalid record in {file}: {record.id} - {str(e)}")
                        continue
        except Exception as e:
            QMessageBox.warning(self, "File Error", f"Failed to process {file}: {str(e)}")
        return variants

    def summarize(self):
        if not self.file_variants or not all(file in self.file_variants for file in self.files):
            QMessageBox.warning(self, "Warning", "Please run 'Start Analysis' before summarizing.")
            return

        all_variants = set().union(*[set(variants.keys()) for variants in self.file_variants.values()])
        total_reads = {variant: sum(self.file_variants[file].get(variant, 0) for file in self.files) 
                       for variant in all_variants}
        sorted_variants = sorted(all_variants, key=lambda v: total_reads[v], reverse=True)
        total_file_reads = {file: sum(self.file_variants[file].values()) for file in self.files}

        group_to_files = defaultdict(list)
        for file in self.files:
            file_name = os.path.basename(file)
            group = self.file_to_group.get(file_name, "Unknown")
            group_to_files[group].append(file)

        sorted_groups = sorted(group_to_files.keys())
        column_types = []
        for group in sorted_groups:
            files = group_to_files[group]
            for idx, file in enumerate(files):
                column_types.append(("file", file, group, idx+1))
            if len(files) > 1:
                column_types.append(("average", group))

        headers = ["Variant Sequence"]
        for col_type in column_types:
            if col_type[0] == "file":
                headers.append(f"{col_type[2]} Rep{col_type[3]}")
            elif col_type[0] == "average":
                headers.append(f"{col_type[1]} Average")

        summary_window = QDialog(self)
        t = self.get_translations()
        summary_window.setWindowTitle(t["summary_title"])
        summary_layout = QVBoxLayout()

        summary_table = QTableWidget()
        summary_table.setColumnCount(len(headers))
        summary_table.setHorizontalHeaderLabels(headers)

        for row, variant in enumerate(sorted_variants):
            summary_table.insertRow(row)
            summary_table.setItem(row, 0, QTableWidgetItem(variant))
            col_idx = 1
            for col_type in column_types:
                if col_type[0] == "file":
                    file = col_type[1]
                    read_count = self.file_variants[file].get(variant, 0)
                    total = total_file_reads[file]
                    percentage = (read_count / total * 100) if total > 0 else 0
                    item = QTableWidgetItem(f"{percentage:.2f}%")
                    summary_table.setItem(row, col_idx, item)
                elif col_type[0] == "average":
                    group = col_type[1]
                    files = group_to_files[group]
                    percentages = []
                    for file in files:
                        read_count = self.file_variants[file].get(variant, 0)
                        total = total_file_reads[file]
                        percentage = (read_count / total * 100) if total > 0 else 0
                        percentages.append(percentage)
                    avg_percentage = sum(percentages) / len(percentages) if percentages else 0
                    item = QTableWidgetItem(f"{avg_percentage:.2f}%")
                    summary_table.setItem(row, col_idx, item)
                col_idx += 1

        light_blue = QBrush(QColor("#ADD8E6"))
        for col_idx, col_type in enumerate(column_types, start=1):
            if col_type[0] == "average":
                for row in range(summary_table.rowCount()):
                    item = summary_table.item(row, col_idx)
                    if item:
                        item.setBackground(light_blue)

        summary_layout.addWidget(summary_table)

        export_button = QPushButton(t["export_summary"])
        export_button.clicked.connect(lambda: self.export_summary_to_csv(summary_table))
        summary_layout.addWidget(export_button)

        group_layout = QHBoxLayout()
        group_label = QLabel(t["sequence_vis"])
        group_combo = QComboBox()
        group_combo.addItems(sorted(set(self.file_to_group.values())))
        group_combo.insertItem(0, "All")
        group_combo.setCurrentIndex(0)
        group_layout.addWidget(group_label)
        group_layout.addWidget(group_combo)
        group_layout.addStretch()
        summary_layout.addLayout(group_layout)

        status_layout = QHBoxLayout()
        self.status_label = QLabel("Plotting status: Not started")
        status_layout.addStretch()
        status_layout.addWidget(self.status_label)
        status_layout.addStretch()
        summary_layout.addLayout(status_layout)

        plot_button = QPushButton(t["plot_variants"])
        plot_button.setStyleSheet("""
            QPushButton { background-color: #28A745; color: white; }
            QPushButton:hover { color: #28A745; }
        """)
        plot_button.clicked.connect(lambda: self.plot_variant_alignment(group_combo.currentText()))
        summary_layout.addWidget(plot_button)

        summary_window.setLayout(summary_layout)
        summary_window.resize(1200, 900)
        
        summary_window.show()
        self.summary_windows.append(summary_window)

    def get_translations(self):
        translations = {
            "en": {
                "export_summary": "Export Summary to CSV",
                "plot_variants": "Visualize Result",
                "summary_title": "Summary",
                "sequence_vis": "Sequence visualization:"
            },
            "zh_TW": {
                "export_summary": "匯出總結至CSV",
                "plot_variants": "視覺化結果",
                "summary_title": "總結",
                "sequence_vis": "序列可視化："
            },
            "zh_CN": {
                "export_summary": "导出总结至CSV",
                "plot_variants": "可视化结果",
                "summary_title": "总结",
                "sequence_vis": "序列可视化："
            },
            "ko": {
                "export_summary": "CSV로 요약 내보내기",
                "plot_variants": "결과 시각화",
                "summary_title": "요약",
                "sequence_vis": "시퀀스 시각화:"
            }
        }
        return translations.get(self.current_language, translations["en"])

    def export_summary_to_csv(self, table):
        file_path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV Files (*.csv)")
        if file_path:
            with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                writer.writerow([table.horizontalHeaderItem(i).text() for i in range(table.columnCount())])
                for row in range(table.rowCount()):
                    row_data = [table.item(row, col).text() for col in range(table.columnCount())]
                    writer.writerow(row_data)
            QMessageBox.information(self, "Success", "CSV export successful!")

    def plot_variant_alignment(self, selected_group):
        if not self.file_variants:
            QMessageBox.warning(self, "Warning", "No variant data to plot.")
            return

        self.status_label.setText("Plotting status: Initiating (0%) Please be patient.")
        QApplication.processEvents()

        # Determine effective range for plotting decision
        if self.anchor_type == "flanking":
            all_variants = set().union(*[set(variants.keys()) for variants in self.file_variants.values()])
            if self.reference_sequence:
                effective_range = len(self.reference_sequence)
            else:
                variant_lengths = [len(v) for v in all_variants]
                effective_range = max(variant_lengths, default=0) if variant_lengths else 0
        else:
            effective_range = self.range

        # Decide plot type based on effective range
        if effective_range <= 3:
            self.plot_stacked_bar(selected_group)
        else:
            self.plot_alignment(selected_group)

    def plot_stacked_bar(self, selected_group):
        # Prepare data
        if selected_group == "All":
            files_to_plot = self.files
        else:
            files_to_plot = [f for f in self.files if self.file_to_group.get(os.path.basename(f), "Unknown") == selected_group]
            if not files_to_plot:
                QMessageBox.warning(self, "Warning", f"No data for group {selected_group}")
                return

        groups = sorted(set(self.file_to_group.get(os.path.basename(f), "Unknown") for f in files_to_plot))  # Sort groups alphabetically
        all_variants = sorted(set().union(*[set(self.file_variants[f].keys()) for f in files_to_plot]), 
                             key=lambda x: sum(self.file_variants[f].get(x, 0) for f in files_to_plot), reverse=True)
        display_limit = self.display_limit_spinbox.value()
        top_variants = all_variants[:display_limit]
        
        # Aggregate percentages by group
        group_percentages = {group: {variant: 0 for variant in top_variants} for group in groups}
        for group in groups:
            group_files = [f for f in files_to_plot if self.file_to_group.get(os.path.basename(f), "Unknown") == group]
            total_reads = sum(sum(self.file_variants[f].values()) for f in group_files)
            if total_reads == 0:
                continue
            variant_counts = {v: sum(self.file_variants[f].get(v, 0) for f in group_files) for v in top_variants}
            for variant in top_variants:
                group_percentages[group][variant] = (variant_counts[variant] / total_reads * 100) if total_reads > 0 else 0

        # Plot setup
        dpi = 100
        fig_width = max(8, len(groups) * 0.8)  # Adjust width based on number of groups
        fig_height = 6
        fig = plt.Figure(figsize=(fig_width, fig_height), dpi=dpi)
        ax = fig.add_subplot(111)

        # Define colors
        single_nucleotide_colors = {
            'A': '#FF5454',  # Red (matches alignment plot)
            'G': '#F5F500',  # Yellow (matches alignment plot)
            'C': '#26A8FF',  # Blue (matches alignment plot)
            'T': '#00D118'   # Green (matches alignment plot)
        }
        # Colorblind-friendly palette (ColorBrewer Set1, adjusted for clarity)
        colorblind_colors = [
            '#E41A1C',  # Red
            '#377EB8',  # Blue
            '#4DAF4A',  # Green
            '#984EA3',  # Purple
            '#FF7F00',  # Orange
            '#FFFF33',  # Yellow
            '#A65628',  # Brown
            '#F781BF',  # Pink
            '#999999'   # Grey
        ]

        # Check if all variants are single nucleotides
        all_single_nucleotides = all(len(v) == 1 and v in 'AGCT' for v in top_variants)
        
        if all_single_nucleotides:
            variant_colors = [single_nucleotide_colors[v] for v in top_variants]
        else:
            # Use colorblind palette, cycling if needed
            variant_colors = [colorblind_colors[i % len(colorblind_colors)] for i in range(len(top_variants))]

        # Plot stacked bars
        x = np.arange(len(groups))
        bottom = np.zeros(len(groups))
        for variant, color in zip(top_variants, variant_colors):
            heights = [group_percentages[group][variant] for group in groups]
            ax.bar(x, heights, bottom=bottom, label=variant, color=color, edgecolor='white')
            bottom += heights

        # Customize plot
        ax.set_xticks(x)
        ax.set_xticklabels(groups, rotation=45, ha='right')
        ax.set_ylabel('Percentage (%)')
        ax.set_ylim(0, 100)
        ax.set_title(f'Variant Distribution ({selected_group})')
        ax.legend(title='Variants', bbox_to_anchor=(1.05, 1), loc='upper left')

        # Finalize
        self.finalize_plot(fig, "Variant Stacked Bar Plot")

    def plot_alignment(self, selected_group):
        colors = {
            'A': '#FF5454', 'C': '#26A8FF', 'G': '#F5F500', 'T': '#00D118',
            '-': '#D3D3D3', 'N': '#B3B3B3', ' ': '#FFFFFF'  # White for padding
        }
        insertion_color = '#43464B'
        default_text_color = '#404340'
        insertion_text_color = '#FFFFFF'

        # Pixel-based constants
        dpi = 100
        box_pixels = 20      # 20x20 pixel boxes
        spacer_pixels = 2    # 2 pixel spacers
        annotation_pixels = 80  # 80 pixels for annotations on both sides
        box_size_inch = box_pixels / dpi  # 0.2 inches
        spacer_inch = spacer_pixels / dpi  # 0.02 inches
        annotation_space_inch = annotation_pixels / dpi  # 0.8 inches

        # Data coordinates
        box_size = 0.48
        spacer = 0.1
        v_spacing = 0.15
        char_font_size = 12.5
        tick_font_size = 8
        reads_font_size = 10
        font_family = 'Courier New'

        if selected_group == "All":
            files_to_plot = self.files
        else:
            files_to_plot = [f for f in self.files if self.file_to_group.get(os.path.basename(f), "Unknown") == selected_group]
            if not files_to_plot:
                QMessageBox.warning(self, "Warning", f"No data for group {selected_group}")
                return

        num_files = len(files_to_plot)
        max_seq_length = max([len(str(seq)) for file in files_to_plot for seq in self.file_variants.get(file, {}).keys()] + 
                             ([len(self.reference_sequence)] if self.reference_sequence else []))
        if max_seq_length > 30:
            cols = 1
        elif max_seq_length < 15:
            cols = 4
        else:
            cols = 2
        rows = math.ceil(num_files / cols)

        # Calculate dimensions
        subplot_width_pixels = (max_seq_length * box_pixels + 
                                (max_seq_length - 1) * spacer_pixels + 
                                2 * annotation_pixels)
        subplot_width_inch = subplot_width_pixels / dpi
        fig_width = subplot_width_inch * cols

        # Figure setup
        fig = plt.Figure(dpi=dpi)
        if num_files > 1:
            axes = [fig.add_subplot(rows, cols, i + 1) for i in range(num_files)]
            fig_height = box_size_inch * (11 * (1 + v_spacing / box_size) + 2) * rows
            fig.set_size_inches(fig_width, fig_height)
        else:
            axes = [fig.add_subplot(111)]

        # Plot each subplot
        for ax_idx, (ax, file, data) in enumerate(zip(axes, files_to_plot, self.process_alignment_data(files_to_plot))):
            if data is None:
                ax.axis('off')
                ax.text(0.5, 0.5, "No Variants" if self.file_variants.get(file, {}) else "No Reads", 
                        ha='center', va='center', fontsize=10, fontfamily=font_family)
                progress = (ax_idx + 1) / num_files * 100
                status = "Initiating" if progress <= 25 else "Completing" if progress < 100 else "Completed"
                self.status_label.setText(f"Plotting status: {status} ({progress:.0f}%)")
                QApplication.processEvents()
                continue

            padded_sequences, variant_list, total_reads, orig_len = data
            num_rows = len(padded_sequences)  # Reference + variants

            # Set limits based on content
            x_max = max_seq_length * (box_size + spacer) + 2 * annotation_space_inch * dpi / box_pixels
            y_max = num_rows * (box_size + v_spacing) + 1.5  # Adjust for title and ticks
            ax.set_xlim(0, x_max)
            ax.set_ylim(0, y_max)
            ax.set_aspect('equal', adjustable='box')  # Enforce square boxes
            ax.axis('off')

            # Adjust figure height for single subplot
            if num_files == 1:
                fig_height = box_size_inch * (num_rows * (1 + v_spacing / box_size) + 2)
                fig.set_size_inches(subplot_width_inch, fig_height)

            x_offset = 0.5
            y_offset = (num_rows - 1) * (box_size + v_spacing) + 0.5
            group_name = self.file_to_group.get(os.path.basename(file), "Unknown")
            file_name = os.path.basename(file)
            title = f"{group_name}\n({file_name})" if sum(1 for f in self.files if self.file_to_group.get(os.path.basename(f), "Unknown") == group_name) > 1 else group_name
            ax.text(x_offset, y_offset + 1.0, title, fontsize=12, fontfamily=font_family, ha='left', color=default_text_color)

            # Numbering only up to original length
            tick_locations = [1, orig_len] + list(range(1, orig_len + 1, 10))[1:]
            tick_locations = sorted(set(tick_locations))
            for x in tick_locations:
                ax.text(x_offset + (x - 1) * (box_size + spacer) + box_size / 2, y_offset + 0.8, str(x),
                        ha='center', fontsize=tick_font_size, fontfamily=font_family, color=default_text_color)

            y = y_offset
            reference_aligned = padded_sequences[0]
            for i, c in enumerate(reference_aligned):
                x = x_offset + i * (box_size + spacer)
                fill_color = colors.get(c, '#B3B3B3')
                ax.add_patch(Rectangle((x, y), box_size, box_size, fill=True, color=fill_color))
                if c != ' ':  # Only add text for non-padding
                    ax.text(x + box_size / 2, y + box_size / 2 - 0.05, c, ha='center', va='center',
                            fontsize=char_font_size, fontfamily=font_family, color=default_text_color)

            reads_x = x_offset + max_seq_length * (box_size + spacer) + 0.5
            ax.text(reads_x, y + box_size / 2, "Reads", fontsize=reads_font_size, fontfamily=font_family, va='center', color=default_text_color)
            ax.axvline(x=reads_x + 1.5 * 1.25, ymin=0, ymax=1, color='grey', lw=0.5)
            ax.text(reads_x + 2.3 * 1.25, y + box_size / 2, "(%)", fontsize=reads_font_size, fontfamily=font_family, va='center', color=default_text_color)

            for line_number, (aligned_seq, (_, reads)) in enumerate(zip(padded_sequences[1:], variant_list), start=1):
                y = y_offset - line_number * (box_size + v_spacing)
                for i, c in enumerate(aligned_seq):
                    x = x_offset + i * (box_size + spacer)
                    ref_char = reference_aligned[i]
                    if c == ref_char and c != '-' and c != ' ':
                        ax.text(x + box_size / 2, y + box_size / 2 - 0.05, '•', ha='center', va='center',
                                fontsize=char_font_size, fontfamily=font_family, color=default_text_color)
                    else:
                        fill_color = colors.get(c, '#B3B3B3')
                        char_color = default_text_color
                        if ref_char == '-' and c != '-' and c != ' ':
                            fill_color = insertion_color
                            char_color = insertion_text_color
                        ax.add_patch(Rectangle((x, y), box_size, box_size, fill=True, color=fill_color))
                        if c != ' ':  # Only add text for non-padding
                            ax.text(x + box_size / 2, y + box_size / 2 - 0.05, c, ha='center', va='center',
                                    fontsize=char_font_size, fontfamily=font_family, color=char_color)

                ax.text(reads_x, y + box_size / 2, f"{reads}", fontsize=reads_font_size, fontfamily=font_family, va='center', color=default_text_color)
                ax.text(reads_x + 1.7 * 1.25, y + box_size / 2, f"{(reads / total_reads * 100):.2f}%",
                        fontsize=reads_font_size, fontfamily=font_family, va='center', color=default_text_color)

            progress = (ax_idx + 1) / num_files * 100
            status = "Initiating" if progress <= 25 else "Completing" if progress < 100 else "Completed"
            self.status_label.setText(f"Plotting status: {status} ({progress:.0f}%)")
            QApplication.processEvents()

        # Finalize
        self.finalize_plot(fig, "Variant Alignment Plot")

    def process_alignment_data(self, files_to_plot):
        aligned_data = []
        original_lengths = []

        for file in files_to_plot:
            variants = self.file_variants.get(file, {})
            if not variants:
                aligned_data.append(None)
                original_lengths.append(0)
                continue

            total_reads = sum(variants.values())
            if total_reads == 0:
                aligned_data.append(None)
                original_lengths.append(0)
                continue

            variant_list = sorted(variants.items(), key=lambda x: x[1], reverse=True)
            reference = self.reference_sequence if self.reference_sequence else variant_list[0][0]
            sequences = [reference] + [seq for seq, _ in variant_list]
            fasta_input = "\n".join([f">seq{i}\n{seq}" for i, seq in enumerate(sequences)])

            with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_input:
                temp_input.write(fasta_input)
                temp_input_path = temp_input.name

            try:
                # Determine MAFFT path
                if getattr(sys, 'frozen', False):
                    # Running as PyInstaller executable, use bundled MAFFT
                    mafft_path = os.path.join(sys._MEIPASS, 'mafft-win', 'mafft.bat')
                    if not os.path.exists(mafft_path):
                        raise FileNotFoundError("Bundled MAFFT not found in the application package. Please ensure MAFFT is correctly bundled.")
                else:
                    # Running from source, try default path or system PATH
                    mafft_path = r"C:\Program Files\MAFFT\mafft-win\mafft.bat"
                    if not os.path.exists(mafft_path):
                        mafft_path = shutil.which("mafft")
                        if not mafft_path:
                            raise FileNotFoundError("MAFFT not found in default path or system PATH. Please install MAFFT manually or bundle it with the application.")

                # Run MAFFT
                mafft_cmd = [mafft_path, "--auto", temp_input_path]
                process = subprocess.Popen(mafft_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, shell=True)
                stdout, stderr = process.communicate()
                os.unlink(temp_input_path)

                if process.returncode != 0 or not stdout.strip():
                    error_msg = f"MAFFT failed for {os.path.basename(file)}:\nSTDERR: {stderr}"
                    raise Exception(error_msg)

                alignment = AlignIO.read(StringIO(stdout), "fasta")
                aligned_sequences = [str(record.seq).upper() for record in alignment]
                original_lengths.append(len(aligned_sequences[0]))
                aligned_data.append((aligned_sequences, variant_list, total_reads))

            except FileNotFoundError as e:
                if os.path.exists(temp_input_path):
                    os.unlink(temp_input_path)
                QMessageBox.critical(self, "MAFFT Not Found", str(e))
                aligned_data.append(None)
                original_lengths.append(0)
                continue
            except Exception as e:
                if os.path.exists(temp_input_path):
                    os.unlink(temp_input_path)
                QMessageBox.critical(self, "Alignment Error", f"Failed to align sequences: {str(e)}")
                aligned_data.append(None)
                original_lengths.append(0)
                continue

        max_aligned_length = max(original_lengths + [0])
        padded_aligned_data = []
        for data, orig_len in zip(aligned_data, original_lengths):
            if data is None:
                padded_aligned_data.append(None)
                continue
            aligned_sequences, variant_list, total_reads = data
            padded_sequences = [seq + ' ' * (max_aligned_length - len(seq)) for seq in aligned_sequences]
            padded_aligned_data.append((padded_sequences, variant_list, total_reads, orig_len))
        
        return padded_aligned_data

    def finalize_plot(self, fig, title):
        # Set up Qt canvas with standard toolbar and scrolling
        scroll_area = QScrollArea()
        canvas = ScrollableFigureCanvas(fig, scroll_area=scroll_area)
        canvas.setMinimumSize(int(fig.get_figwidth() * fig.dpi), int(fig.get_figheight() * fig.dpi))
        canvas.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        toolbar = NavigationToolbar(canvas, self)
        layout = QVBoxLayout()
        layout.addWidget(toolbar)
        
        scroll_area.setWidget(canvas)
        scroll_area.setWidgetResizable(True)
        scroll_area.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll_area.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        
        # Apply thinner blue scrollbar stylesheet
        scroll_area.setStyleSheet("""
            QScrollBar:vertical {
                background: #ADD8E6;
                width: 8px;
                margin: 0px;
            }
            QScrollBar::handle:vertical {
                background: #007BFF;
                min-height: 20px;
            }
            QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
                height: 0px;
                background: none;
            }
            QScrollBar:horizontal {
                background: #ADD8E6;
                height: 8px;
                margin: 0px;
            }
            QScrollBar::handle:horizontal {
                background: #007BFF;
                min-width: 20px;
            }
            QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {
                width: 0px;
                background: none;
            }
        """)
        layout.addWidget(scroll_area)

        plot_window = QMainWindow(self)
        plot_window.setWindowTitle(title)
        central_widget = QWidget()
        plot_window.setCentralWidget(central_widget)
        central_widget.setLayout(layout)

        fig.tight_layout()
        plot_window.showMaximized()

        self.status_label.setText("Plotting status: Completed (100%)")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = EditingAnalyzer()
    window.show()
    sys.exit(app.exec())