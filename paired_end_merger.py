from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox, QTableWidget, QTableWidgetItem,
    QFileDialog, QProgressBar, QMessageBox, QApplication, QTextEdit
)
from PyQt6.QtCore import pyqtSignal, QThread, QObject
from pathlib import Path
import os
import re
from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import time

class MergeWorker(QObject):
    progress = pyqtSignal(int)
    finished = pyqtSignal(list)
    status = pyqtSignal(str)
    report = pyqtSignal(str)
    update_status = pyqtSignal(str)

    def __init__(self, grouped_files, min_overlap):
        super().__init__()
        self.grouped_files = grouped_files
        self.min_overlap = min_overlap
        self.aborted = False

    def run(self):
        merged_files = []
        total_groups = len(self.grouped_files)
        total_reads_processed = 0

        # Calculate total reads across all files for progress estimation
        total_reads = 0
        for group, files in self.grouped_files.items():
            r1_file = files["R1"]
            total_reads += sum(1 for _ in SeqIO.parse(gzip.open(r1_file, "rt") if r1_file.endswith(".gz") else r1_file, "fastq"))
        
        self.progress.emit(0)

        # Create output directory
        first_file_path = next(iter(self.grouped_files.values()))["R1"]
        output_dir = Path(first_file_path).parent / "Merged sequence files"
        output_dir.mkdir(exist_ok=True)

        for i, (group, files) in enumerate(self.grouped_files.items()):
            if self.aborted:
                self.status.emit("Status: Aborted")
                self.update_status.emit("Merging aborted by user.")
                return
            
            self.status.emit(f"Merging group {group} ({i+1}/{total_groups})")
            r1_file, r2_file = files["R1"], files["R2"]
            # Extract base name before "_R1_001" or "_R2_001" pattern
            base_name = re.sub(r"_R[12]_001\.(fastq|fq)(\.gz)?$", "", os.path.basename(r1_file))
            output_file = output_dir / f"{base_name}.fastq"
            start_time = time.time()
            merged_count, reads_processed = self.merge_files(r1_file, r2_file, str(output_file))
            end_time = time.time()
            total_reads_processed += reads_processed
            
            finish_time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
            self.report.emit(
                f"Group {group}: Merged {merged_count} reads from {os.path.basename(r1_file)} and {os.path.basename(r2_file)}\n"
                f"Processing finished at: {finish_time_str}"
            )
            self.update_status.emit(f"Processed {total_reads_processed}/{total_reads} reads across all groups")
            self.progress.emit(int((total_reads_processed / total_reads) * 100))
            merged_files.append(str(output_file))

        if not self.aborted:
            self.status.emit("Status: Merge Complete")
            self.finished.emit(merged_files)
        else:
            self.status.emit("Status: Aborted")

    def merge_files(self, r1_file, r2_file, output_file):
        self.update_status.emit(f"Loading {os.path.basename(r1_file)} and {os.path.basename(r2_file)}...")
        r1_records = list(SeqIO.parse(gzip.open(r1_file, "rt") if r1_file.endswith(".gz") else r1_file, "fastq"))
        r2_records = list(SeqIO.parse(gzip.open(r2_file, "rt") if r2_file.endswith(".gz") else r2_file, "fastq"))
        total_reads = min(len(r1_records), len(r2_records))
        merged_count = 0
        reads_processed = 0

        self.update_status.emit(f"Merging {total_reads} read pairs from {os.path.basename(r1_file)} and {os.path.basename(r2_file)}")

        with open(output_file, "w") as output_handle:
            for r1, r2 in zip(r1_records, r2_records):
                if self.aborted:
                    break
                reads_processed += 1
                seq1, seq2 = str(r1.seq), str(r2.seq.reverse_complement())
                overlap = self.find_overlap(seq1, seq2)
                if overlap >= self.min_overlap:
                    merged_seq = seq1 + seq2[overlap:]
                    merged_qual = r1.letter_annotations["phred_quality"] + r2.letter_annotations["phred_quality"][overlap:]
                    record = SeqIO.SeqRecord(Seq(merged_seq), id=r1.id, description="", letter_annotations={"phred_quality": merged_qual})
                    SeqIO.write(record, output_handle, "fastq")
                    merged_count += 1
                
                if reads_processed % 1000 == 0:
                    self.update_status.emit(f"Merging {os.path.basename(r1_file)}: {reads_processed}/{total_reads} reads processed ({merged_count} merged)")
                    QApplication.processEvents()

        self.update_status.emit(f"Completed merging {os.path.basename(r1_file)}: {merged_count} reads merged out of {total_reads}")
        return merged_count, reads_processed

    def find_overlap(self, seq1, seq2):
        for i in range(min(len(seq1), len(seq2)), self.min_overlap - 1, -1):
            if seq1[-i:] == seq2[:i]:
                return i
        return 0

    def abort(self):
        self.aborted = True

class PairedEndMerger(QWidget):
    merge_finished = pyqtSignal(list)
    
    def __init__(self):
        super().__init__()
        self.files = []
        self.grouped_files = {}
        self.min_overlap = 20
        self.current_language = "en"
        self.worker = None
        self.thread = None
        self.translations = {
            "en": {
                "step1": "Step 1: Choose minimal overlap length",
                "selected_files": "Selected Files: None",
                "selected_files_count": "Selected {0} files",
                "step2": "Step 2: Select paired files",
                "group": "Group",
                "r1_file": "R1 File",
                "r2_file": "R2 File",
                "status_idle": "Status: Idle",
                "step3": "Step 3: Merge files",
                "abort": "Abort"
            },
            "zh_TW": {
                "step1": "步驟一：選擇最小重疊長度",
                "selected_files": "已選檔案：無",
                "selected_files_count": "已選 {0} 個檔案",
                "step2": "步驟二：選擇配對檔案",
                "group": "組別",
                "r1_file": "R1 檔案",
                "r2_file": "R2 檔案",
                "status_idle": "狀態：閒置",
                "step3": "步驟三：合併檔案",
                "abort": "中止"
            },
            "zh_CN": {
                "step1": "步骤一：选择最小重叠长度",
                "selected_files": "已选文件：无",
                "selected_files_count": "已选 {0} 个文件",
                "step2": "步骤二：选择配对文件",
                "group": "组别",
                "r1_file": "R1 文件",
                "r2_file": "R2 文件",
                "status_idle": "状态：空闲",
                "step3": "步骤三：合并文件",
                "abort": "中止"
            },
            "ja": {
                "step1": "ステップ1：最小オーバーラップ長を選択",
                "selected_files": "選択されたファイル：なし",
                "selected_files_count": "{0} ファイルを選択",
                "step2": "ステップ2：ペアファイルを選択",
                "group": "グループ",
                "r1_file": "R1 ファイル",
                "r2_file": "R2 ファイル",
                "status_idle": "状態：待機中",
                "step3": "ステップ3：ファイルを結合",
                "abort": "中止"
            },
            "ko": {
                "step1": "1단계: 최소 오버랩 길이 선택",
                "selected_files": "선택된 파일: 없음",
                "selected_files_count": "{0} 파일 선택됨",
                "step2": "2단계: 페어 파일 선택",
                "group": "그룹",
                "r1_file": "R1 파일",
                "r2_file": "R2 파일",
                "status_idle": "상태: 대기 중",
                "step3": "3단계: 파일 병합",
                "abort": "중단"
            },
            "fr": {
                "step1": "Étape 1 : Choisir la longueur minimale de chevauchement",
                "selected_files": "Fichiers sélectionnés : Aucun",
                "selected_files_count": "{0} fichiers sélectionnés",
                "step2": "Étape 2 : Sélectionner les fichiers appariés",
                "group": "Groupe",
                "r1_file": "Fichier R1",
                "r2_file": "Fichier R2",
                "status_idle": "Statut : Inactif",
                "step3": "Étape 3 : Fusionner les fichiers",
                "abort": "Annuler"
            },
            "de": {
                "step1": "Schritt 1: Minimale Überlappungslänge wählen",
                "selected_files": "Ausgewählte Dateien: Keine",
                "selected_files_count": "{0} Dateien ausgewählt",
                "step2": "Schritt 2: Gepairte Dateien auswählen",
                "group": "Gruppe",
                "r1_file": "R1 Datei",
                "r2_file": "R2 Datei",
                "status_idle": "Status: Inaktiv",
                "step3": "Schritt 3: Dateien zusammenführen",
                "abort": "Abbrechen"
            },
            "es": {
                "step1": "Paso 1: Elegir longitud mínima de solapamiento",
                "selected_files": "Archivos seleccionados: Ninguno",
                "selected_files_count": "{0} archivos seleccionados",
                "step2": "Paso 2: Seleccionar archivos emparejados",
                "group": "Grupo",
                "r1_file": "Archivo R1",
                "r2_file": "Archivo R2",
                "status_idle": "Estado: Inactivo",
                "step3": "Paso 3: Fusionar archivos",
                "abort": "Abortar"
            },
            "pt": {
                "step1": "Passo 1: Escolher comprimento mínimo de sobreposição",
                "selected_files": "Arquivos selecionados: Nenhum",
                "selected_files_count": "{0} arquivos selecionados",
                "step2": "Passo 2: Selecionar arquivos pareados",
                "group": "Grupo",
                "r1_file": "Arquivo R1",
                "r2_file": "Arquivo R2",
                "status_idle": "Status: Inativo",
                "step3": "Passo 3: Mesclar arquivos",
                "abort": "Abortar"
            }
        }
        self.setup_ui()
    
    def setup_ui(self):
        layout = QVBoxLayout()
        
        overlap_layout = QHBoxLayout()
        self.step1_label = QLabel(self.translations["en"]["step1"])
        self.step1_label.setStyleSheet("color: #007BFF;")
        self.overlap_combo = QComboBox()
        self.overlap_combo.addItems(["10", "20", "50"])
        self.overlap_combo.setCurrentText("20")
        self.overlap_combo.setFixedWidth(60)
        self.overlap_combo.currentTextChanged.connect(self.update_overlap)
        overlap_layout.addWidget(self.step1_label)
        overlap_layout.addWidget(self.overlap_combo)
        overlap_layout.addStretch()
        layout.addLayout(overlap_layout)
        
        file_layout = QHBoxLayout()
        self.file_label = QLabel(self.translations["en"]["selected_files"])
        self.select_file_btn = QPushButton(self.translations["en"]["step2"])
        self.select_file_btn.setStyleSheet("QPushButton { background-color: #007BFF; color: white; } QPushButton:hover { color: blue; }")
        self.select_file_btn.clicked.connect(self.select_files)
        file_layout.addWidget(self.file_label)
        file_layout.addWidget(self.select_file_btn)
        
        self.file_table = QTableWidget()
        self.file_table.setColumnCount(3)
        self.file_table.setHorizontalHeaderLabels([self.translations["en"]["group"], self.translations["en"]["r1_file"], self.translations["en"]["r2_file"]])
        self.file_table.setColumnWidth(0, 200)
        self.file_table.setColumnWidth(1, 300)
        self.file_table.setColumnWidth(2, 300)
        layout.addWidget(self.file_table)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setMinimum(0)
        self.progress_bar.setMaximum(100)
        self.progress_bar.setVisible(False)
        
        self.status_label = QLabel(self.translations["en"]["status_idle"])
        self.detailed_status_label = QLabel("")
        self.progress_info_label = QLabel("")
        
        button_layout = QHBoxLayout()
        self.merge_btn = QPushButton(self.translations["en"]["step3"])
        self.merge_btn.setStyleSheet("QPushButton { background-color: #28A745; color: white; } QPushButton:hover { color: #28A745; }")
        self.merge_btn.clicked.connect(self.start_merge)
        self.merge_btn.setEnabled(False)
        
        self.abort_btn = QPushButton(self.translations["en"]["abort"])
        self.abort_btn.setStyleSheet("QPushButton { background-color: #FFC1CC; color: black; } QPushButton:hover { color: #FF0000; }")
        self.abort_btn.clicked.connect(self.abort_merge)
        self.abort_btn.setEnabled(False)
        
        button_layout.addWidget(self.merge_btn)
        button_layout.addWidget(self.abort_btn)
        
        self.merge_report = QTextEdit()
        self.merge_report.setReadOnly(True)
        layout.addWidget(self.merge_report)
        
        layout.addLayout(file_layout)
        layout.addWidget(self.progress_bar)
        layout.addWidget(self.status_label)
        layout.addWidget(self.detailed_status_label)
        layout.addWidget(self.progress_info_label)
        layout.addLayout(button_layout)
        
        self.setLayout(layout)
        self.set_language(self.current_language)

    def set_language(self, language):
        self.current_language = language
        self.translations[self.current_language] = self.translations.get(self.current_language, self.translations["en"])
        t = self.translations[self.current_language]
        self.step1_label.setText(t["step1"])
        self.file_label.setText(t["selected_files"])
        self.select_file_btn.setText(t["step2"])
        self.file_table.setHorizontalHeaderLabels([t["group"], t["r1_file"], t["r2_file"]])
        self.status_label.setText(t["status_idle"])
        self.merge_btn.setText(t["step3"])
        self.abort_btn.setText(t["abort"])

    def update_overlap(self, value):
        self.min_overlap = int(value)

    def select_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select Paired FASTQ Files", "", "FASTQ Files (*.fastq *.fq *.fastq.gz *.fq.gz)")
        if files:
            self.files = files
            self.group_files()
            self.update_file_table()
            self.file_label.setText(self.translations[self.current_language]["selected_files_count"].format(len(self.files)))
            self.merge_btn.setEnabled(True)

    def group_files(self):
        self.grouped_files.clear()
        for file in self.files:
            file_name = os.path.basename(file)
            group = re.sub(r"(_R[12]).*$", "", file_name)
            if group not in self.grouped_files:
                self.grouped_files[group] = {"R1": None, "R2": None}
            if "_R1" in file_name:
                self.grouped_files[group]["R1"] = file
            elif "_R2" in file_name:
                self.grouped_files[group]["R2"] = file

    def update_file_table(self):
        self.file_table.setRowCount(0)
        for group, files in self.grouped_files.items():
            row = self.file_table.rowCount()
            self.file_table.insertRow(row)
            self.file_table.setItem(row, 0, QTableWidgetItem(group))
            self.file_table.setItem(row, 1, QTableWidgetItem(os.path.basename(files["R1"]) if files["R1"] else "N/A"))
            self.file_table.setItem(row, 2, QTableWidgetItem(os.path.basename(files["R2"]) if files["R2"] else "N/A"))

    def start_merge(self):
        if not self.grouped_files:
            QMessageBox.warning(self, "Error", "Please select paired files first.")
            return

        self.progress_bar.setValue(0)
        self.progress_bar.setVisible(True)
        self.merge_btn.setEnabled(False)
        self.abort_btn.setEnabled(True)
        self.merge_report.clear()
        self.detailed_status_label.setText("Starting merge process...")

        self.thread = QThread()
        self.worker = MergeWorker(self.grouped_files, self.min_overlap)
        self.worker.moveToThread(self.thread)
        self.thread.started.connect(self.worker.run)
        self.worker.progress.connect(self.update_progress)
        self.worker.finished.connect(self.on_merge_finished)
        self.worker.status.connect(self.status_label.setText)
        self.worker.update_status.connect(self.detailed_status_label.setText)
        self.worker.report.connect(self.merge_report.append)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)
        self.thread.start()

    def update_progress(self, value):
        self.progress_bar.setValue(value)
        QApplication.processEvents()

    def on_merge_finished(self, merged_files):
        self.merge_finished.emit(merged_files)
        self.progress_bar.setVisible(False)
        self.merge_btn.setEnabled(True)
        self.abort_btn.setEnabled(False)
        self.detailed_status_label.setText("Merge process completed.")
        QMessageBox.information(self, "Success", "Merging completed successfully.")

    def abort_merge(self):
        if self.worker:
            self.worker.abort()
            self.status_label.setText("Status: Aborting...")
            self.detailed_status_label.setText("Aborting merge process...")
            self.abort_btn.setEnabled(False)