from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
    QPushButton, QTableWidget, QTableWidgetItem, 
    QFileDialog, QTextEdit, QProgressBar, QMessageBox, QRadioButton, QButtonGroup, QApplication
)
from PyQt6.QtCore import pyqtSignal
from Bio import SeqIO
from Bio.Seq import Seq
from pathlib import Path
import os
import pandas as pd
import time
import gzip
import re

class SequenceExtractor(QWidget):
    extraction_finished = pyqtSignal(list)
    
    def __init__(self):
        super().__init__()
        self.files = []
        self.barcodes = []
        self.primer = ""
        self.project_name = ""
        self.current_language = "en"
        self.aborted = False
        self.barcode_type = "replicates"  # Default
        self.barcode_to_label = {}  # Maps barcodes to "RepX" or group names
        
        self.translations = {
            "current_project": "Current Project: None",
            "selected_files": "Selected Files: None",
            "selected_files_count": "Selected {0} files",
            "step2": "Step 2: Select files",
            "file_name": "File Name",
            "total_reads": "Total Reads",
            "step3": "Step 3: Extract Sequences",
            "abort": "Abort",
            "status_idle": "Status: Idle",
            "status_loading": "Status: Loading file {0} ({1}/{2})",
            "step1": "Step 1: Choose extraction mode",
            "both_directions": "Extract both forward and reverse (Recommended for NGS using Nanopore)",
            "forward_only": "Extract from one direction only (Recommended for NGS using Illumina. Use either R1+R2 merged file or R1 alone)",
            "extraction_complete": "Extraction complete"
        }
        
        self.setup_ui()
        
    def setup_ui(self):
        layout = QVBoxLayout()
        
        self.step1_label = QLabel(self.translations["step1"])
        self.step1_label.setStyleSheet("color: #007BFF;")
        layout.addWidget(self.step1_label)
        
        direction_layout = QHBoxLayout()
        self.both_directions_radio = QRadioButton(self.translations["both_directions"])
        self.forward_only_radio = QRadioButton(self.translations["forward_only"])
        self.both_directions_radio.setChecked(True)
        direction_group = QButtonGroup(self)
        direction_group.addButton(self.both_directions_radio)
        direction_group.addButton(self.forward_only_radio)
        direction_group.setExclusive(True)
        direction_layout.addWidget(self.both_directions_radio)
        direction_layout.addWidget(self.forward_only_radio)
        layout.addLayout(direction_layout)
        
        self.project_label = QLabel(self.translations["current_project"])
        layout.addWidget(self.project_label)
        
        file_layout = QHBoxLayout()
        self.file_label = QLabel(self.translations["selected_files"])
        self.select_file_btn = QPushButton(self.translations["step2"])
        self.select_file_btn.setStyleSheet("""
            QPushButton { background-color: #007BFF; color: white; }
            QPushButton:hover { color: blue; }
        """)
        self.select_file_btn.clicked.connect(self.select_files)
        file_layout.addWidget(self.file_label)
        file_layout.addWidget(self.select_file_btn)
        
        self.file_table = QTableWidget()
        self.file_table.setColumnCount(2)
        self.file_table.setHorizontalHeaderLabels([self.translations["file_name"], self.translations["total_reads"]])
        self.file_table.setColumnWidth(0, 600)
        self.file_table.setColumnWidth(1, 300)
        layout.addWidget(self.file_table)
        
        self.progress_bar = QProgressBar()
        self.progress_bar.setVisible(False)
        
        self.status_label = QLabel(self.translations["status_idle"])
        self.progress_info_label = QLabel("")
        
        button_layout = QHBoxLayout()
        self.extract_btn = QPushButton(self.translations["step3"])
        self.extract_btn.setStyleSheet("""
            QPushButton { background-color: #28A745; color: white; }
            QPushButton:hover { color: #28A745; }
        """)
        self.extract_btn.clicked.connect(self.start_extraction)
        self.extract_btn.setEnabled(False)
        
        self.abort_btn = QPushButton(self.translations["abort"])
        self.abort_btn.setStyleSheet("""
            QPushButton { background-color: #FFC1CC; color: black; }
            QPushButton:hover { color: #FF0000; }
        """)
        self.abort_btn.clicked.connect(self.abort_extraction)
        self.abort_btn.setEnabled(False)
        
        button_layout.addWidget(self.extract_btn)
        button_layout.addWidget(self.abort_btn)
        
        self.extraction_report = QTextEdit()
        self.extraction_report.setReadOnly(True)
        
        layout.addLayout(file_layout)
        layout.addWidget(self.progress_bar)
        layout.addWidget(self.status_label)
        layout.addWidget(self.progress_info_label)
        layout.addLayout(button_layout)
        layout.addWidget(self.extraction_report)
        
        self.setLayout(layout)
    
    def set_language(self, language):
        self.current_language = language
        translations = {
            "en": {
                "current_project": "Current Project: None",
                "selected_files": "Selected Files: None",
                "selected_files_count": "Selected {0} files",
                "step2": "Step 2: Select files",
                "file_name": "File Name",
                "total_reads": "Total Reads",
                "step3": "Step 3: Extract Sequences",
                "abort": "Abort",
                "status_idle": "Status: Idle",
                "status_loading": "Status: Loading file {0} ({1}/{2})",
                "step1": "Step 1: Choose extraction mode",
                "both_directions": "Extract both forward and reverse (Recommended for NGS using Nanopore)",
                "forward_only": "Extract from one direction only (Recommended for NGS using Illumina. Use either R1+R2 merged file or R1 alone)",
                "extraction_complete": "Extraction complete"
            },
            "zh_TW": {
                "current_project": "當前專案：無",
                "selected_files": "已選檔案：無",
                "selected_files_count": "已選 {0} 個檔案",
                "step2": "步驟二：選擇檔案",
                "file_name": "檔案名稱",
                "total_reads": "總讀數",
                "step3": "步驟三：提取序列",
                "abort": "中止",
                "status_idle": "狀態：閒置",
                "status_loading": "狀態：正在載入檔案 {0} ({1}/{2})",
                "step1": "步驟一：選擇提取模式",
                "both_directions": "從正向和反向提取（推薦用於使用Nanopore的NGS）",
                "forward_only": "僅從正向提取（推薦用於使用Illumina的NGS，僅需R1或R1+R2合併的檔案）",
                "extraction_complete": "提取完成"
            },
            "zh_CN": {
                "current_project": "当前项目：无",
                "selected_files": "已选文件：无",
                "selected_files_count": "已选 {0} 个文件",
                "step2": "步骤二：选择文件",
                "file_name": "文件名",
                "total_reads": "总读数",
                "step3": "步骤三：提取序列",
                "abort": "中止",
                "status_idle": "状态：空闲",
                "status_loading": "状态：正在加载文件 {0} ({1}/{2})",
                "step1": "步骤一：选择提取模式",
                "both_directions": "从正向和反向提取（推荐用于使用Nanopore的NGS）",
                "forward_only": "仅从正向提取（推荐用于使用Illumina的NGS，仅需R1文件）",
                "extraction_complete": "提取完成"
            },
            "ja": {
                "current_project": "現在のプロジェクト：なし",
                "selected_files": "選択されたファイル：なし",
                "selected_files_count": "{0} 個のファイルを選択",
                "step2": "ステップ2：ファイルを選択",
                "file_name": "ファイル名",
                "total_reads": "総リード数",
                "step3": "ステップ3：シーケンスを抽出",
                "abort": "中止",
                "status_idle": "状態：待機中",
                "status_loading": "状態：ファイル {0} を読み込み中 ({1}/{2})",
                "step1": "ステップ1：抽出モードを選択",
                "both_directions": "正方向と逆方向の両方を抽出（Nanoporeを使用したNGSに推奨）",
                "forward_only": "一方向のみ抽出（Illuminaを使用したNGSに推奨。R1+R2統合ファイルまたはR1のみを使用）",
                "extraction_complete": "抽出完了"
            },
            "ko": {
                "current_project": "현재 프로젝트: 없음",
                "selected_files": "선택된 파일: 없음",
                "selected_files_count": "{0} 개 파일 선택됨",
                "step2": "2단계: 파일 선택",
                "file_name": "파일 이름",
                "total_reads": "총 리드 수",
                "step3": "3단계: 시퀀스 추출",
                "abort": "중단",
                "status_idle": "상태: 대기 중",
                "status_loading": "상태: 파일 {0} 로드 중 ({1}/{2})",
                "step1": "1단계: 추출 모드 선택",
                "both_directions": "정방향과 역방향 모두 추출 (Nanopore를 사용한 NGS에 권장)",
                "forward_only": "한 방향만 추출 (Illumina를 사용한 NGS에 권장. R1+R2 병합 파일 또는 R1 단독 사용)",
                "extraction_complete": "추출 완료"
            },
            "fr": {
                "current_project": "Projet actuel : Aucun",
                "selected_files": "Fichiers sélectionnés : Aucun",
                "selected_files_count": "{0} fichiers sélectionnés",
                "step2": "Étape 2 : Sélectionner des fichiers",
                "file_name": "Nom du fichier",
                "total_reads": "Lectures totales",
                "step3": "Étape 3 : Extraire les séquences",
                "abort": "Abandonner",
                "status_idle": "Statut : Inactif",
                "status_loading": "Statut : Chargement du fichier {0} ({1}/{2})",
                "step1": "Étape 1 : Choisir le mode d'extraction",
                "both_directions": "Extraire dans les deux sens (Recommandé pour NGS avec Nanopore)",
                "forward_only": "Extraire dans une seule direction (Recommandé pour NGS avec Illumina. Utiliser un fichier R1+R2 fusionné ou R1 seul)",
                "extraction_complete": "Extraction terminée"
            },
            "de": {
                "current_project": "Aktuelles Projekt: Keines",
                "selected_files": "Ausgewählte Dateien: Keine",
                "selected_files_count": "{0} Dateien ausgewählt",
                "step2": "Schritt 2: Dateien auswählen",
                "file_name": "Dateiname",
                "total_reads": "Gesamte Reads",
                "step3": "Schritt 3: Sequenzen extrahieren",
                "abort": "Abbrechen",
                "status_idle": "Status: Inaktiv",
                "status_loading": "Status: Lade Datei {0} ({1}/{2})",
                "step1": "Schritt 1: Extraktionsmodus wählen",
                "both_directions": "In beide Richtungen extrahieren (Empfohlen für NGS mit Nanopore)",
                "forward_only": "Nur in eine Richtung extrahieren (Empfohlen für NGS mit Illumina. Verwenden Sie entweder eine R1+R2 kombinierte Datei oder nur R1)",
                "extraction_complete": "Extraktion abgeschlossen"
            },
            "es": {
                "current_project": "Proyecto actual: Ninguno",
                "selected_files": "Archivos seleccionados: Ninguno",
                "selected_files_count": "{0} archivos seleccionados",
                "step2": "Paso 2: Seleccionar archivos",
                "file_name": "Nombre del archivo",
                "total_reads": "Lecturas totales",
                "step3": "Paso 3: Extraer secuencias",
                "abort": "Abortar",
                "status_idle": "Estado: Inactivo",
                "status_loading": "Estado: Cargando archivo {0} ({1}/{2})",
                "step1": "Paso 1: Elegir modo de extracción",
                "both_directions": "Extraer en ambas direcciones (Recomendado para NGS con Nanopore)",
                "forward_only": "Extraer en una sola dirección (Recomendado para NGS con Illumina. Use un archivo combinado R1+R2 o solo R1)",
                "extraction_complete": "Extracción completada"
            },
            "pt": {
                "current_project": "Projeto atual: Nenhum",
                "selected_files": "Arquivos selecionados: Nenhum",
                "selected_files_count": "{0} arquivos selecionados",
                "step2": "Passo 2: Selecionar arquivos",
                "file_name": "Nome do arquivo",
                "total_reads": "Leituras totais",
                "step3": "Passo 3: Extrair sequências",
                "abort": "Abortar",
                "status_idle": "Status: Inativo",
                "status_loading": "Status: Carregando arquivo {0} ({1}/{2})",
                "step1": "Passo 1: Escolher modo de extração",
                "both_directions": "Extrair em ambas as direções (Recomendado para NGS com Nanopore)",
                "forward_only": "Extrair em uma única direção (Recomendado para NGS com Illumina. Use um arquivo mesclado R1+R2 ou apenas R1)",
                "extraction_complete": "Extração concluída"
            }
        }
        
        self.translations = translations.get(self.current_language, translations["en"])
        self.project_label.setText(self.translations["current_project"])
        self.file_label.setText(self.translations["selected_files"])
        self.select_file_btn.setText(self.translations["step2"])
        self.file_table.setHorizontalHeaderLabels([self.translations["file_name"], self.translations["total_reads"]])
        self.extract_btn.setText(self.translations["step3"])
        self.abort_btn.setText(self.translations["abort"])
        self.status_label.setText(self.translations["status_idle"])
        self.step1_label.setText(self.translations["step1"])
        self.both_directions_radio.setText(self.translations["both_directions"])
        self.forward_only_radio.setText(self.translations["forward_only"])
    
    def select_files(self):
        files, _ = QFileDialog.getOpenFileNames(self, "Select Files", "", "FASTQ Files (*.fastq *.fastq.gz)")
        if files:
            self.files = files
            self.file_label.setText(self.translations["selected_files_count"].format(len(self.files)))
            self.file_table.setRowCount(0)
            
            for idx, file in enumerate(self.files, 1):
                file_name = os.path.basename(file)
                self.status_label.setText(self.translations["status_loading"].format(file_name, idx, len(self.files)))
                QApplication.processEvents()
                
                total_reads = self.get_read_count(file)
                row = self.file_table.rowCount()
                self.file_table.insertRow(row)
                self.file_table.setItem(row, 0, QTableWidgetItem(file_name))
                self.file_table.setItem(row, 1, QTableWidgetItem(str(total_reads)))
                QApplication.processEvents()
            
            self.status_label.setText(self.translations["status_idle"])
            self.extract_btn.setEnabled(True)
    
    def get_read_count(self, file):
        try:
            opener = gzip.open if file.endswith('.gz') else open
            with opener(file, "rt") as handle:
                count = sum(1 for _ in SeqIO.parse(handle, "fastq"))
            print(f"Counted {count} records in {file}")
            return count
        except Exception as e:
            print(f"Error counting records in {file}: {str(e)}")
            return 0
    
    def set_sequences(self, barcodes, primer, project_name):
        self.barcodes = barcodes if barcodes else []
        self.primer = primer
        self.project_name = project_name
        self.project_label.setText(self.translations["current_project"].replace("None", project_name))
    
    def set_project(self, project_info):
        self.project_name = project_info.get('name', '')
        sequences = project_info.get('sequences', [])
        self.barcodes = [seq.get('barcode', '') for seq in sequences if seq.get('barcode', '').strip()]
        self.primer = sequences[0].get('primer', '') if sequences else ''
        self.barcode_type = project_info.get('barcode_type', 'replicates')
        self.barcode_to_label = {}
        for seq in sequences:
            barcode = seq.get('barcode', '').strip()
            if barcode:
                if self.barcode_type == "replicates":
                    replicate = seq.get('replicate', '')
                    if replicate:
                        match = re.match(r"Replicate\s+(\d+)", replicate)
                        if match:
                            rep_number = match.group(1)
                            self.barcode_to_label[barcode] = f"Rep{rep_number}"
                        else:
                            self.barcode_to_label[barcode] = replicate
                else:  # groups
                    group = seq.get('group', '')
                    self.barcode_to_label[barcode] = group
        self.project_label.setText(self.translations["current_project"].replace("None", self.project_name))
    
    def start_extraction(self):
        if not self.files:
            QMessageBox.warning(self, "Warning", "Please select files first")
            return
        
        if not self.primer:
            QMessageBox.warning(self, "Warning", "Please set primer sequence first")
            return
        
        self.aborted = False
        self.progress_bar.setVisible(True)
        self.extract_btn.setEnabled(False)
        self.abort_btn.setEnabled(True)
        self.status_label.setText(f"Status: Starting extraction for {len(self.files)} files...")
        self.progress_info_label.setText(f"Processing 0/{len(self.files)} files")
        QApplication.processEvents()
        
        output_dir = Path(self.files[0]).parent / "Barcode-segregated files"
        output_dir.mkdir(exist_ok=True)
        
        extraction_report = {}
        extracted_files = []
        start_time = time.time()
        
        for file_idx, input_file in enumerate(self.files, 1):
            if self.aborted:
                self.status_label.setText("Status: Extraction aborted by user.")
                break
            
            elapsed_time = time.time() - start_time
            files_processed = file_idx - 1
            files_remaining = len(self.files) - files_processed
            if files_processed > 0:
                time_per_file = elapsed_time / files_processed
                est_remaining_time = time_per_file * files_remaining
                est_minutes = int(est_remaining_time // 60)
                est_seconds = int(est_remaining_time % 60)
                time_str = f" (Est. remaining: {est_minutes}m {est_seconds}s)"
            else:
                time_str = ""
            
            self.progress_info_label.setText(f"Processing {file_idx}/{len(self.files)} files{time_str}")
            self.status_label.setText(f"Status: Processing file {file_idx} of {len(self.files)}: {os.path.basename(input_file)}")
            self.progress_bar.setValue(int((file_idx - 1) / len(self.files) * 100))
            QApplication.processEvents()
            
            extracted_files.extend(self.split_fastq(input_file, self.barcodes, self.primer, output_dir, extraction_report))
        
        if not self.aborted:
            self.progress_bar.setValue(100)
            self.status_label.setText(f"Status: Extraction completed for {len(self.files)} files.")
            report_message = self.translations["extraction_complete"] + "\n\n"
            for input_file, counts in extraction_report.items():
                report_message += f"File {os.path.basename(input_file)}:\n"
                for key, count_info in counts.items():
                    report_message += f"{' ' * 40}{'Barcode ' + key if self.barcodes else 'Filtered'}: {count_info['total']} sequences\n"
                    report_message += f"{' ' * 60}Forward: {count_info['forward']} sequences\n"
                    if self.both_directions_radio.isChecked():
                        report_message += f"{' ' * 60}Opposite: {count_info['opposite']} sequences\n"
            self.extraction_report.setText(report_message)
            self.save_extracted_sequence_info(extracted_files)
            self.extraction_finished.emit(extracted_files)
        else:
            self.extraction_report.setText("Extraction was aborted.")
            self.progress_info_label.setText(f"Processing {file_idx-1}/{len(self.files)} files (Aborted)")
        
        self.progress_bar.setVisible(False)
        self.extract_btn.setEnabled(True)
        self.abort_btn.setEnabled(False)
    
    def abort_extraction(self):
        self.aborted = True
        self.status_label.setText("Status: Aborting extraction... Please wait for current file to finish.")
        QApplication.processEvents()
    
    def reverse_complement(self, sequence):
        return str(Seq(sequence).reverse_complement())
    
    def split_fastq(self, input_file, barcodes, primer, output_folder, extraction_report):
        extracted_files = []
        
        if barcodes:
            extraction_report[input_file] = {barcode: {'total': 0, 'forward': 0, 'opposite': 0} for barcode in barcodes}
        else:
            extraction_report[input_file] = {'filtered': {'total': 0, 'forward': 0, 'opposite': 0}}
        
        try:
            input_dir = Path(input_file).parent
            input_base = Path(input_file).stem
            if barcodes:
                extracted_records = {barcode: [] for barcode in barcodes}
            else:
                extracted_records = {'filtered': []}
            
            opener = gzip.open if input_file.endswith('.gz') else open
            
            if self.both_directions_radio.isChecked():
                rc_dir = input_dir / "ReverseComplement"
                rc_dir.mkdir(exist_ok=True)
                rc_file = rc_dir / f"{input_base}-ReverseComplement.fastq"
                
                self.status_label.setText(f"Status: Generating reverse complement for {os.path.basename(input_file)}")
                QApplication.processEvents()
                
                records_written = 0
                with opener(input_file, "rt") as in_handle, open(rc_file, "w") as rc_handle:
                    for record in SeqIO.parse(in_handle, "fastq"):
                        if self.aborted:
                            break
                        rc_seq = self.reverse_complement(str(record.seq))
                        record.seq = Seq(rc_seq)
                        if 'phred_quality' in record.letter_annotations:
                            record.letter_annotations['phred_quality'] = record.letter_annotations['phred_quality'][::-1]
                        else:
                            print(f"Warning: No quality scores for record {record.id}")
                        SeqIO.write(record, rc_handle, "fastq")
                        records_written += 1
                        QApplication.processEvents()
                print(f"Wrote {records_written} reverse complement records to {rc_file}")
                
                if not self.aborted and not rc_file.exists():
                    raise FileNotFoundError(f"Failed to create {rc_file}")
                file_size = rc_file.stat().st_size if rc_file.exists() else 0
                print(f"Verified {rc_file} exists with size: {file_size} bytes")
                if file_size == 0:
                    raise ValueError(f"{rc_file} is empty")
                
                total_records = self.get_read_count(input_file) + (self.get_read_count(str(rc_file)) if not self.aborted else 0)
                if total_records == 0:
                    print(f"Warning: No records found in {input_file} or {rc_file}")
                    return extracted_files
                
                processed_records = 0
                
                self.status_label.setText(f"Status: Screening forward reads for {os.path.basename(input_file)}")
                QApplication.processEvents()
                
                with opener(input_file, "rt") as handle:
                    for record in SeqIO.parse(handle, "fastq"):
                        if self.aborted:
                            break
                        seq = str(record.seq)
                        if barcodes:
                            for barcode in barcodes:
                                pattern = barcode + primer
                                if pattern in seq:
                                    extracted_records[barcode].append(record)
                                    extraction_report[input_file][barcode]['forward'] += 1
                                    extraction_report[input_file][barcode]['total'] += 1
                                    print(f"Match for barcode {barcode} in forward: {seq[:20]}...")
                        else:
                            if primer in seq:
                                extracted_records['filtered'].append(record)
                                extraction_report[input_file]['filtered']['forward'] += 1
                                extraction_report[input_file]['filtered']['total'] += 1
                                print(f"Match for primer only in forward: {seq[:20]}...")
                        
                        processed_records += 1
                        self.progress_bar.setValue(int((processed_records / total_records) * 100))
                        QApplication.processEvents()
                
                if not self.aborted:
                    self.status_label.setText(f"Status: Screening reverse complement reads for {os.path.basename(input_file)}")
                    QApplication.processEvents()
                    
                    with open(rc_file, "r") as handle:
                        for record in SeqIO.parse(handle, "fastq"):
                            if self.aborted:
                                break
                            seq = str(record.seq)
                            if barcodes:
                                for barcode in barcodes:
                                    pattern = barcode + primer
                                    if pattern in seq:
                                        extracted_records[barcode].append(record)
                                        extraction_report[input_file][barcode]['opposite'] += 1
                                        extraction_report[input_file][barcode]['total'] += 1
                                        print(f"Match for barcode {barcode} in opposite: {seq[:20]}...")
                            else:
                                if primer in seq:
                                    extracted_records['filtered'].append(record)
                                    extraction_report[input_file]['filtered']['opposite'] += 1
                                    extraction_report[input_file]['filtered']['total'] += 1
                                    print(f"Match for primer only in opposite: {seq[:20]}...")
                            
                            processed_records += 1
                            self.progress_bar.setValue(int((processed_records / total_records) * 100))
                            QApplication.processEvents()
                
                if self.aborted and rc_file.exists():
                    rc_file.unlink()
                    print(f"Cleaned up {rc_file} after abort")
            else:
                total_records = self.get_read_count(input_file)
                if total_records == 0:
                    print(f"Warning: No records found in {input_file}")
                    return extracted_files
                
                processed_records = 0
                
                self.status_label.setText(f"Status: Screening forward reads for {os.path.basename(input_file)}")
                QApplication.processEvents()
                
                with opener(input_file, "rt") as handle:
                    for record in SeqIO.parse(handle, "fastq"):
                        if self.aborted:
                            break
                        seq = str(record.seq)
                        if barcodes:
                            for barcode in barcodes:
                                pattern = barcode + primer
                                if pattern in seq:
                                    extracted_records[barcode].append(record)
                                    extraction_report[input_file][barcode]['forward'] += 1
                                    extraction_report[input_file][barcode]['total'] += 1
                                    print(f"Match for barcode {barcode} in forward: {seq[:20]}...")
                        else:
                            if primer in seq:
                                extracted_records['filtered'].append(record)
                                extraction_report[input_file]['filtered']['forward'] += 1
                                extraction_report[input_file]['filtered']['total'] += 1
                                print(f"Match for primer only in forward: {seq[:20]}...")
                        
                        processed_records += 1
                        self.progress_bar.setValue(int((processed_records / total_records) * 100))
                        QApplication.processEvents()
            
            if not self.aborted:
                self.status_label.setText(f"Status: Writing extracted reads for {os.path.basename(input_file)}")
                QApplication.processEvents()
                
                for key, records in extracted_records.items():
                    if records:
                        if key == 'filtered':
                            output_file = output_folder / f"{input_base}-filtered.fastq"
                        else:
                            if self.barcode_type == "replicates":
                                rep_label = self.barcode_to_label.get(key, "Unknown")
                                output_file = output_folder / f"{input_base}-[{rep_label}]-{key}.fastq"
                            else:  # groups
                                output_file = output_folder / f"{input_base}-{key}.fastq"
                        with open(output_file, "w") as handle:
                            SeqIO.write(records, handle, "fastq")
                        extracted_files.append(str(output_file))
                        print(f"Wrote {len(records)} records to {output_file}")
                    else:
                        print(f"No records extracted for {'barcode ' + key if barcodes else 'filtered'}")
        
        except Exception as e:
            print(f"Error processing {input_file}: {str(e)}")
            if self.both_directions_radio.isChecked() and 'rc_file' in locals() and rc_file.exists():
                rc_file.unlink()
                print(f"Cleaned up {rc_file} after error")
        
        return extracted_files
    
    def save_extracted_sequence_info(self, extracted_files):
        data = []
        for extracted_file in extracted_files:
            base_name = os.path.basename(extracted_file)
            data.append({"File Name": base_name})
        if data:
            temp_dir = Path("temp")
            temp_dir.mkdir(exist_ok=True)
            df = pd.DataFrame(data)
            df.to_csv(temp_dir / "ExtractedSequenceInfo.csv", index=False)