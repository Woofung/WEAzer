from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QLabel,
    QPushButton, QLineEdit, QFileDialog, QMessageBox,
    QTableWidget, QTableWidgetItem, QRadioButton, QStackedWidget, QFrame, QButtonGroup, QListWidget
)
from PyQt6.QtCore import pyqtSignal, Qt
from PyQt6.QtGui import QPalette, QColor
import json
import os
import re

class ProjectManager(QWidget):
    """
    ProjectManager widget for managing project settings, including barcode sequences and anchor configurations.
    """
    project_loaded = pyqtSignal(dict)
    anchors_updated = pyqtSignal(list)

    def __init__(self, sequence_extractor, editing_analyzer):
        super().__init__()
        self.sequence_extractor = sequence_extractor
        self.editing_analyzer = editing_analyzer
        self.current_language = "en"
        self.setup_ui()

    def browse_project_path(self):
        path = QFileDialog.getExistingDirectory(self, "Select Project Folder")
        if path:
            self.project_path.setText(path)

    def new_project(self):
        self.project_name.clear()
        self.project_path.clear()
        self.sequence_table.setRowCount(0)
        self.single_anchor_input.clear()
        self.single_anchor_input.setPlaceholderText("Left anchor sequence")
        self.left_anchor_input.clear()
        self.left_anchor_input.setPlaceholderText("Left anchor sequence")
        self.right_anchor_input.clear()
        self.right_anchor_input.setPlaceholderText("Right anchor sequence")
        self.range_input.clear()
        self.range_input.setPlaceholderText("Enter the range of nucleotides for analysis in bp (For example: 20)")
        self.range_input.setStyleSheet("")  # Reset style to default for enabled state
        self.reference_sequence_input.clear()  # Clear reference sequence input

    def validate_project_data(self):
        if not self.project_name.text().strip():
            QMessageBox.warning(self, "Warning", "Please enter a project name")
            return False
        if not self.project_path.text().strip():
            QMessageBox.warning(self, "Warning", "Please select a project path")
            return False
        if self.sequence_table.rowCount() == 0:
            QMessageBox.warning(self, "Warning", "Please add at least one sequence")
            return False
        for row in range(self.sequence_table.rowCount()):
            replicate_item = self.sequence_table.item(row, 0)
            primer_item = self.sequence_table.item(row, 2)
            if not (replicate_item and replicate_item.text().strip() and primer_item and primer_item.text().strip()):
                QMessageBox.warning(self, "Warning", "Please fill in both Replicate/Group and Primer Sequence")
                return False
        if self.flanking_anchor_radio.isChecked():
            if not (self.left_anchor_input.text().strip() and self.right_anchor_input.text().strip()):
                QMessageBox.warning(self, "Warning", "Please enter both Left and Right Anchor sequences")
                return False
        else:
            if not self.single_anchor_input.text().strip():
                QMessageBox.warning(self, "Warning", "Please enter an Anchor sequence")
                return False
            if not self.range_input.text().strip() or not self.range_input.text().isdigit() or int(self.range_input.text()) <= 0:
                QMessageBox.warning(self, "Warning", "Please enter a valid positive integer for Desired Analyzing Range")
                return False
        return True

    def save_project(self):
        if not self.validate_project_data():
            return
        project_data = {
            "name": self.project_name.text(),
            "path": self.project_path.text(),
            "sequences": [],
            "barcode_type": "replicates" if self.replicates_radio.isChecked() else "groups",
            "anchor_type": "flanking" if self.flanking_anchor_radio.isChecked() else "left" if self.left_anchor_radio.isChecked() else "right",
            "anchors": [],
            "range": self.range_input.text() if not self.flanking_anchor_radio.isChecked() else "",
            "reference_sequence": self.reference_sequence_input.text().strip().upper()  # Save reference sequence (uppercase)
        }
        for row in range(self.sequence_table.rowCount()):
            replicate_item = self.sequence_table.item(row, 0)
            barcode_item = self.sequence_table.item(row, 1)
            primer_item = self.sequence_table.item(row, 2)
            barcode_text = barcode_item.text().strip().upper() if barcode_item else ""  # Convert to uppercase
            primer_text = primer_item.text().strip().upper()  # Convert to uppercase
            project_data["sequences"].append({
                "replicate": replicate_item.text() if self.replicates_radio.isChecked() else "",
                "group": replicate_item.text() if not self.replicates_radio.isChecked() else "",
                "barcode": barcode_text,
                "primer": primer_text
            })
        if self.flanking_anchor_radio.isChecked():
            project_data["anchors"] = [self.left_anchor_input.text().strip().upper(), self.right_anchor_input.text().strip().upper()]  # Convert to uppercase
        else:
            project_data["anchors"] = [self.single_anchor_input.text().strip().upper()]  # Convert to uppercase

        project_file = os.path.join(self.project_path.text(), f"{self.project_name.text()}.json")
        with open(project_file, 'w') as f:
            json.dump(project_data, f, indent=4)

        # Pass barcodes as a list of strings to SequenceExtractor
        barcodes = [seq["barcode"] for seq in project_data["sequences"]]
        primer = project_data["sequences"][0]["primer"] if project_data["sequences"] else ""
        self.sequence_extractor.set_sequences(barcodes, primer, project_data["name"])
        self.editing_analyzer.set_project(project_data)
        self.project_loaded.emit(project_data)
        self.anchors_updated.emit(project_data["anchors"])
        # Update recent projects
        self.update_recent_projects(project_file)
        QMessageBox.information(self, "Success", "Project saved and applied successfully")

    def load_project(self):
        project_file, _ = QFileDialog.getOpenFileName(self, "Load Project", "", "JSON Files (*.json)")
        if project_file:
            self.load_project_file(project_file)

    def add_sequence(self):
        row = self.sequence_table.rowCount()
        self.sequence_table.insertRow(row)
        replicate_or_group = f"Replicate {row + 1}" if self.replicates_radio.isChecked() else f"Group {row + 1}"
        self.sequence_table.setItem(row, 0, QTableWidgetItem(replicate_or_group))
        barcode_item = QTableWidgetItem("(GGAA)")  # Uppercase default
        barcode_item.setForeground(QColor("grey"))
        self.sequence_table.setItem(row, 1, barcode_item)
        primer_item = QTableWidgetItem("(GTGATTGCTGGCCGGAACA)")  # Uppercase default
        primer_item.setForeground(QColor("grey"))
        self.sequence_table.setItem(row, 2, primer_item)

    def remove_sequence(self):
        selected_rows = self.sequence_table.selectionModel().selectedRows()
        for row in sorted(selected_rows, reverse=True):
            self.sequence_table.removeRow(row.row())

    def load_project_file(self, project_file):
        try:
            with open(project_file, 'r') as f:
                project_data = json.load(f)
            self.project_name.setText(project_data.get("name", ""))
            # Update project_path to the directory of the loaded file
            new_project_path = os.path.dirname(project_file)
            self.project_path.setText(new_project_path)
            # Update the project_data to reflect the new path
            project_data["path"] = new_project_path
            self.sequence_table.setRowCount(0)
            barcode_type = project_data.get("barcode_type", "replicates")
            self.replicates_radio.setChecked(barcode_type == "replicates")
            self.groups_radio.setChecked(barcode_type == "groups")
            self.update_sequence_table_headers()
            for seq in project_data.get("sequences", []):
                row = self.sequence_table.rowCount()
                self.sequence_table.insertRow(row)
                identifier = seq.get("replicate", "") if barcode_type == "replicates" else seq.get("group", "")
                self.sequence_table.setItem(row, 0, QTableWidgetItem(identifier))
                self.sequence_table.setItem(row, 1, QTableWidgetItem(seq.get("barcode", "")))
                self.sequence_table.setItem(row, 2, QTableWidgetItem(seq.get("primer", "")))
            anchor_type = project_data.get("anchor_type", "left")
            if anchor_type == "flanking":
                self.flanking_anchor_radio.setChecked(True)
                self.anchor_stack.setCurrentWidget(self.flanking_anchor_widget)
                self.left_anchor_input.setText(project_data["anchors"][0] if len(project_data["anchors"]) > 0 else "")
                self.right_anchor_input.setText(project_data["anchors"][1] if len(project_data["anchors"]) > 1 else "")
                self.range_input.setText("Flanked sequence")
                self.range_input.setEnabled(False)
                self.range_input.setStyleSheet("color: grey;")
            else:
                self.left_anchor_radio.setChecked(anchor_type == "left")
                self.right_anchor_radio.setChecked(anchor_type == "right")
                self.anchor_stack.setCurrentWidget(self.single_anchor_widget)
                self.single_anchor_input.setText(project_data["anchors"][0] if project_data["anchors"] else "")
                self.range_input.setText(project_data.get("range", ""))
                self.range_input.setEnabled(True)
                self.range_input.setStyleSheet("")
            self.reference_sequence_input.setText(project_data.get("reference_sequence", ""))
            self.update_anchor_input()  # Ensure UI reflects loaded state
            barcodes = [seq["barcode"] for seq in project_data["sequences"]]
            primer = project_data["sequences"][0]["primer"] if project_data["sequences"] else ""
            self.sequence_extractor.set_sequences(barcodes, primer, project_data["name"])
            self.editing_analyzer.set_project(project_data)
            self.project_loaded.emit(project_data)
            self.anchors_updated.emit(project_data["anchors"])
            # Update recent projects
            self.update_recent_projects(project_file)
            QMessageBox.information(self, "Success", "Project loaded successfully")
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to load project: {str(e)}")

    def update_anchor_input(self):
        if self.flanking_anchor_radio.isChecked():
            self.anchor_stack.setCurrentWidget(self.flanking_anchor_widget)
            self.range_input.setText("Flanked sequence")
            self.range_input.setEnabled(False)
            self.range_input.setStyleSheet("color: grey;")
        else:
            self.anchor_stack.setCurrentWidget(self.single_anchor_widget)
            self.range_input.setEnabled(True)
            self.range_input.setPlaceholderText("Enter the range of nucleotides for analysis in bp (For example: 20)")
            self.range_input.setStyleSheet("")
            if self.right_anchor_radio.isChecked():
                self.single_anchor_layout.removeWidget(self.single_anchor_input)
                self.single_anchor_layout.removeWidget(self.single_range_label)
                self.single_anchor_layout.addWidget(self.single_range_label)
                self.single_anchor_layout.addWidget(self.single_anchor_input)
                self.single_anchor_input.setPlaceholderText("Right anchor sequence")
                self.single_anchor_input.setAlignment(Qt.AlignmentFlag.AlignLeft)
            else:  # Left Anchor
                self.single_anchor_layout.removeWidget(self.single_anchor_input)
                self.single_anchor_layout.removeWidget(self.single_range_label)
                self.single_anchor_layout.addWidget(self.single_anchor_input)
                self.single_anchor_layout.addWidget(self.single_range_label)
                self.single_anchor_input.setPlaceholderText("Left anchor sequence")
                self.single_anchor_input.setAlignment(Qt.AlignmentFlag.AlignRight)

    def update_sequence_table_headers(self):
        old_type = "Group" if self.replicates_radio.isChecked() else "Replicate"
        new_type = "Replicate" if self.replicates_radio.isChecked() else "Group"
        self.sequence_table.setHorizontalHeaderLabels([new_type, "Barcode", "Primer Sequence"])
        for row in range(self.sequence_table.rowCount()):
            item = self.sequence_table.item(row, 0)
            if item:
                current_text = item.text()
                match = re.match(rf"^{old_type}\s+(\d+)$", current_text)
                if match:
                    number = match.group(1)
                    new_text = f"{new_type} {number}"
                    self.sequence_table.setItem(row, 0, QTableWidgetItem(new_text))

    def set_language(self, language):
        self.current_language = language
        translations = {
            "en": {
                "project_name": "Project Name:",
                "project_path": "Project Path:",
                "browse": "Browse...",
                "new_project": "New Project",
                "save_project": "Save Project",
                "load_project": "Load Project",
                "add_sequence": "Add Barcode Sequence",
                "remove_sequence": "Remove Barcode Sequence",
                "anchor_sequence": "Anchor Sequence",
                "left_anchor": "Left Anchor",
                "right_anchor": "Right Anchor",
                "flanking_anchors": "Flanking Anchors",
                "desired_analyzing_range": "Desired Analyzing Range:",
                "replicate": "Replicate",
                "barcode": "Barcode",
                "primer_sequence": "Primer Sequence",
                "step1": "Step 1: Enter 5'-barcode and primers for sequence extraction",
                "step2": "Step 2: Choose your anchor sequence and analyzing range",
                "barcodes_replicates": "Barcodes for different experimental replicates",
                "barcodes_groups": "Barcodes for different experimental groups",
                "reference_sequence": "Reference Sequence:",
                "recent_projects": "Recent Projects:"  # Added for recent projects
            },
            "zh_TW": {
                "project_name": "專案名稱：",
                "project_path": "專案路徑：",
                "browse": "瀏覽...",
                "new_project": "新建專案",
                "save_project": "儲存專案",
                "load_project": "載入專案",
                "add_sequence": "新增條碼序列",
                "remove_sequence": "移除條碼序列",
                "anchor_sequence": "錨定序列",
                "left_anchor": "左錨定",
                "right_anchor": "右錨定",
                "flanking_anchors": "雙側錨定",
                "desired_analyzing_range": "期望分析範圍：",
                "replicate": "重複",
                "barcode": "條碼",
                "primer_sequence": "引物序列",
                "step1": "步驟 1：輸入 5'-條碼和引物以進行序列提取",
                "step2": "步驟 2：選擇您的錨定序列和分析範圍",
                "barcodes_replicates": "用條碼區分實驗重複",
                "barcodes_groups": "用條碼區分不同實驗組別",
                "reference_sequence": "參考序列：",
                "recent_projects": "最近專案："  # Added for recent projects
            },
            "zh_CN": {
                "project_name": "项目名称：",
                "project_path": "项目路径：",
                "browse": "浏览...",
                "new_project": "新建项目",
                "save_project": "保存项目",
                "load_project": "加载项目",
                "add_sequence": "添加条码序列",
                "remove_sequence": "移除条码序列",
                "anchor_sequence": "锚定序列",
                "left_anchor": "左锚定",
                "right_anchor": "右锚定",
                "flanking_anchors": "双侧锚定",
                "desired_analyzing_range": "期望分析范围：",
                "replicate": "重复",
                "barcode": "条码",
                "primer_sequence": "引物序列",
                "step1": "步骤 1：输入 5'-条码和引物以进行序列提取",
                "step2": "步骤 2：选择您的锚定序列和分析范围",
                "barcodes_replicates": "用条码区分实验重复",
                "barcodes_groups": "用条码区分不同实验组别",
                "reference_sequence": "参考序列：",
                "recent_projects": "最近项目："  # Added for recent projects
            },
            "ja": {
                "project_name": "プロジェクト名：",
                "project_path": "プロジェクトパス：",
                "browse": "参照...",
                "new_project": "新規プロジェクト",
                "save_project": "プロジェクトを保存",
                "load_project": "プロジェクトを読み込む",
                "add_sequence": "バーコードシーケンスを追加",
                "remove_sequence": "バーコードシーケンスを削除",
                "anchor_sequence": "アンカーシーケンス",
                "left_anchor": "左アンカー",
                "right_anchor": "右アンカー",
                "flanking_anchors": "両側アンカー",
                "desired_analyzing_range": "希望する解析範囲：",
                "replicate": "レプリケート",
                "barcode": "バーコード",
                "primer_sequence": "プライマーシーケンス",
                "step1": "ステップ 1：シーケンス抽出のために 5'-バーコードとプライマーを入力",
                "step2": "ステップ 2：アンカーシーケンスと解析範囲を選択",
                "barcodes_replicates": "異なる実験レプリケート用のバーコード",
                "barcodes_groups": "異なる実験グループ用のバーコード",
                "reference_sequence": "参照シーケンス：",
                "recent_projects": "最近のプロジェクト："  # Added for recent projects
            },
            "ko": {
                "project_name": "프로젝트 이름:",
                "project_path": "프로젝트 경로:",
                "browse": "탐색...",
                "new_project": "새 프로젝트",
                "save_project": "프로젝트 저장",
                "load_project": "프로젝트 불러오기",
                "add_sequence": "바코드 시퀀스 추가",
                "remove_sequence": "바코드 시퀀스 제거",
                "anchor_sequence": "앵커 시퀀스",
                "left_anchor": "왼쪽 앵커",
                "right_anchor": "오른쪽 앵커",
                "flanking_anchors": "양측 앵커",
                "desired_analyzing_range": "원하는 분석 범위:",
                "replicate": "반복",
                "barcode": "바코드",
                "primer_sequence": "프라이머 시퀀스",
                "step1": "1단계: 시퀀스 추출을 위해 5'-바코드와 프라이머 입력",
                "step2": "2단계: 앵커 시퀀스와 분석 범위 선택",
                "barcodes_replicates": "다른 실험 반복용 바코드",
                "barcodes_groups": "다른 실험 그룹용 바코드",
                "reference_sequence": "참조 시퀀스:",
                "recent_projects": "최근 프로젝트:"  # Added for recent projects
            },
            "fr": {
                "project_name": "Nom du projet :",
                "project_path": "Chemin du projet :",
                "browse": "Parcourir...",
                "new_project": "Nouveau projet",
                "save_project": "Enregistrer le projet",
                "load_project": "Charger le projet",
                "add_sequence": "Ajouter une séquence de code-barres",
                "remove_sequence": "Supprimer une séquence de code-barres",
                "anchor_sequence": "Séquence d'ancrage",
                "left_anchor": "Ancrage à gauche",
                "right_anchor": "Ancrage à droite",
                "flanking_anchors": "Ancrages flanquants",
                "desired_analyzing_range": "Plage d'analyse souhaitée :",
                "replicate": "Réplicat",
                "barcode": "Code-barres",
                "primer_sequence": "Séquence de l'amorce",
                "step1": "Étape 1 : Entrez le code-barres 5' et les amorces pour l'extraction de séquences",
                "step2": "Étape 2 : Choisissez votre séquence d'ancrage et la plage d'analyse",
                "barcodes_replicates": "Codes-barres pour différents réplicats expérimentaux",
                "barcodes_groups": "Codes-barres pour différents groupes expérimentaux",
                "reference_sequence": "Séquence de référence :",
                "recent_projects": "Projets récents :"  # Added for recent projects
            },
            "de": {
                "project_name": "Projektname:",
                "project_path": "Projektpfad:",
                "browse": "Durchsuchen...",
                "new_project": "Neues Projekt",
                "save_project": "Projekt speichern",
                "load_project": "Projekt laden",
                "add_sequence": "Barcode-Sequenz hinzufügen",
                "remove_sequence": "Barcode-Sequenz entfernen",
                "anchor_sequence": "Ankersequenz",
                "left_anchor": "Linker Anker",
                "right_anchor": "Rechter Anker",
                "flanking_anchors": "Flankierende Anker",
                "desired_analyzing_range": "Gewünschter Analysebereich:",
                "replicate": "Replikat",
                "barcode": "Barcode",
                "primer_sequence": "Primer-Sequenz",
                "step1": "Schritt 1: Geben Sie 5'-Barcode und Primer für die Sequenzextraktion ein",
                "step2": "Schritt 2: Wählen Sie Ihre Ankersequenz und den Analysebereich",
                "barcodes_replicates": "Barcodes für verschiedene experimentelle Replikate",
                "barcodes_groups": "Barcodes für verschiedene experimentelle Gruppen",
                "reference_sequence": "Referenzsequenz:",
                "recent_projects": "Letzte Projekte:"  # Added for recent projects
            },
            "es": {
                "project_name": "Nombre del proyecto:",
                "project_path": "Ruta del proyecto:",
                "browse": "Explorar...",
                "new_project": "Nuevo proyecto",
                "save_project": "Guardar proyecto",
                "load_project": "Cargar proyecto",
                "add_sequence": "Agregar secuencia de código de barras",
                "remove_sequence": "Eliminar secuencia de código de barras",
                "anchor_sequence": "Secuencia de anclaje",
                "left_anchor": "Anclaje izquierdo",
                "right_anchor": "Anclaje derecho",
                "flanking_anchors": "Anclajes flanqueantes",
                "desired_analyzing_range": "Rango de análisis deseado:",
                "replicate": "Réplica",
                "barcode": "Código de barras",
                "primer_sequence": "Secuencia de cebador",
                "step1": "Paso 1: Introduzca el código de barras 5' y los cebadores para la extracción de secuencias",
                "step2": "Paso 2: Elija su secuencia de anclaje y el rango de análisis",
                "barcodes_replicates": "Códigos de barras para diferentes réplicas experimentales",
                "barcodes_groups": "Códigos de barras para diferentes grupos experimentales",
                "reference_sequence": "Secuencia de referencia:",
                "recent_projects": "Proyectos recientes:"  # Added for recent projects
            },
            "pt": {
                "project_name": "Nome do projeto:",
                "project_path": "Caminho do projeto:",
                "browse": "Navegar...",
                "new_project": "Novo projeto",
                "save_project": "Salvar projeto",
                "load_project": "Carregar projeto",
                "add_sequence": "Adicionar sequência de código de barras",
                "remove_sequence": "Remover sequência de código de barras",
                "anchor_sequence": "Sequência de ancoragem",
                "left_anchor": "Ancoragem à esquerda",
                "right_anchor": "Ancoragem à direita",
                "flanking_anchors": "Ancoragens flanqueantes",
                "desired_analyzing_range": "Faixa de análise desejada:",
                "replicate": "Réplica",
                "barcode": "Código de barras",
                "primer_sequence": "Sequência de primer",
                "step1": "Passo 1: Insira o código de barras 5' e os primers para extração de sequências",
                "step2": "Passo 2: Escolha sua sequência de ancoragem e faixa de análise",
                "barcodes_replicates": "Códigos de barras para diferentes réplicas experimentais",
                "barcodes_groups": "Códigos de barras para diferentes grupos experimentais",
                "reference_sequence": "Sequência de referência:",
                "recent_projects": "Projetos recentes:"  # Added for recent projects
            }
        }
        t = translations.get(self.current_language, translations["en"])
        self.name_label.setText(t["project_name"])
        self.path_label.setText(t["project_path"])
        self.browse_btn.setText(t["browse"])
        self.new_btn.setText(t["new_project"])
        self.save_btn.setText(t["save_project"])
        self.load_btn.setText(t["load_project"])
        self.add_sequence_btn.setText(t["add_sequence"])
        self.remove_sequence_btn.setText(t["remove_sequence"])
        self.left_anchor_radio.setText(t["left_anchor"])
        self.right_anchor_radio.setText(t["right_anchor"])
        self.flanking_anchor_radio.setText(t["flanking_anchors"])
        self.range_label.setText(t["desired_analyzing_range"])
        self.step1_label.setText(t["step1"])
        self.step2_label.setText(t["step2"])
        self.replicates_radio.setText(t["barcodes_replicates"])
        self.groups_radio.setText(t["barcodes_groups"])
        self.reference_sequence_label.setText(t["reference_sequence"])
        self.recent_label.setText(t["recent_projects"])  # Added for recent projects
        self.update_sequence_table_headers()

    def setup_ui(self):
        layout = QVBoxLayout()

        # Project information area
        form_layout = QVBoxLayout()
        name_layout = QHBoxLayout()
        self.name_label = QLabel("Project Name:")
        self.project_name = QLineEdit()
        name_layout.addWidget(self.name_label)
        name_layout.addWidget(self.project_name)

        path_layout = QHBoxLayout()
        self.path_label = QLabel("Project Path:")
        self.project_path = QLineEdit()
        self.project_path.setReadOnly(True)
        self.browse_btn = QPushButton("Browse...")
        self.browse_btn.clicked.connect(self.browse_project_path)
        path_layout.addWidget(self.path_label)
        path_layout.addWidget(self.project_path)
        path_layout.addWidget(self.browse_btn)

        form_layout.addLayout(name_layout)
        form_layout.addLayout(path_layout)

        # Project operation buttons
        project_btn_layout = QHBoxLayout()
        self.new_btn = QPushButton("New Project")
        self.save_btn = QPushButton("Save Project")
        self.save_btn.setStyleSheet("QPushButton { background-color: #28A745; color: white; } QPushButton:hover { color: #28A745; }")
        self.load_btn = QPushButton("Load Project")
        self.load_btn.setStyleSheet("QPushButton { background-color: #007BFF; color: white; } QPushButton:hover { color: #007BFF; }")
        project_btn_layout.addWidget(self.new_btn)
        project_btn_layout.addWidget(self.save_btn)
        project_btn_layout.addWidget(self.load_btn)

        # Barcode sequence table
        self.step1_label = QLabel("Step 1: Enter 5'-barcode and primers for sequence extraction")
        self.step1_label.setStyleSheet("color: #007BFF;")
        barcode_type_layout = QHBoxLayout()
        self.replicates_radio = QRadioButton("Barcodes for different experimental replicates")
        self.groups_radio = QRadioButton("Barcodes for different experimental groups")
        self.replicates_radio.setChecked(True)
        barcode_group = QButtonGroup(self)
        barcode_group.addButton(self.replicates_radio)
        barcode_group.addButton(self.groups_radio)
        barcode_group.setExclusive(True)
        barcode_type_layout.addWidget(self.replicates_radio)
        barcode_type_layout.addWidget(self.groups_radio)

        self.sequence_table = QTableWidget()
        self.sequence_table.setColumnCount(3)
        self.sequence_table.setHorizontalHeaderLabels(["Replicate", "Barcode", "Primer Sequence"])
        self.sequence_table.setColumnWidth(0, 100)
        self.sequence_table.setColumnWidth(1, 100)
        self.sequence_table.setColumnWidth(2, 400)
        self.sequence_table.setFixedHeight(200)

        sequence_btn_layout = QHBoxLayout()
        self.add_sequence_btn = QPushButton("Add Barcode Sequence")
        self.remove_sequence_btn = QPushButton("Remove Barcode Sequence")
        sequence_btn_layout.addWidget(self.add_sequence_btn)
        sequence_btn_layout.addWidget(self.remove_sequence_btn)

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.Shape.HLine)
        separator.setFrameShadow(QFrame.Shadow.Sunken)

        # Anchor sequence area
        anchor_layout = QVBoxLayout()
        self.step2_label = QLabel("Step 2: Choose your anchor sequence and analyzing range")
        self.step2_label.setStyleSheet("color: #007BFF;")
        anchor_type_layout = QHBoxLayout()
        self.left_anchor_radio = QRadioButton("Left Anchor")
        self.right_anchor_radio = QRadioButton("Right Anchor")
        self.flanking_anchor_radio = QRadioButton("Flanking Anchors")
        self.left_anchor_radio.setChecked(True)
        anchor_group = QButtonGroup(self)
        anchor_group.addButton(self.left_anchor_radio)
        anchor_group.addButton(self.right_anchor_radio)
        anchor_group.addButton(self.flanking_anchor_radio)
        anchor_group.setExclusive(True)
        anchor_type_layout.addWidget(self.left_anchor_radio)
        anchor_type_layout.addWidget(self.right_anchor_radio)
        anchor_type_layout.addWidget(self.flanking_anchor_radio)

        range_layout = QHBoxLayout()
        self.range_label = QLabel("Desired Analyzing Range:")
        self.range_input = QLineEdit()
        self.range_input.setPlaceholderText("Enter the range of nucleotides for analysis in bp (For example: 20)")
        range_layout.addWidget(self.range_label)
        range_layout.addWidget(self.range_input)

        self.anchor_stack = QStackedWidget()

        # Single anchor widget (Left or Right)
        self.single_anchor_widget = QWidget()
        self.single_anchor_layout = QHBoxLayout(self.single_anchor_widget)
        self.single_anchor_input = QLineEdit()
        self.single_anchor_input.setPlaceholderText("Left anchor sequence")
        self.single_anchor_input.setFixedHeight(30)
        self.single_anchor_input.setAlignment(Qt.AlignmentFlag.AlignRight)
        self.single_range_label = QLabel("|<------- Analyzing Range -------->|")
        self.single_range_label.setStyleSheet("color: #007BFF;")
        self.single_range_label.setFixedHeight(30)
        self.single_range_label.setFixedWidth(200)
        self.single_anchor_layout.addWidget(self.single_anchor_input)
        self.single_anchor_layout.addWidget(self.single_range_label)

        # Flanking anchor widget
        self.flanking_anchor_widget = QWidget()
        flanking_anchor_layout = QHBoxLayout(self.flanking_anchor_widget)
        self.left_anchor_input = QLineEdit()
        self.left_anchor_input.setPlaceholderText("Left anchor sequence")
        self.left_anchor_input.setFixedHeight(30)
        self.left_anchor_input.setAlignment(Qt.AlignmentFlag.AlignRight)
        self.flanking_range_label = QLabel("|<------- Analyzing Range -------->|")
        self.flanking_range_label.setStyleSheet("color: #007BFF;")
        self.flanking_range_label.setFixedHeight(30)
        self.flanking_range_label.setFixedWidth(200)
        self.right_anchor_input = QLineEdit()
        self.right_anchor_input.setPlaceholderText("Right anchor sequence")
        self.right_anchor_input.setFixedHeight(30)
        self.right_anchor_input.setAlignment(Qt.AlignmentFlag.AlignLeft)
        flanking_anchor_layout.addWidget(self.left_anchor_input)
        flanking_anchor_layout.addWidget(self.flanking_range_label)
        flanking_anchor_layout.addWidget(self.right_anchor_input)

        self.anchor_stack.addWidget(self.single_anchor_widget)
        self.anchor_stack.addWidget(self.flanking_anchor_widget)

        anchor_layout.addWidget(self.step2_label)
        anchor_layout.addLayout(anchor_type_layout)
        anchor_layout.addLayout(range_layout)
        anchor_layout.addWidget(self.anchor_stack)

        # Reference sequence input
        reference_layout = QHBoxLayout()
        self.reference_sequence_label = QLabel("Reference Sequence:")
        self.reference_sequence_input = QLineEdit()
        self.reference_sequence_input.setPlaceholderText("Enter reference sequence (optional)")
        reference_layout.addWidget(self.reference_sequence_label)
        reference_layout.addWidget(self.reference_sequence_input)

        # Recent projects list
        recent_layout = QVBoxLayout()
        self.recent_label = QLabel("Recent Projects:")
        recent_layout.addWidget(self.recent_label)
        self.recent_projects = QListWidget()
        recent_layout.addWidget(self.recent_projects)
        self.recent_projects.setFixedHeight(100)

        layout.addLayout(form_layout)
        layout.addLayout(project_btn_layout)
        layout.addWidget(self.step1_label)
        layout.addLayout(barcode_type_layout)
        layout.addWidget(self.sequence_table)
        layout.addLayout(sequence_btn_layout)
        layout.addWidget(separator)
        layout.addLayout(anchor_layout)
        layout.addLayout(reference_layout)
        layout.addLayout(recent_layout)

        self.setLayout(layout)

        # Connect signals
        self.new_btn.clicked.connect(self.new_project)
        self.save_btn.clicked.connect(self.save_project)
        self.load_btn.clicked.connect(self.load_project)
        self.add_sequence_btn.clicked.connect(self.add_sequence)
        self.remove_sequence_btn.clicked.connect(self.remove_sequence)
        self.left_anchor_radio.toggled.connect(self.update_anchor_input)
        self.right_anchor_radio.toggled.connect(self.update_anchor_input)
        self.flanking_anchor_radio.toggled.connect(self.update_anchor_input)
        self.replicates_radio.toggled.connect(self.update_sequence_table_headers)
        self.groups_radio.toggled.connect(self.update_sequence_table_headers)
        self.recent_projects.itemDoubleClicked.connect(self.load_recent_project)  # Connect recent projects

        # Load recent projects
        self.load_recent_projects()

        # Ensure initial state is correct
        self.update_anchor_input()

    def load_recent_project(self, item):
        """Load a project when double-clicking it in the recent projects list"""
        project_file = item.text()
        if os.path.exists(project_file):
            self.load_project_file(project_file)
        else:
            QMessageBox.warning(self, "Warning", "Project file not found")
            self.load_recent_projects()

    def update_recent_projects(self, project_file):
        """Update the list of recent projects"""
        recent_projects = self.load_recent_projects_list()
        if project_file not in recent_projects:
            recent_projects.insert(0, project_file)
            if len(recent_projects) > 5:
                recent_projects.pop()
        
        with open('recent_projects.json', 'w') as f:
            json.dump(recent_projects, f)
        
        self.load_recent_projects()

    def load_recent_projects(self):
        """Load and display recent projects"""
        self.recent_projects.clear()
        recent_projects = self.load_recent_projects_list()
        self.recent_projects.addItems(recent_projects)

    def load_recent_projects_list(self):
        """Load recent projects from file"""
        if os.path.exists('recent_projects.json'):
            try:
                with open('recent_projects.json', 'r') as f:
                    return json.load(f)
            except:
                pass
        return []