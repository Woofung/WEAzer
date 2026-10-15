from PyQt6.QtWidgets import QApplication, QDialog, QVBoxLayout, QTextEdit, QPushButton, QHBoxLayout, QMessageBox
from main_window import MainWindow
import os
import sys

class LicenseDialog(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("WEAzer License Agreement")
        self.setMinimumSize(600, 400)
        self.setup_ui()

    def setup_ui(self):
        layout = QVBoxLayout()

        # License text
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

        # Buttons
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

class EditingAnalyzer:
    def __init__(self):
        self.check_license_agreement()  # Check license before proceeding
        self.main_window = MainWindow()
        self.main_window.show()

    def check_license_agreement(self):
        # Path to a flag file indicating license acceptance
        flag_file = os.path.expanduser("~/.weazer_license_accepted")
        if not os.path.exists(flag_file):
            app = QApplication.instance() or QApplication(sys.argv)  # Ensure QApplication exists
            dialog = LicenseDialog()
            if dialog.exec() == QDialog.DialogCode.Accepted:  # Corrected to use DialogCode.Accepted
                # User agreed, create the flag file
                with open(flag_file, 'w') as f:
                    f.write("License accepted on {}".format(__import__('datetime').datetime.now()))
            else:
                # User disagreed, show message and exit
                QMessageBox.critical(None, "License Not Accepted", 
                                     "You must agree to the license to use WEAzer. The application will now close.")
                sys.exit(1)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    analyzer = EditingAnalyzer()
    sys.exit(app.exec())