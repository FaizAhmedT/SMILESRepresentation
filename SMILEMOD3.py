import sys
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QLabel, QLineEdit, QPushButton,
    QVBoxLayout, QWidget, QTextEdit, QDialog, QDialogButtonBox,
    QAction, QMenuBar, QFileDialog, QStatusBar, QHBoxLayout, QComboBox, QCompleter
)
from PyQt5.QtGui import QPixmap
from PIL import Image


class MoleculeDrawer(QMainWindow):

    compound_mapping = {
    'Aspirin': 'CHEMBL25',
    'Paracetamol': 'CHEMBL112',
    'Caffeine': 'CHEMBL113',
    'Ibuprofen': 'CHEMBL585',
    'Loratadine': 'CHEMBL628',
    'Simvastatin': 'CHEMBL1487',
    'Metformin': 'CHEMBL1431',
    'Amoxicillin': 'CHEMBL62',
    'Omeprazole': 'CHEMBL1543',
    'Acetaminophen': 'CHEMBL112',
    'Atorvastatin': 'CHEMBL1489',
    'Diazepam': 'CHEMBL607',
    'Citalopram': 'CHEMBL599',
    'Warfarin': 'CHEMBL95',
    'Morphine': 'CHEMBL70',
    'Sildenafil': 'CHEMBL192',
    'Metoprolol': 'CHEMBL15',
    'Lisinopril': 'CHEMBL554',
    'Hydrochlorothiazide': 'CHEMBL155',
    'Montelukast': 'CHEMBL828',
    'Clopidogrel': 'CHEMBL1771',
    'Levothyroxine': 'CHEMBL398982',
    'Pregabalin': 'CHEMBL1059',
    'Quetiapine': 'CHEMBL715',
    'Carbamazepine': 'CHEMBL108',
    'Tramadol': 'CHEMBL122',
    'Fluoxetine': 'CHEMBL41',
    'Ranitidine': 'CHEMBL1900',
    'Lansoprazole': 'CHEMBL109',
    'Ceftriaxone': 'CHEMBL1264',
    'Digoxin': 'CHEMBL1531',
    'Tamsulosin': 'CHEMBL1231',
    'Furosemide': 'CHEMBL70',
    'Venlafaxine': 'CHEMBL551',
    'Allopurinol': 'CHEMBL59',
    'Doxycycline': 'CHEMBL1433',
    'Gabapentin': 'CHEMBL37',
    'Methotrexate': 'CHEMBL688',
    'Enalapril': 'CHEMBL165',
    'Losartan': 'CHEMBL191',
    'Clonazepam': 'CHEMBL452',
    'Cefixime': 'CHEMBL1755',
    'Esomeprazole': 'CHEMBL53',
    'Amitriptyline': 'CHEMBL799',
    'Hydroxyzine': 'CHEMBL707',
    'Bupropion': 'CHEMBL203',
    'Valsartan': 'CHEMBL1069',
    # Add more compounds as needed
}

    compound_smiles = {
    'Aspirin': 'O=C(C)Oc1ccccc1C(=O)O',
    'Paracetamol': 'CC(=O)Nc1ccc(O)cc1',
    'Caffeine': 'CN1C=NC2=C1C(=O)N(C(=O)N2C)C',
    'Ibuprofen': 'CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O',
    'Loratadine': 'CC(C)(C)OC(=O)N(c1ccc(Cl)cc1)c2ccc(Cl)cc2',
    'Simvastatin': 'CC(CC(=O)O[C@H]1CC(C=C2[C@@H](C)C=C(C)C2=O)OC1)C(C)C',
    'Metformin': 'CN(C)C(=N)N',
    'Amoxicillin': 'CC1([C@@H](N(C(=O)C1=O)C)O)C',
    'Omeprazole': 'CC1=CN=C(C(=C1)[S](=O)(=O)N2CCCC2)OC',
    'Acetaminophen': 'CC(=O)Nc1ccc(O)cc1',
    'Atorvastatin': 'CC(C)(C)OC(=O)C(C)OC(=O)C(C)C(C)C(C)C(C)C(CC(C)O)OC=CC(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C',
    'Diazepam': 'CN=C(C)C(=O)N1CC(=O)N(C(c2ccccc2)(c2ccccc2)C2=CC=CC=C2)C1=O',
    'Citalopram': 'CN(C)CCC(=O)c1ccc(NC(=O)CCC2=CC=CC=C2)cc1',
    'Warfarin': 'CC(=O)CC(c1ccccc1)C1=C(O)C(=O)C(c2ccccc2)C1=O',
    'Morphine': 'CC(=O)N(CC1=CC2=C(C=C1)OCO2)C(C)C',
    'Sildenafil': 'CN(C)C(=O)c1c(C)nc2ccccc2n1C',
    'Metoprolol': 'CC(CC(=O)OCC)NCC(C)CO',
    'Lisinopril': 'CC(C)C(C(=O)N1CCCC1C(=O)N(C)C)N(C)C',
    'Hydrochlorothiazide': 'CN(C)S(=O)(=O)c1ccc(C(Cl)Cl)cc1',
    'Montelukast': 'CC1=CC(=O)C2=C(C1=O)N=C(N2)CC1=CN=C(S1)C1=CC=C(N=C1)C(F)(F)F',
    'Clopidogrel': 'CC(=O)C1=C(C=CS1(=O)=O)C1=CC=CC=C1',
    'Levothyroxine': 'CC(C)NCC(C(=O)NC1=CC=C(O)C=C1)C1=CC=CC=C1',
    'Pregabalin': 'CC(C)(C(=O)O)NCCN1C=CC(=C1)C2=CC=CC=N2',
    'Quetiapine': 'O=C(CN1CCN(CC1)C2=CC=C(C=C2)C2=CC=CC=C2)NCC',
    'Carbamazepine': 'CC(=O)Nc1ccc(-c2ccccc2)c(-c2ccccc2)n1',
    'Tramadol': 'CN(C)CC(C1=CC=CC=C1)O',
    'Fluoxetine': 'CN(C)CCOCC1=CC=C(C=C1)C1=CC=CC=C1',
    'Ranitidine': 'CC(NC(NCC1=CC=CC=C1)=O)NCC1=CC=CC=C1',
    'Lansoprazole': 'CC(C)(C)OC(=O)Nc1nnc(-c2ccc(Cl)cc2)n1',
    'Ceftriaxone': 'CC1([C@@H]2NC(=O)[C@H](N)C(=O)[C@H](S)C(=O)[C@@H](N)C(=O)[C@H](C)NC(=O)[C@@H](N)C(=O)O)SC2(C)C1=O',
    'Digoxin': 'CC12CCC3C(C1CC=C2C(=O)C1CC2(C)CCC3(O)CC=C4C(C(=O)C(C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)C(C(=O)O1)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O)O',
    'Tamsulosin': 'CC(C)OC1=CC(=O)N(C)C=C(C2=CC=C(O)C=C2)C1=O',
    'Furosemide': 'CC(CC(C)C)C1=CC=C(C=C1)S(=O)(=O)NC1=NC(=O)NC(=O)C1',
    'Venlafaxine': 'CC(C)(C)OC(=O)C1=CC(=C(C=C1)Cl)NCC(C)C',
    'Allopurinol': 'CC1=CN=CN=C1NCC(O)CO',
    'Doxycycline': 'CC(C)C1=C(C(=C(N(C1=O)C)C2=CC=C(C=C2)O)C(=O)O)O',
    'Gabapentin': 'CC(C(=O)O)NCC(=O)NCC(=O)NCC1=CC=CC=C1',
    'Methotrexate': 'CC(C)(C(=O)O)C1=CC(=O)C2=C1C=CN2',
    'Enalapril': 'CC(C)NCC(O)=O',
    'Losartan': 'CC(C)COC(=O)CN(C)C(=O)C1=CC=C(C=C1)C1=CC=CC=C1',
    'Clonazepam': 'CC1=CC(=O)N(C2=CC=CC=C2)C2=CC=CC=C21',
    'Cefixime': 'CC(C)(C)OC(=O)N1C(=C(C2=C(C(=O)[O-])N(C2=O)C2=CC=C(C=C2)O)S1(=O)=O)C(=O)O',
    'Esomeprazole': 'CC1=CN=C(C(=C1)[S](=O)(=O)N2CCC(C2)C(=O)O)OC',
    'Amitriptyline': 'CC(C)(C)NCCC=C1C2=CC=CC=C2SC3=CC=CC=C31',
    'Hydroxyzine': 'CC(C)(CC1=CC2=C(C=C1)OCO2)O',
    'Bupropion': 'CC(CC(=O)NC1=CC=C(C=C1)Cl)C1=CC=CC=C1',
    'Valsartan': 'CC(C)OC(=O)C1=CC(=CC=C1)NC(=O)N(C)C1=CC=C(C=C1)C1=CC=C(C=C1)C1=CC=CC=C1',
}

    def __init__(self):
        super().__init__()

        self.setWindowTitle("Molecule Drawer")
        self.setGeometry(100, 100, 800, 600)

        self.create_menu()
        self.create_main_layout()
        self.create_status_bar()

        self.compound_combobox.currentIndexChanged.connect(self.update_fields_from_combobox)

        # Variable to store the currently drawn molecule
        self.current_mol = None

    def create_menu(self):
        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu('File')

        open_action = QAction('Open', self)
        open_action.triggered.connect(self.open_molecule_file)
        file_menu.addAction(open_action)

        exit_action = QAction('Exit', self)
        exit_action.triggered.connect(self.close)
        file_menu.addAction(exit_action)

    def create_main_layout(self):
        central_widget = QWidget(self)
        self.setCentralWidget(central_widget)

        layout = QVBoxLayout(central_widget)

        # Compound Selection
        self.compound_label = QLabel("Select Compound:")
        layout.addWidget(self.compound_label)

        self.compound_combobox = QComboBox()
        self.compound_combobox.addItems(self.compound_mapping.keys())
        self.compound_combobox.setEditable(True)
        layout.addWidget(self.compound_combobox)

        self.smiles_label = QLabel("Enter SMILES:")
        layout.addWidget(self.smiles_label)

        self.smiles_entry = QLineEdit()
        layout.addWidget(self.smiles_entry)

        self.chembl_id_label = QLabel("Enter CHEMBL ID:")
        layout.addWidget(self.chembl_id_label)

        self.chembl_id_entry = QLineEdit()
        layout.addWidget(self.chembl_id_entry)

        self.draw_button = QPushButton("Draw Molecule", self)
        self.draw_button.clicked.connect(self.draw_molecule)
        layout.addWidget(self.draw_button)

        self.info_text = QTextEdit()
        layout.addWidget(self.info_text)

        # Additional Functionality Buttons
        button_layout = QHBoxLayout()

        self.clear_button = QPushButton("Clear", self)
        self.clear_button.clicked.connect(self.clear_fields)
        button_layout.addWidget(self.clear_button)

        self.show_descriptors_button = QPushButton("Show Descriptors", self)
        self.show_descriptors_button.clicked.connect(self.show_additional_descriptors)
        button_layout.addWidget(self.show_descriptors_button)

        self.export_info_button = QPushButton("Export Info", self)
        self.export_info_button.clicked.connect(self.export_molecule_info)
        button_layout.addWidget(self.export_info_button)

        self.open_image_button = QPushButton("Open Image", self)
        self.open_image_button.clicked.connect(self.open_image)
        button_layout.addWidget(self.open_image_button)

        self.calculate_molecular_weight_button = QPushButton("Calculate Molecular Weight", self)
        self.calculate_molecular_weight_button.clicked.connect(self.calculate_molecular_weight)
        button_layout.addWidget(self.calculate_molecular_weight_button)

        self.molecular_info_button = QPushButton("Molecular Info", self)
        self.molecular_info_button.clicked.connect(self.show_molecular_info)
        button_layout.addWidget(self.molecular_info_button)

        layout.addLayout(button_layout)
    
    def update_fields_from_combobox(self):
        selected_compound = self.compound_combobox.currentText()
        chembl_id = self.compound_mapping.get(selected_compound, '')
        smi = self.compound_smiles.get(selected_compound, '')

        self.smiles_entry.setText(smi)
        self.chembl_id_entry.setText(chembl_id)

    def create_status_bar(self):
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)

    def open_molecule_file(self):
        file_dialog = QFileDialog()
        file_path, _ = file_dialog.getOpenFileName(self, 'Open Molecule File', '', 'Molecule Files (*.mol *.sdf *.smi)')

        if file_path:
            with open(file_path, 'r') as file:
                mol_data = file.read()
                self.smiles_entry.setText(mol_data)

    def draw_molecule(self):
        # Use the selected compound name to get the corresponding CHEMBLID
        selected_compound = self.compound_combobox.currentText()
        chembl_id = self.compound_mapping.get(selected_compound, '')

        # Use the SMILES associated with the selected compound
        smi = self.compound_smiles.get(selected_compound, '')

        self.current_mol = Chem.MolFromSmiles(smi)
        smi = self.smiles_entry.text()
        chembl_id = self.chembl_id_entry.text()

        if self.current_mol:
            img = Draw.MolToImage(self.current_mol)
            img.save(f'{chembl_id}.png')
            pixmap = QPixmap(f'{chembl_id}.png')

            molecule_dialog = QDialog(self)
            molecule_dialog.setWindowTitle(f'Molecule: {chembl_id}')

            label = QLabel()
            label.setPixmap(pixmap)

            info_text = QTextEdit()
            info_text.setPlainText(f'Molecule: {chembl_id}\nSMILES: {smi}')

            button_box = QDialogButtonBox(QDialogButtonBox.Save | QDialogButtonBox.Close)
            button_box.accepted.connect(molecule_dialog.accept)
            button_box.rejected.connect(molecule_dialog.reject)
            button_box.button(QDialogButtonBox.Save).clicked.connect(self.save_image)

            layout = QVBoxLayout(molecule_dialog)
            layout.addWidget(label)
            layout.addWidget(info_text)
            layout.addWidget(button_box)

            molecule_dialog.exec_()

            # Update status bar
            self.statusBar.showMessage("Molecule drawn successfully", 3000)
        else:
            self.statusBar.showMessage("Invalid SMILES string. Please try again.", 3000)

    def save_image(self):
        chembl_id = self.chembl_id_entry.text()
        img = Image.open(f'{chembl_id}.png')
        img.save(f'{chembl_id}_pyqt.png')
        self.statusBar.showMessage("Image saved successfully", 3000)

    def clear_fields(self):
        self.smiles_entry.clear()
        self.chembl_id_entry.clear()
        self.info_text.clear()

    def show_additional_descriptors(self):
        smi = self.smiles_entry.text()
        chembl_id = self.chembl_id_entry.text()

        mol = Chem.MolFromSmiles(smi)

        if mol:
            descriptor_window = QDialog(self)
            descriptor_window.setWindowTitle("Additional Molecular Descriptors")

            # Create a QVBoxLayout for the descriptor_window
            descriptor_layout = QVBoxLayout(descriptor_window)

            descriptors = self.get_additional_descriptors(mol)
            for key, value in descriptors.items():
                descriptor_label = QLabel(f'{key}: {value}')
                descriptor_layout.addWidget(descriptor_label)

            descriptor_window.exec_()
        else:
            self.statusBar.showMessage("Invalid SMILES string for descriptors. Please try again.", 3000)

    def get_additional_descriptors(self, mol):
        descriptors = {
            'TPSA': Descriptors.TPSA(mol),
            'NumRotatableBonds': Descriptors.NumRotatableBonds(mol),
            'NumHeavyAtoms': Descriptors.HeavyAtomCount(mol),
            'NumRings': Descriptors.RingCount(mol),
            'NumHeteroatoms': Descriptors.NumHeteroatoms(mol),
        }
        return descriptors

    def export_molecule_info(self):
        smi = self.smiles_entry.text()
        chembl_id = self.chembl_id_entry.text()

        mol = Chem.MolFromSmiles(smi)

        if mol:
            file_dialog = QFileDialog()
            file_path, _ = file_dialog.getSaveFileName(self, 'Export Molecule Info', '', 'Text Files (*.txt)')

            if file_path:
                with open(file_path, 'w') as file:
                    info = f'Molecule: {chembl_id}\nSMILES: {smi}\n\n'
                    mol_info = self.get_molecular_info(mol)
                    for key, value in mol_info.items():
                        info += f'{key}: {value}\n'
                    file.write(info)
                self.statusBar.showMessage("Molecule info exported successfully", 3000)
        else:
            self.statusBar.showMessage("Invalid SMILES string for export. Please try again.", 3000)

    def get_molecular_info(self, mol):
        info = {
            'Number of Atoms': mol.GetNumAtoms(),
            'Molecular Weight': Descriptors.MolWt(mol),
            'LogP': Descriptors.MolLogP(mol),
            'HBA (Acceptors)': Descriptors.NumHAcceptors(mol),
            'HBD (Donors)': Descriptors.NumHDonors(mol),
        }
        return info

    def open_image(self):
        image_path, _ = QFileDialog.getOpenFileName(self, 'Open Image File', '', 'Image Files (*.png *.jpg *.jpeg *.bmp)')

        if image_path:
            pixmap = QPixmap(image_path)

            # Display the image in a QDialog
            image_dialog = QDialog(self)
            image_dialog.setWindowTitle('Opened Image')

            label = QLabel()
            label.setPixmap(pixmap)

            layout = QVBoxLayout(image_dialog)
            layout.addWidget(label)

            image_dialog.exec_()
            self.statusBar.showMessage("Image opened successfully", 3000)

    def calculate_molecular_weight(self):
        smi = self.smiles_entry.text()
        chembl_id = self.chembl_id_entry.text()

        mol = Chem.MolFromSmiles(smi)

        if mol:
            molecular_weight = Descriptors.MolWt(mol)
            self.statusBar.showMessage(f"Molecular Weight: {molecular_weight:.4f}", 3000)
        else:
            self.statusBar.showMessage("Invalid SMILES string for calculating molecular weight. Please try again.", 3000)

    def show_molecular_info(self):
        smi = self.smiles_entry.text()
        chembl_id = self.chembl_id_entry.text()

        mol = Chem.MolFromSmiles(smi)

        if mol:
            molecular_info_dialog = QDialog(self)
            molecular_info_dialog.setWindowTitle(f'Molecular Info: {chembl_id}')

            info_text = QTextEdit()
            mol_info = self.get_molecular_info(mol)
            info_text.setPlainText(f'Molecule: {chembl_id}\nSMILES: {smi}\n\n')
            for key, value in mol_info.items():
                info_text.append(f'{key}: {value}')

            layout = QVBoxLayout(molecular_info_dialog)
            layout.addWidget(info_text)

            molecular_info_dialog.exec_()
        else:
            self.statusBar.showMessage("Invalid SMILES string for molecular info. Please try again.", 3000)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MoleculeDrawer()
    window.show()
    sys.exit(app.exec_())
