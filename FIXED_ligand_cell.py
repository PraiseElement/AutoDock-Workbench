# ==================== FIXED LIGAND PREPARATION ====================
# 
# This cell replaces the original ligand preparation section.
# It fixes:
#   1. SDF files now work (auto-converted via RDKit)
#   2. 3D preview displays correctly
#   3. Removed confusing "Inactivate" field
#   4. Better error messages

ligand_upload = widgets.FileUpload(
    accept='.pdb,.mol2,.sdf,.pdbqt',
    multiple=False,
    description='Ligand File',
    layout=Layout(width='auto')
)

# Simplified repairs dropdown with helpful descriptions
ligand_repairs = widgets.Dropdown(
    options=[
        ('Add bonds and hydrogens (recommended)', 'bonds_hydrogens'),
        ('Add hydrogens only', 'hydrogens'),
        ('Check and add missing hydrogens', 'checkhydrogens'),
        ('No repairs', '')
    ],
    value='bonds_hydrogens',
    description='Repairs:',
    style={'description_width': '80px'}
)

ligand_output = widgets.Output()
ligand_3d_view = widgets.Output()
ligand_2d_view = widgets.Output()

prepare_ligand_btn = widgets.Button(
    description='Prepare Ligand',
    button_style='primary',
    icon='check'
)

def prepare_ligand(btn):
    with ligand_output:
        clear_output()
        if not ligand_upload.value:
            print('Please upload a ligand file first!')
            return
        
        # Save uploaded file
        uploaded_file = ligand_upload.value[0]
        filename = uploaded_file['name']
        content = uploaded_file['content']
        ligand_path = state.work_dir / filename
        
        with open(ligand_path, 'wb') as f:
            f.write(content)
        
        state.ligand_file = ligand_path
        print(f'Saved: {ligand_path}')
        
        output_pdbqt = state.work_dir / f'{ligand_path.stem}_prepared.pdbqt'
        repairs = ligand_repairs.value
        
        # For SDF files, use RDKit to convert to PDB first (MolKit can't read SDF)
        if str(ligand_path).lower().endswith('.sdf'):
            print('Converting SDF via RDKit...')
            try:
                rdkit_mol = Chem.SDMolSupplier(str(ligand_path))[0]
                if rdkit_mol is None:
                    print('Error: Could not read SDF file')
                    return
                rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
                # Generate 3D coordinates if needed
                if rdkit_mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(rdkit_mol, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(rdkit_mol)
                # Save as PDB for MolKit
                pdb_path = state.work_dir / f'{ligand_path.stem}_temp.pdb'
                Chem.MolToPDBFile(rdkit_mol, str(pdb_path))
                ligand_path = pdb_path
                print(f'Converted to: {pdb_path.name}')
            except Exception as e:
                print(f'Error converting SDF: {e}')
                return
        
        # Prepare ligand with MolKit
        try:
            mol = Read(str(ligand_path))[0]
            print(f'Loaded: {len(mol.allAtoms)} atoms')
            
            prep = AD4LigandPreparation(
                mol,
                mode='automatic',
                repairs=repairs,
                charges_to_add='gasteiger',
                cleanup='nphs_lps',
                outputfilename=str(output_pdbqt)
            )
            
            state.ligand_pdbqt = output_pdbqt
            print(f'Ligand prepared: {output_pdbqt.name}')
            # Show rotatable bonds if available
            if hasattr(prep, 'ndihe'):
                print(f'   Rotatable bonds: {prep.ndihe}')
            elif hasattr(prep, 'TORSDOF'):
                print(f'   Rotatable bonds: {prep.TORSDOF}')
            
        except Exception as e:
            print(f'Error: {e}')
            import traceback
            traceback.print_exc()
            return
    
    # 3D visualization - Convert PDBQT to simple PDB for py3Dmol
    with ligand_3d_view:
        clear_output()
        if state.ligand_pdbqt and state.ligand_pdbqt.exists():
            try:
                # Read PDBQT and convert to simple PDB (py3Dmol doesn't handle BRANCH/ROOT well)
                pdb_lines = []
                with open(state.ligand_pdbqt, 'r') as f:
                    for line in f:
                        if line.startswith(('HETATM', 'ATOM')):
                            # Convert PDBQT line to PDB (remove last column - atom type)
                            pdb_lines.append(line[:66].rstrip() + '\n')
                pdb_lines.append('END\n')
                pdb_content = ''.join(pdb_lines)
                
                view = py3Dmol.view(width=400, height=300)
                view.addModel(pdb_content, 'pdb')
                view.setStyle({'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.15}})
                view.addStyle({'elem': 'Cl'}, {'stick': {'color': 'green', 'radius': 0.2}})
                view.addStyle({'elem': 'O'}, {'stick': {'color': 'red', 'radius': 0.18}})
                view.setBackgroundColor('0x1a1a2e')
                view.zoomTo()
                view.show()
            except Exception as e:
                print(f'3D view error: {e}')
    
    # 2D visualization with RDKit
    with ligand_2d_view:
        clear_output()
        try:
            orig_path = state.ligand_file
            if str(orig_path).lower().endswith('.mol2'):
                rdkit_mol = Chem.MolFromMol2File(str(orig_path), sanitize=False, removeHs=False)
                if rdkit_mol:
                    try:
                        Chem.SanitizeMol(rdkit_mol)
                    except:
                        pass  # Use unsanitized mol for display
            elif str(orig_path).lower().endswith('.sdf'):
                rdkit_mol = Chem.SDMolSupplier(str(orig_path), sanitize=False, removeHs=False)[0]
                if rdkit_mol:
                    try:
                        Chem.SanitizeMol(rdkit_mol)
                    except:
                        pass
            else:
                # PDB files often have non-standard atom names causing "Element 'X' not found"
                # Use permissive parsing to avoid this error
                rdkit_mol = Chem.MolFromPDBFile(str(orig_path), sanitize=False, removeHs=False)
                if rdkit_mol:
                    try:
                        Chem.SanitizeMol(rdkit_mol)
                    except:
                        pass  # Use unsanitized mol for display
            
            if rdkit_mol:
                img = Draw.MolToImage(rdkit_mol, size=(300, 250))
                display(img)
                print(f'MW: {Descriptors.MolWt(rdkit_mol):.2f}')
                print(f'LogP: {Descriptors.MolLogP(rdkit_mol):.2f}')
                print(f'HBD: {Descriptors.NumHDonors(rdkit_mol)}')
                print(f'HBA: {Descriptors.NumHAcceptors(rdkit_mol)}')
        except Exception as e:
            print(f'2D view not available: {e}')

prepare_ligand_btn.on_click(prepare_ligand)

# Help text
help_text = widgets.HTML('''
<div style="background:#2d2d4a; padding:10px; border-radius:5px; margin:5px 0;">
<b>Supported formats:</b> PDB, MOL2, SDF, PDBQT<br>
<b>Tip:</b> SDF files from PubChem/ChEMBL work great!
</div>
''')

ligand_panel = VBox([
    widgets.HTML('<h3>Upload Ligand</h3>'),
    help_text,
    ligand_upload,
    ligand_repairs,
    prepare_ligand_btn,
    ligand_output,
    HBox([ligand_3d_view, ligand_2d_view])
])

display(ligand_panel)
