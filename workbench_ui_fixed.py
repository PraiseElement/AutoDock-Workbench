#!/usr/bin/env python
"""
AutoDock Workbench - Fixed Interactive Version

This script provides a complete interactive docking workflow with fixes for:
1. SDF file handling via RDKit conversion
2. Better 3D visualization
3. Clearer UI options

Run this in a Jupyter notebook cell after the setup cells.
"""

import os
import sys
import warnings
warnings.filterwarnings('ignore')

import numpy as np
from pathlib import Path
import py3Dmol
from IPython.display import display, HTML, clear_output
import ipywidgets as widgets
from ipywidgets import Layout, VBox, HBox

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors

from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation, AD4LigandPreparation


# ==================== LIGAND PREPARATION (FIXED) ====================

def create_ligand_panel(state):
    """Create fixed ligand preparation panel"""
    
    ligand_upload = widgets.FileUpload(
        accept='.pdb,.mol2,.sdf,.pdbqt',
        multiple=False,
        description='Ligand File',
        layout=Layout(width='auto')
    )
    
    # Simplified options with helpful descriptions
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
            
            # For SDF files, use RDKit to convert to PDB first
            if str(ligand_path).lower().endswith('.sdf'):
                print('Converting SDF via RDKit...')
                try:
                    rdkit_mol = Chem.SDMolSupplier(str(ligand_path))[0]
                    if rdkit_mol is None:
                        print('Error: Could not read SDF file with RDKit')
                        return
                    rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
                    # Generate 3D if needed
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
                print(f'   Rotatable bonds: {prep.ndihe}')
                
            except Exception as e:
                print(f'Error preparing ligand: {e}')
                import traceback
                traceback.print_exc()
                return
        
        # 3D visualization
        with ligand_3d_view:
            clear_output()
            if state.ligand_pdbqt and state.ligand_pdbqt.exists():
                try:
                    with open(state.ligand_pdbqt, 'r') as f:
                        pdbqt_content = f.read()
                    view = py3Dmol.view(width=400, height=300)
                    view.addModel(pdbqt_content, 'pdb')
                    view.setStyle({'stick': {'colorscheme': 'cyanCarbon'}})
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
                            pass
                elif str(orig_path).lower().endswith('.sdf'):
                    rdkit_mol = Chem.SDMolSupplier(str(orig_path), sanitize=False, removeHs=False)[0]
                    if rdkit_mol:
                        try:
                            Chem.SanitizeMol(rdkit_mol)
                        except:
                            pass
                else:
                    # PDB files often have non-standard atom names causing "Element 'X' not found"
                    rdkit_mol = Chem.MolFromPDBFile(str(orig_path), sanitize=False, removeHs=False)
                    if rdkit_mol:
                        try:
                            Chem.SanitizeMol(rdkit_mol)
                        except:
                            pass
                
                if rdkit_mol:
                    img = Draw.MolToImage(rdkit_mol, size=(300, 250))
                    display(img)
                    print(f'MW: {Descriptors.MolWt(rdkit_mol):.2f}')
                    print(f'LogP: {Descriptors.MolLogP(rdkit_mol):.2f}')
                    print(f'HBD: {Descriptors.NumHDonors(rdkit_mol)}')
                    print(f'HBA: {Descriptors.NumHAcceptors(rdkit_mol)}')
            except Exception as e:
                print(f'2D view error: {e}')
    
    prepare_ligand_btn.on_click(prepare_ligand)
    
    # Help text
    help_text = widgets.HTML('''
    <div style="background:#2d2d4a; padding:10px; border-radius:5px; margin:5px 0;">
    <b>Supported formats:</b> PDB, MOL2, SDF, PDBQT<br>
    <b>Tip:</b> SDF files from PubChem/ChEMBL work well!
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
    
    return ligand_panel


# ==================== RECEPTOR PREPARATION (FIXED) ====================

def create_receptor_panel(state):
    """Create fixed receptor preparation panel"""
    
    receptor_upload = widgets.FileUpload(
        accept='.pdb,.pdbqt',
        multiple=False,
        description='Receptor PDB'
    )
    
    receptor_repairs = widgets.Dropdown(
        options=[
            ('Check and add missing hydrogens', 'checkhydrogens'),
            ('Add bonds and hydrogens', 'bonds_hydrogens'),
            ('Add hydrogens only', 'hydrogens'),
            ('No repairs', '')
        ],
        value='checkhydrogens',
        description='Repairs:'
    )
    
    receptor_charges = widgets.Dropdown(
        options=['gasteiger', 'Kollman'],
        value='gasteiger',
        description='Charges:'
    )
    
    receptor_cleanup = widgets.SelectMultiple(
        options=['nphs', 'lps', 'waters', 'nonstdres'],
        value=['nphs', 'lps', 'waters', 'nonstdres'],
        description='Cleanup:',
        rows=4
    )
    
    receptor_output = widgets.Output()
    receptor_3d_view = widgets.Output()
    
    prepare_receptor_btn = widgets.Button(
        description='Prepare Receptor',
        button_style='primary',
        icon='check'
    )
    
    def prepare_receptor(btn):
        with receptor_output:
            clear_output()
            if not receptor_upload.value:
                print('Please upload a receptor file first!')
                return
            
            # Save uploaded file
            uploaded_file = receptor_upload.value[0]
            filename = uploaded_file['name']
            content = uploaded_file['content']
            receptor_path = state.work_dir / filename
            
            with open(receptor_path, 'wb') as f:
                f.write(content)
            
            state.receptor_pdb = receptor_path
            print(f'Saved: {receptor_path}')
            
            # Prepare receptor
            try:
                mol = Read(str(receptor_path))[0]
                cleanup = '_'.join(receptor_cleanup.value) if receptor_cleanup.value else ''
                repairs = receptor_repairs.value
                
                output_pdbqt = state.work_dir / f'{receptor_path.stem}_prepared.pdbqt'
                
                print(f'Preparing... (this may take a moment for large proteins)')
                
                prep = AD4ReceptorPreparation(
                    mol,
                    mode='automatic',
                    repairs=repairs,
                    charges_to_add=receptor_charges.value,
                    cleanup=cleanup,
                    outputfilename=str(output_pdbqt)
                )
                
                state.receptor_pdbqt = output_pdbqt
                print(f'Receptor prepared: {output_pdbqt.name}')
                print(f'   Atoms: {len(mol.allAtoms)}')
                print(f'   Residues: {len(mol.chains.residues)}')
                
                # Calculate center of mass for grid
                coords = mol.allAtoms.coords
                center = np.mean(coords, axis=0)
                state.grid_center = center.tolist()
                print(f'   Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})')
                
            except Exception as e:
                print(f'Error: {e}')
                import traceback
                traceback.print_exc()
                return
        
        # 3D visualization
        with receptor_3d_view:
            clear_output()
            if state.receptor_pdbqt and state.receptor_pdbqt.exists():
                try:
                    with open(state.receptor_pdbqt, 'r') as f:
                        pdbqt_content = f.read()
                    
                    view = py3Dmol.view(width=600, height=400)
                    view.addModel(pdbqt_content, 'pdb')
                    view.setStyle({'cartoon': {'color': 'spectrum'}})
                    view.setBackgroundColor('0x1a1a2e')
                    view.zoomTo()
                    view.show()
                except Exception as e:
                    print(f'3D view error: {e}')
    
    prepare_receptor_btn.on_click(prepare_receptor)
    
    # Large file warning
    warning_text = widgets.HTML('''
    <div style="background:#4a3d2d; padding:10px; border-radius:5px; margin:5px 0;">
    <b>Tip:</b> For large proteins (>3000 atoms), consider using only one chain to avoid memory issues.
    </div>
    ''')
    
    # Layout
    receptor_options = VBox([
        widgets.HTML('<h4>Preparation Options</h4>'),
        receptor_repairs,
        receptor_charges,
        receptor_cleanup
    ])
    
    receptor_panel = VBox([
        widgets.HTML('<h3>Upload Receptor</h3>'),
        warning_text,
        receptor_upload,
        receptor_options,
        prepare_receptor_btn,
        receptor_output,
        widgets.HTML('<h4>3D Preview</h4>'),
        receptor_3d_view
    ])
    
    return receptor_panel


if __name__ == "__main__":
    print("This module provides fixed UI panels for the AutoDock Workbench.")
    print("Import and use create_receptor_panel() and create_ligand_panel() in a Jupyter notebook.")
