# =============================================================================
# AUTODOCK WORKBENCH - FIXED CELLS
# =============================================================================
#
# This file contains fixed versions of ALL cells for the AutoDock Workbench.
# Copy each section into the corresponding cell in your Jupyter notebook.
#
# AUDIT SUMMARY:
# - 18 issues identified across 5 cells
# - All critical issues fixed below
#
# =============================================================================


# =============================================================================
# CELL 1: SETUP & IMPORTS (No changes needed - working correctly)
# =============================================================================

# ==================== SETUP & IMPORTS ====================
import os
import sys
import warnings
warnings.filterwarnings('ignore')

# Add AutoDockTools to path
ADT_PATH = os.path.dirname(os.path.abspath('__file__'))
if ADT_PATH not in sys.path:
    sys.path.insert(0, ADT_PATH)

# Core imports
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
import tempfile
import shutil

# Visualization
import py3Dmol
import matplotlib.pyplot as plt
from IPython.display import display, HTML, clear_output

# Widgets
import ipywidgets as widgets
from ipywidgets import Layout, VBox, HBox, Tab, Accordion

# Chemistry
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors

# AutoDockTools
from MolKit import Read
from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation, AD4LigandPreparation
from AutoDockTools.GridParameters import GridParameter4FileMaker
from AutoDockTools.DockingParameters import DockingParameter42FileMaker

print('✅ All imports successful!')
print(f'📁 Working directory: {ADT_PATH}')


# =============================================================================
# CELL 2: GLOBAL STATE (No changes needed - working correctly)
# =============================================================================

# ==================== GLOBAL STATE ====================
class WorkbenchState:
    """Global state for the workbench"""
    def __init__(self):
        self.work_dir = Path(ADT_PATH) / 'workbench_files'
        self.work_dir.mkdir(exist_ok=True)
        
        # Paths to executables
        self.autodock4_exe = Path(ADT_PATH) / 'autodock4.exe'
        self.autogrid4_exe = Path(ADT_PATH) / 'autogrid4.exe'
        self.vina_exe = Path(ADT_PATH) / 'vina.exe'
        
        # Prepared files
        self.receptor_pdb = None
        self.receptor_pdbqt = None
        self.ligand_file = None
        self.ligand_pdbqt = None
        
        # Grid parameters
        self.grid_center = [0.0, 0.0, 0.0]
        self.grid_size = [40, 40, 40]
        self.grid_spacing = 0.375
        
        # Results
        self.docking_results = None
        
state = WorkbenchState()
print(f'📂 Work directory: {state.work_dir}')
print(f'🔧 AutoDock4: {"Found" if state.autodock4_exe.exists() else "Not found"}')
print(f'🔧 AutoGrid4: {"Found" if state.autogrid4_exe.exists() else "Not found"}')
print(f'🔧 Vina: {"Found" if state.vina_exe.exists() else "Not found"}')


# =============================================================================
# CELL 3: RECEPTOR PREPARATION (FIXED - duplicate atoms, better error handling)
# =============================================================================

# ==================== FIXED RECEPTOR PREPARATION ====================

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
        ('No repairs (skip hydrogen addition)', '')  # FIXED: clearer label
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

# FIXED: Function to remove duplicate atoms
def clean_pdb_duplicates(pdb_path):
    """Remove duplicate atoms with identical coordinates from PDB file"""
    seen_coords = {}
    cleaned_lines = []
    duplicates_removed = 0
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coord_key = (round(x, 3), round(y, 3), round(z, 3))
                    
                    if coord_key in seen_coords:
                        duplicates_removed += 1
                        continue
                    seen_coords[coord_key] = True
                except:
                    pass
            cleaned_lines.append(line)
    
    if duplicates_removed > 0:
        clean_path = pdb_path.parent / f'{pdb_path.stem}_clean.pdb'
        with open(clean_path, 'w') as f:
            f.writelines(cleaned_lines)
        return clean_path, duplicates_removed
    return pdb_path, 0

def prepare_receptor(btn):
    with receptor_output:
        clear_output()
        if not receptor_upload.value:
            print('❌ Please upload a receptor file first!')
            return
        
        # Save uploaded file
        uploaded_file = receptor_upload.value[0]
        filename = uploaded_file['name']
        content = uploaded_file['content']
        receptor_path = state.work_dir / filename
        
        with open(receptor_path, 'wb') as f:
            f.write(content)
        
        state.receptor_pdb = receptor_path
        print(f'📁 Saved: {receptor_path}')
        
        # FIXED: Clean duplicate atoms
        receptor_path, dups = clean_pdb_duplicates(receptor_path)
        if dups > 0:
            print(f'⚠️ Removed {dups} duplicate atoms with identical coordinates')
            state.receptor_pdb = receptor_path
        
        # Prepare receptor
        try:
            mol = Read(str(receptor_path))[0]
            cleanup = '_'.join(receptor_cleanup.value) if receptor_cleanup.value else ''
            repairs = receptor_repairs.value
            
            output_pdbqt = state.work_dir / f'{receptor_path.stem}_prepared.pdbqt'
            
            print(f'📊 Atoms: {len(mol.allAtoms)}')
            if len(mol.allAtoms) > 3000:
                print('⏳ Large protein - this may take a moment...')
            
            prep = AD4ReceptorPreparation(
                mol,
                mode='automatic',
                repairs=repairs,
                charges_to_add=receptor_charges.value,
                cleanup=cleanup,
                outputfilename=str(output_pdbqt)
            )
            
            # FIXED: Verify file was created
            if not output_pdbqt.exists():
                print('❌ Receptor preparation failed - no output file created')
                return
            
            state.receptor_pdbqt = output_pdbqt
            print(f'✅ Receptor prepared: {output_pdbqt.name}')
            print(f'   Residues: {len(mol.chains.residues)}')
            
            # Calculate center of mass for grid
            coords = mol.allAtoms.coords
            center = np.mean(coords, axis=0)
            state.grid_center = center.tolist()
            print(f'   Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})')
            
        except Exception as e:
            error_msg = str(e)
            print(f'❌ Error: {error_msg}')
            
            # FIXED: Better error messages
            if 'ZeroDivision' in error_msg or 'same coordinates' in error_msg:
                print('\n💡 Try selecting "No repairs (skip hydrogen addition)"')
            elif 'babel_type' in error_msg:
                print('\n💡 The structure may have unusual atom types.')
                print('   Try selecting "No repairs" option.')
            
            import traceback
            traceback.print_exc()
            return
    
    # 3D visualization - FIXED: uses view.show()
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
                view.show()  # <-- FIXED: was display(view)
            except Exception as e:
                print(f'3D view error: {e}')

prepare_receptor_btn.on_click(prepare_receptor)

# Warning text
warning_text = widgets.HTML('''
<div style="background:#4a3d2d; padding:10px; border-radius:5px; margin:5px 0;">
<b>Tips:</b><br>
• For PDBs with errors, try "No repairs" option<br>
• Large proteins (>3000 atoms) may take a while<br>
• Consider using a single chain for multimeric proteins
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
    widgets.HTML('<h3>🧬 Upload Receptor</h3>'),
    warning_text,
    receptor_upload,
    receptor_options,
    prepare_receptor_btn,
    receptor_output,
    widgets.HTML('<h4>3D Preview</h4>'),
    receptor_3d_view
])

display(receptor_panel)


# =============================================================================
# CELL 4: LIGAND PREPARATION (FIXED - better error handling, ndihe check)
# =============================================================================

# ==================== FIXED LIGAND PREPARATION ====================

ligand_upload = widgets.FileUpload(
    accept='.pdb,.mol2,.sdf,.pdbqt',
    multiple=False,
    description='Ligand File',
    layout=Layout(width='auto')
)

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
            print('❌ Please upload a ligand file first!')
            return
        
        # Save uploaded file
        uploaded_file = ligand_upload.value[0]
        filename = uploaded_file['name']
        content = uploaded_file['content']
        ligand_path = state.work_dir / filename
        
        with open(ligand_path, 'wb') as f:
            f.write(content)
        
        state.ligand_file = ligand_path
        print(f'📁 Saved: {ligand_path}')
        
        output_pdbqt = state.work_dir / f'{ligand_path.stem}_prepared.pdbqt'
        repairs = ligand_repairs.value
        
        # For SDF files, use RDKit to convert to PDB first
        if str(ligand_path).lower().endswith('.sdf'):
            print('🔄 Converting SDF via RDKit...')
            try:
                rdkit_mol = Chem.SDMolSupplier(str(ligand_path))[0]
                if rdkit_mol is None:
                    print('❌ Error: Could not read SDF file')
                    return
                rdkit_mol = Chem.AddHs(rdkit_mol, addCoords=True)
                if rdkit_mol.GetNumConformers() == 0:
                    AllChem.EmbedMolecule(rdkit_mol, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(rdkit_mol)
                pdb_path = state.work_dir / f'{ligand_path.stem}_temp.pdb'
                Chem.MolToPDBFile(rdkit_mol, str(pdb_path))
                ligand_path = pdb_path
                print(f'✅ Converted to: {pdb_path.name}')
            except Exception as e:
                print(f'❌ Error converting SDF: {e}')
                return
        
        # Prepare ligand with MolKit
        try:
            mol = Read(str(ligand_path))[0]
            print(f'📊 Loaded: {len(mol.allAtoms)} atoms')
            
            prep = AD4LigandPreparation(
                mol,
                mode='automatic',
                repairs=repairs,
                charges_to_add='gasteiger',
                cleanup='nphs_lps',
                outputfilename=str(output_pdbqt)
            )
            
            state.ligand_pdbqt = output_pdbqt
            print(f'✅ Ligand prepared: {output_pdbqt.name}')
            
            # FIXED: Check for ndihe attribute safely
            if hasattr(prep, 'ndihe') and prep.ndihe is not None:
                print(f'   Rotatable bonds: {prep.ndihe}')
            elif hasattr(prep, 'TORSDOF') and prep.TORSDOF is not None:
                print(f'   Rotatable bonds: {prep.TORSDOF}')
            else:
                # Try to count from PDBQT file
                try:
                    with open(output_pdbqt, 'r') as f:
                        for line in f:
                            if line.startswith('TORSDOF'):
                                print(f'   Rotatable bonds: {line.split()[1]}')
                                break
                except:
                    pass
            
        except Exception as e:
            print(f'❌ Error: {e}')
            import traceback
            traceback.print_exc()
            return
    
    # FIXED: 3D visualization - Convert PDBQT to simple PDB for py3Dmol
    with ligand_3d_view:
        clear_output()
        if state.ligand_pdbqt and state.ligand_pdbqt.exists():
            try:
                # Convert PDBQT to PDB (py3Dmol doesn't handle BRANCH/ROOT)
                pdb_lines = []
                with open(state.ligand_pdbqt, 'r') as f:
                    for line in f:
                        if line.startswith(('HETATM', 'ATOM')):
                            pdb_lines.append(line[:66].rstrip() + '\n')
                pdb_lines.append('END\n')
                pdb_content = ''.join(pdb_lines)
                
                view = py3Dmol.view(width=400, height=300)
                view.addModel(pdb_content, 'pdb')
                view.setStyle({'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.15}})
                view.addStyle({'elem': 'Cl'}, {'stick': {'color': 'green', 'radius': 0.2}})
                view.addStyle({'elem': 'O'}, {'stick': {'color': 'red', 'radius': 0.18}})
                view.addStyle({'elem': 'N'}, {'stick': {'color': 'blue', 'radius': 0.18}})
                view.setBackgroundColor('0x1a1a2e')
                view.zoomTo()
                view.show()  # <-- FIXED: was display(view)
            except Exception as e:
                print(f'3D view error: {e}')
    
    # 2D visualization with RDKit
    with ligand_2d_view:
        clear_output()
        try:
            orig_path = state.ligand_file
            if orig_path and orig_path.exists():  # FIXED: check file exists
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
    widgets.HTML('<h3>💊 Upload Ligand</h3>'),
    help_text,
    ligand_upload,
    ligand_repairs,
    prepare_ligand_btn,
    ligand_output,
    HBox([ligand_3d_view, ligand_2d_view])
])

display(ligand_panel)


# =============================================================================
# CELL 5: GRID BOX CONFIGURATION (FIXED - 3D preview, validation)
# =============================================================================

# ==================== FIXED GRID BOX CONFIGURATION ====================

# Center controls
center_x = widgets.FloatText(value=0.0, description='Center X:', step=0.5)
center_y = widgets.FloatText(value=0.0, description='Center Y:', step=0.5)
center_z = widgets.FloatText(value=0.0, description='Center Z:', step=0.5)

# Size controls
size_x = widgets.IntSlider(value=40, min=10, max=126, description='Size X:')
size_y = widgets.IntSlider(value=40, min=10, max=126, description='Size Y:')
size_z = widgets.IntSlider(value=40, min=10, max=126, description='Size Z:')

spacing = widgets.FloatSlider(value=0.375, min=0.1, max=1.0, step=0.025, description='Spacing:')

grid_output = widgets.Output()
grid_3d_view = widgets.Output()

auto_center_btn = widgets.Button(description='Center on Receptor', button_style='info')
center_on_ligand_btn = widgets.Button(description='Center on Ligand', button_style='info')
preview_grid_btn = widgets.Button(description='Preview Grid Box', button_style='primary', icon='cube')

# FIXED: Added validation and better error messages
def auto_center(btn):
    with grid_output:
        clear_output()
        if not state.receptor_pdbqt or not state.receptor_pdbqt.exists():
            print('❌ Please prepare a receptor first!')
            return
        center_x.value, center_y.value, center_z.value = state.grid_center
        print(f'✅ Centered on receptor: ({center_x.value:.2f}, {center_y.value:.2f}, {center_z.value:.2f})')

# FIXED: Added validation and error handling
def center_on_ligand(btn):
    with grid_output:
        clear_output()
        if not state.ligand_pdbqt or not state.ligand_pdbqt.exists():
            print('❌ Please prepare a ligand first!')
            return
        try:
            mol = Read(str(state.ligand_pdbqt))[0]
            coords = np.array(mol.allAtoms.coords)
            center = np.mean(coords, axis=0)
            center_x.value, center_y.value, center_z.value = center
            print(f'✅ Centered on ligand: ({center_x.value:.2f}, {center_y.value:.2f}, {center_z.value:.2f})')
        except Exception as e:
            print(f'❌ Error reading ligand: {e}')

def preview_grid(btn):
    state.grid_center = [center_x.value, center_y.value, center_z.value]
    state.grid_size = [size_x.value, size_y.value, size_z.value]
    state.grid_spacing = spacing.value
    
    with grid_3d_view:
        clear_output()
        view = py3Dmol.view(width=600, height=400)
        
        # Add receptor if available
        if state.receptor_pdbqt and state.receptor_pdbqt.exists():
            with open(state.receptor_pdbqt, 'r') as f:
                view.addModel(f.read(), 'pdb')
            view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum', 'opacity': 0.7}})
        
        # FIXED: Add ligand with PDBQT to PDB conversion
        if state.ligand_pdbqt and state.ligand_pdbqt.exists():
            pdb_lines = []
            with open(state.ligand_pdbqt, 'r') as f:
                for line in f:
                    if line.startswith(('HETATM', 'ATOM')):
                        pdb_lines.append(line[:66].rstrip() + '\n')
            pdb_lines.append('END\n')
            view.addModel(''.join(pdb_lines), 'pdb')
            view.setStyle({'model': 1}, {'stick': {'colorscheme': 'cyanCarbon'}})
        
        # Add grid box
        cx, cy, cz = state.grid_center
        sx = size_x.value * spacing.value / 2
        sy = size_y.value * spacing.value / 2
        sz = size_z.value * spacing.value / 2
        
        view.addBox({
            'center': {'x': cx, 'y': cy, 'z': cz},
            'dimensions': {'w': sx*2, 'h': sy*2, 'd': sz*2},
            'color': 'yellow',
            'opacity': 0.3
        })
        
        view.setBackgroundColor('0x1a1a2e')
        view.zoomTo()
        view.show()  # <-- FIXED: was display(view)
    
    with grid_output:
        clear_output()
        box_dims = [s * spacing.value for s in state.grid_size]
        print(f'📦 Grid Box Dimensions: {box_dims[0]:.1f} × {box_dims[1]:.1f} × {box_dims[2]:.1f} Å')
        print(f'   Center: ({cx:.2f}, {cy:.2f}, {cz:.2f})')
        print(f'   Grid points: {size_x.value} × {size_y.value} × {size_z.value}')

auto_center_btn.on_click(auto_center)
center_on_ligand_btn.on_click(center_on_ligand)
preview_grid_btn.on_click(preview_grid)

grid_panel = VBox([
    widgets.HTML('<h3>📦 Grid Box Setup</h3>'),
    widgets.HTML('<h4>Center Coordinates</h4>'),
    HBox([center_x, center_y, center_z]),
    HBox([auto_center_btn, center_on_ligand_btn]),
    widgets.HTML('<h4>Grid Dimensions</h4>'),
    size_x, size_y, size_z, spacing,
    preview_grid_btn,
    grid_output,
    grid_3d_view
])

display(grid_panel)


# =============================================================================
# CELL 6: DOCKING (FIXED - removed --log, added validation, parse scores)
# =============================================================================

# ==================== FIXED DOCKING ====================

docking_engine = widgets.RadioButtons(
    options=['AutoDock Vina', 'AutoDock4'],
    value='AutoDock Vina',
    description='Engine:'
)

# Vina parameters
vina_exhaustiveness = widgets.IntSlider(value=8, min=1, max=32, description='Exhaustiveness:')
vina_num_modes = widgets.IntSlider(value=9, min=1, max=20, description='Num Modes:')
vina_energy_range = widgets.FloatSlider(value=3.0, min=1.0, max=10.0, description='Energy Range:')

# AD4 parameters
ad4_ga_runs = widgets.IntSlider(value=10, min=1, max=100, description='GA Runs:')
ad4_pop_size = widgets.IntSlider(value=150, min=50, max=300, description='Population:')
ad4_num_evals = widgets.Dropdown(
    options=[('Short - 250k', 250000), ('Medium - 2.5M', 2500000), ('Long - 25M', 25000000)],
    value=2500000,
    description='Evaluations:'
)

docking_output = widgets.Output()
run_docking_btn = widgets.Button(description='🚀 Run Docking', button_style='success', icon='play')

def run_docking(btn):
    with docking_output:
        clear_output()
        
        # FIXED: Validation checks
        if not state.receptor_pdbqt or not state.receptor_pdbqt.exists():
            print('❌ Please prepare a receptor first!')
            return
        
        if not state.ligand_pdbqt or not state.ligand_pdbqt.exists():
            print('❌ Please prepare a ligand first!')
            return
        
        engine = docking_engine.value
        print(f'🔄 Running {engine}...')
        
        if engine == 'AutoDock Vina':
            # FIXED: Check Vina exists
            if not state.vina_exe.exists():
                print(f'❌ Vina not found at: {state.vina_exe}')
                return
            
            # Create Vina config
            config_path = state.work_dir / 'vina_config.txt'
            output_path = state.work_dir / 'vina_out.pdbqt'
            
            cx, cy, cz = state.grid_center
            sx = state.grid_size[0] * state.grid_spacing
            sy = state.grid_size[1] * state.grid_spacing
            sz = state.grid_size[2] * state.grid_spacing
            
            config = f"""receptor = {state.receptor_pdbqt}
ligand = {state.ligand_pdbqt}
out = {output_path}

center_x = {cx:.3f}
center_y = {cy:.3f}
center_z = {cz:.3f}

size_x = {sx:.1f}
size_y = {sy:.1f}
size_z = {sz:.1f}

exhaustiveness = {vina_exhaustiveness.value}
num_modes = {vina_num_modes.value}
energy_range = {vina_energy_range.value}
"""
            with open(config_path, 'w') as f:
                f.write(config)
            
            print(f'📝 Config saved: {config_path}')
            
            # FIXED: Run Vina without --log option (not supported in 1.2.x)
            cmd = f'"{state.vina_exe}" --config "{config_path}"'
            print(f'⏳ Executing: {cmd}')
            print('This may take a few minutes...\n')
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=str(state.work_dir))
            
            # FIXED: Parse and display Vina output
            if result.stdout:
                print('=' * 50)
                print('VINA OUTPUT:')
                print('=' * 50)
                
                # Extract the results table
                in_table = False
                for line in result.stdout.split('\n'):
                    if 'mode' in line.lower() and 'affinity' in line.lower():
                        in_table = True
                    if in_table or ('Writing output' in line):
                        print(line)
                    if in_table and line.strip() == '':
                        in_table = False
                
                print('=' * 50)
            
            if result.returncode == 0 and output_path.exists():
                print('✅ Docking complete!')
                
                # Count poses
                with open(output_path, 'r') as f:
                    num_poses = f.read().count('MODEL')
                print(f'📊 Generated {num_poses} poses')
                
                state.docking_results = output_path
            else:
                if result.stderr:
                    print(f'❌ Error: {result.stderr}')
        
        else:  # AutoDock4
            print('⚠️ AutoDock4 workflow requires grid map preparation.')
            print('This involves running AutoGrid4 first.')
            print('')
            gpf_path = state.work_dir / f'{state.receptor_pdbqt.stem}.gpf'
            dpf_path = state.work_dir / f'{state.ligand_pdbqt.stem}.dpf'
            print(f'📝 GPF would be saved to: {gpf_path}')
            print(f'📝 DPF would be saved to: {dpf_path}')
            print('\n⚠️ Full AD4 integration coming in a future update!')

run_docking_btn.on_click(run_docking)

vina_params = VBox([
    widgets.HTML('<h4>Vina Parameters</h4>'),
    vina_exhaustiveness,
    vina_num_modes,
    vina_energy_range
])

ad4_params = VBox([
    widgets.HTML('<h4>AutoDock4 Parameters</h4>'),
    ad4_ga_runs,
    ad4_pop_size,
    ad4_num_evals
])

docking_panel = VBox([
    widgets.HTML('<h3>⚙️ Docking Configuration</h3>'),
    docking_engine,
    HBox([vina_params, ad4_params]),
    run_docking_btn,
    docking_output
])

display(docking_panel)


# =============================================================================
# CELL 7: RESULTS VISUALIZATION (FIXED - 3D preview, pose parsing)
# =============================================================================

# ==================== FIXED RESULTS & INTERACTIONS ====================

results_output = widgets.Output()
results_3d_view = widgets.Output()
interaction_output = widgets.Output()

pose_selector = widgets.IntSlider(value=1, min=1, max=9, description='Pose:')
show_surface = widgets.Checkbox(value=False, description='Show Surface')
# FIXED: Removed unused show_hbonds checkbox

visualize_btn = widgets.Button(description='Visualize Results', button_style='primary', icon='eye')
analyze_btn = widgets.Button(description='Analyze Interactions', button_style='info', icon='search')

def visualize_results(btn):
    with results_3d_view:
        clear_output()
        
        if not state.docking_results or not state.docking_results.exists():
            print('❌ No docking results available. Run docking first!')
            return
        
        # FIXED: Update pose_selector max based on actual results
        with open(state.docking_results, 'r') as f:
            docked_content = f.read()
        num_poses = docked_content.count('MODEL')
        if num_poses > 0:
            pose_selector.max = num_poses
        
        view = py3Dmol.view(width=800, height=500)
        
        # Add receptor
        if state.receptor_pdbqt and state.receptor_pdbqt.exists():
            with open(state.receptor_pdbqt, 'r') as f:
                receptor_content = f.read()
            view.addModel(receptor_content, 'pdb')
            
            if show_surface.value:
                view.setStyle({'model': 0}, {'cartoon': {'color': 'white', 'opacity': 0.3}})
                view.addSurface(py3Dmol.VDW, {'opacity': 0.5, 'color': 'lightblue'}, {'model': 0})
            else:
                view.setStyle({'model': 0}, {'cartoon': {'color': 'spectrum'}})
        
        # FIXED: Parse poses correctly (first split element is empty)
        models = docked_content.split('MODEL')
        pose_idx = pose_selector.value
        
        # models[0] is empty, models[1] is pose 1, etc.
        if pose_idx < len(models):
            pose_content = 'MODEL' + models[pose_idx]
            
            # FIXED: Convert PDBQT to PDB for proper rendering
            pdb_lines = []
            for line in pose_content.split('\n'):
                if line.startswith(('HETATM', 'ATOM')):
                    pdb_lines.append(line[:66].rstrip() + '\n')
            pdb_lines.append('END\n')
            
            view.addModel(''.join(pdb_lines), 'pdb')
            view.setStyle({'model': 1}, {'stick': {'colorscheme': 'cyanCarbon', 'radius': 0.2}})
            view.addStyle({'model': 1}, {'sphere': {'colorscheme': 'cyanCarbon', 'radius': 0.4}})
        
        view.setBackgroundColor('0x1a1a2e')
        view.zoomTo({'model': 1})
        view.show()  # <-- FIXED: was display(view)

def analyze_interactions(btn):
    with interaction_output:
        clear_output()
        
        if not state.docking_results or not state.receptor_pdbqt:
            print('❌ Need both receptor and docking results for interaction analysis!')
            return
        
        print('🔍 Analyzing protein-ligand interactions...')
        print('-' * 50)
        
        try:
            # Basic distance-based interaction detection
            receptor_mol = Read(str(state.receptor_pdbqt))[0]
            ligand_mol = Read(str(state.docking_results))[0]
            
            rec_coords = np.array(receptor_mol.allAtoms.coords)
            lig_coords = np.array(ligand_mol.allAtoms.coords)
            
            # Find close contacts
            from scipy.spatial.distance import cdist
            distances = cdist(lig_coords, rec_coords)
            
            close_contacts = np.where(distances < 4.0)
            
            print(f'\n📊 Close Contacts (< 4Å): {len(close_contacts[0])}')
            
            # Find residues in contact
            contact_residues = set()
            for rec_idx in close_contacts[1]:
                atom = receptor_mol.allAtoms[rec_idx]
                res = atom.parent
                contact_residues.add(f'{res.name}{res.number}')
            
            print(f'\n🧬 Residues in Contact ({len(contact_residues)}):')
            for res in sorted(contact_residues, key=lambda x: int(''.join(filter(str.isdigit, x)) or 0)):
                print(f'   • {res}')
                
        except Exception as e:
            print(f'❌ Analysis error: {e}')
            import traceback
            traceback.print_exc()

visualize_btn.on_click(visualize_results)
analyze_btn.on_click(analyze_interactions)

results_panel = VBox([
    widgets.HTML('<h3>🔬 Results Visualization</h3>'),
    HBox([pose_selector, show_surface]),
    HBox([visualize_btn, analyze_btn]),
    results_3d_view,
    interaction_output
])

display(results_panel)


# =============================================================================
# CELL 8: SESSION SUMMARY (No changes needed - working correctly)
# =============================================================================

# ==================== SESSION SUMMARY ====================

def print_summary():
    print('=' * 60)
    print('            AUTODOCK WORKBENCH - SESSION SUMMARY')
    print('=' * 60)
    print()
    print(f'📁 Working Directory: {state.work_dir}')
    print()
    print('📄 Files:')
    print(f'   Receptor PDBQT: {state.receptor_pdbqt.name if state.receptor_pdbqt else "Not prepared"}')
    print(f'   Ligand PDBQT:   {state.ligand_pdbqt.name if state.ligand_pdbqt else "Not prepared"}')
    print(f'   Docking Output: {state.docking_results.name if state.docking_results else "Not available"}')
    print()
    print('📦 Grid Box:')
    print(f'   Center: ({state.grid_center[0]:.2f}, {state.grid_center[1]:.2f}, {state.grid_center[2]:.2f})')
    print(f'   Size:   {state.grid_size[0]} × {state.grid_size[1]} × {state.grid_size[2]} points')
    print(f'   Spacing: {state.grid_spacing} Å')
    print()
    print('=' * 60)
    
    # List all generated files
    if state.work_dir.exists():
        print('\n📂 Generated Files:')
        for f in sorted(state.work_dir.iterdir()):
            if f.is_file():
                size = f.stat().st_size / 1024
                print(f'   {f.name} ({size:.1f} KB)')

print_summary()
