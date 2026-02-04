#!/usr/bin/env python
"""
Test Docking Script for AutoDock Workbench

This script tests the complete docking workflow:
1. Receptor preparation
2. Ligand preparation (using RDKit + meeko for PDBQT conversion)
3. Grid box setup
4. Vina docking
"""

import os
import sys
import numpy as np
from pathlib import Path

# Fix Windows console encoding
if sys.platform == 'win32':
    import io
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding='utf-8', errors='replace')
    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding='utf-8', errors='replace')

# Add AutoDockTools to path
ADT_PATH = Path(__file__).parent
if str(ADT_PATH) not in sys.path:
    sys.path.insert(0, str(ADT_PATH))

print("=" * 60)
print("AutoDock Workbench - Docking Test Script")
print("=" * 60)

# Working directory
WORK_DIR = ADT_PATH / 'workbench_files'
WORK_DIR.mkdir(exist_ok=True)

# Executables
VINA_EXE = ADT_PATH / 'vina.exe'

print(f"\n[DIR] Working directory: {WORK_DIR}")
print(f"[VINA] {'Found' if VINA_EXE.exists() else 'Not found'}")

# ============= STEP 1: Check/Prepare Receptor =============
print("\n" + "=" * 60)
print("STEP 1: Receptor Preparation")
print("=" * 60)

# Use Chain A only (smaller, avoids memory issues)
receptor_pdb = WORK_DIR / '2B35_A.pdb'
receptor_pdbqt = WORK_DIR / '2B35_A_prepared.pdbqt'

# If receptor_pdbqt doesn't exist, create it
if not receptor_pdbqt.exists():
    if not receptor_pdb.exists():
        # Extract Chain A from full PDB
        full_pdb = WORK_DIR / '2B35.pdb'
        if full_pdb.exists():
            print(f"[FILE] Extracting Chain A from {full_pdb.name}...")
            with open(full_pdb, 'r') as f, open(receptor_pdb, 'w') as out:
                for line in f:
                    if line.startswith(('ATOM', 'HETATM')):
                        chain = line[21] if len(line) > 21 else ''
                        if chain == 'A':
                            out.write(line)
                    elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'END')):
                        out.write(line)
            print(f"[OK] Created: {receptor_pdb}")
        else:
            print(f"[ERROR] No receptor PDB found!")
            sys.exit(1)
    
    print(f"[FILE] Preparing receptor: {receptor_pdb.name}")
    try:
        from MolKit import Read
        from AutoDockTools.MoleculePreparation import AD4ReceptorPreparation
        
        mol = Read(str(receptor_pdb))[0]
        print(f"   Loaded: {len(mol.allAtoms)} atoms")
        
        prep = AD4ReceptorPreparation(
            mol,
            mode='automatic',
            repairs='checkhydrogens',
            charges_to_add='gasteiger',
            cleanup='nphs_lps_waters_nonstdres',
            outputfilename=str(receptor_pdbqt)
        )
        
        print(f"[OK] Receptor prepared: {receptor_pdbqt.name}")
        
        # Calculate center for grid
        coords = np.array(mol.allAtoms.coords)
        center = np.mean(coords, axis=0)
        print(f"   Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")
        
    except Exception as e:
        print(f"[ERROR] Receptor preparation failed: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
else:
    print(f"[OK] Receptor already prepared: {receptor_pdbqt.name}")
    # Calculate center from existing file
    coords = []
    with open(receptor_pdbqt, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except:
                    pass
    center = np.mean(coords, axis=0)
    print(f"   Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})")

# ============= STEP 2: Prepare Ligand =============
print("\n" + "=" * 60)
print("STEP 2: Ligand Preparation")
print("=" * 60)

ligand_sdf = WORK_DIR / 'TCL_ideal (1).sdf'
ligand_pdbqt = WORK_DIR / 'TCL_prepared.pdbqt'

if not ligand_sdf.exists():
    print(f"[ERROR] Ligand SDF not found: {ligand_sdf}")
    sys.exit(1)

print(f"[FILE] Preparing ligand: {ligand_sdf.name}")

# Try using RDKit + meeko for ligand preparation (more reliable than MolKit for small molecules)
try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
    
    # Read ligand
    mol = Chem.SDMolSupplier(str(ligand_sdf))[0]
    if mol is None:
        raise ValueError("Could not read SDF file")
    
    # Add hydrogens if needed
    mol = Chem.AddHs(mol, addCoords=True)
    
    # Generate 3D coordinates if not present
    if mol.GetNumConformers() == 0:
        AllChem.EmbedMolecule(mol, randomSeed=42)
        AllChem.MMFFOptimizeMolecule(mol)
    
    print(f"   Loaded: {mol.GetNumAtoms()} atoms")
    
    # Try meeko for PDBQT conversion
    try:
        from meeko import MoleculePreparation, PDBQTWriterLegacy
        
        preparator = MoleculePreparation()
        mol_prepared = preparator.prepare(mol)[0]
        
        pdbqt_string, _, _ = PDBQTWriterLegacy.write_string(mol_prepared)
        
        with open(ligand_pdbqt, 'w') as f:
            f.write(pdbqt_string)
        
        print(f"[OK] Ligand prepared (meeko): {ligand_pdbqt.name}")
        
    except ImportError:
        print("   meeko not available, using alternative method...")
        
        # Fallback: Convert to PDB first, then use MolKit
        ligand_pdb = WORK_DIR / 'TCL_temp.pdb'
        Chem.MolToPDBFile(mol, str(ligand_pdb))
        
        from MolKit import Read
        from AutoDockTools.MoleculePreparation import AD4LigandPreparation
        
        lig_mol = Read(str(ligand_pdb))[0]
        
        prep = AD4LigandPreparation(
            lig_mol,
            mode='automatic',
            repairs='bonds_hydrogens',
            charges_to_add='gasteiger',
            cleanup='nphs_lps',
            outputfilename=str(ligand_pdbqt)
        )
        
        print(f"[OK] Ligand prepared (MolKit): {ligand_pdbqt.name}")
        
except Exception as e:
    print(f"[ERROR] RDKit-based preparation failed: {e}")
    print("   Trying MolKit direct approach...")
    
    try:
        from MolKit import Read
        from AutoDockTools.MoleculePreparation import AD4LigandPreparation
        
        # MolKit can read SDF via mol2 format
        mol = Read(str(ligand_sdf))[0]
        
        prep = AD4LigandPreparation(
            mol,
            mode='automatic',
            repairs='bonds_hydrogens',
            charges_to_add='gasteiger',
            cleanup='nphs_lps',
            outputfilename=str(ligand_pdbqt)
        )
        
        print(f"[OK] Ligand prepared (MolKit): {ligand_pdbqt.name}")
        
    except Exception as e2:
        print(f"[ERROR] MolKit preparation also failed: {e2}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

# Get ligand center
lig_coords = []
with open(ligand_pdbqt, 'r') as f:
    for line in f:
        if line.startswith(('ATOM', 'HETATM')):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                lig_coords.append([x, y, z])
            except:
                pass
lig_center = np.mean(lig_coords, axis=0) if lig_coords else center

# ============= STEP 3: Grid Box Setup =============
print("\n" + "=" * 60)
print("STEP 3: Grid Box Setup")
print("=" * 60)

# Use the RECEPTOR center (binding site location) not the ligand's original position
# The crystallographic ligand in 2B35 is at the binding site!
grid_center = center  # Use receptor center for binding site
grid_size = [25, 25, 25]  # Angstroms

print(f"   Center: ({grid_center[0]:.2f}, {grid_center[1]:.2f}, {grid_center[2]:.2f})")
print(f"   Size: {grid_size[0]} x {grid_size[1]} x {grid_size[2]} Angstrom")

# ============= STEP 4: Run Vina =============
print("\n" + "=" * 60)
print("STEP 4: Vina Docking")
print("=" * 60)

# Create Vina config
vina_config = WORK_DIR / 'vina_config.txt'
vina_output = WORK_DIR / 'vina_out.pdbqt'

config_content = f"""receptor = {receptor_pdbqt}
ligand = {ligand_pdbqt}
out = {vina_output}

center_x = {grid_center[0]:.3f}
center_y = {grid_center[1]:.3f}
center_z = {grid_center[2]:.3f}

size_x = {grid_size[0]}
size_y = {grid_size[1]}
size_z = {grid_size[2]}

exhaustiveness = 8
num_modes = 9
energy_range = 3
"""

with open(vina_config, 'w') as f:
    f.write(config_content)

print(f"[NOTE] Config saved: {vina_config.name}")

# Run Vina
import subprocess

cmd = f'"{VINA_EXE}" --config "{vina_config}"'
print(f"[WAIT] Running: vina.exe --config vina_config.txt")
print("   This may take a few minutes...")

result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=str(WORK_DIR))

# Vina outputs results to stdout
print("\n" + "-" * 50)
print("[RESULT] VINA OUTPUT")
print("-" * 50)
print(result.stdout)

if result.returncode == 0 and vina_output.exists():
    print("\n[OK] Docking complete!")
    print(f"[FILE] Output: {vina_output.name}")
    
    # Count poses
    with open(vina_output, 'r') as f:
        poses = f.read().count('MODEL')
    print(f"   Generated {poses} poses")
else:
    print(f"\n[ERROR] Docking may have issues!")
    if result.stderr:
        print(f"   Error: {result.stderr}")

print("\n" + "=" * 60)
print("DONE!")
print("=" * 60)


