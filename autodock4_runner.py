# autodock4_runner.py
"""
Complete AutoDock4 workflow implementation.
Import this in the notebook and call run_autodock4() to execute.
"""

import subprocess
from pathlib import Path
import numpy as np


def get_ligand_atom_types(pdbqt_path):
    """Extract unique atom types from ligand PDBQT file."""
    atom_types = set()
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                parts = line.split()
                if len(parts) >= 12:
                    atom_type = parts[-1]
                    # Only include valid AutoDock atom types
                    valid_types = ['A', 'C', 'N', 'O', 'S', 'H', 'HD', 'HS', 'P', 'F', 'Cl', 'Br', 'I', 
                                   'NA', 'OA', 'SA', 'NS', 'OS', 'Mg', 'Mn', 'Zn', 'Ca', 'Fe']
                    if atom_type in valid_types:
                        atom_types.add(atom_type)
    return sorted(list(atom_types))


def run_autodock4(state, ga_runs=10, pop_size=150, num_evals=2500000):
    """
    Run complete AutoDock4 docking workflow.
    
    Args:
        state: WorkbenchState object with receptor_pdbqt, ligand_pdbqt, grid_center, etc.
        ga_runs: Number of genetic algorithm runs (default: 10)
        pop_size: GA population size (default: 150)
        num_evals: Number of energy evaluations (default: 2.5M)
    
    Returns:
        Path to docked output file, or None if failed.
    """
    print('=' * 50)
    print('AUTODOCK4 WORKFLOW')
    print('=' * 50)
    
    # Check executables exist
    if not state.autodock4_exe.exists():
        print(f'[ERROR] AutoDock4 not found at: {state.autodock4_exe}')
        return None
    if not state.autogrid4_exe.exists():
        print(f'[ERROR] AutoGrid4 not found at: {state.autogrid4_exe}')
        return None
    
    # Get ligand atom types for grid maps
    ligand_types = get_ligand_atom_types(state.ligand_pdbqt)
    print(f'Ligand atom types: {" ".join(ligand_types)}')
    
    # File paths
    receptor_stem = state.receptor_pdbqt.stem
    ligand_stem = state.ligand_pdbqt.stem
    gpf_path = state.work_dir / f'{receptor_stem}.gpf'
    glg_path = state.work_dir / f'{receptor_stem}.glg'
    dpf_path = state.work_dir / f'{ligand_stem}_{receptor_stem}.dpf'
    dlg_path = state.work_dir / f'{ligand_stem}_{receptor_stem}.dlg'
    
    # Grid parameters
    cx, cy, cz = state.grid_center
    npts = state.grid_size
    spacing = state.grid_spacing
    
    # ============== STEP 1: Create GPF (Grid Parameter File) ==============
    print('\n[1/4] Creating Grid Parameter File (GPF)...')
    
    # Build ligand type maps
    ligand_maps = '\n'.join([f'map {receptor_stem}.{atype}.map' for atype in ligand_types])
    
    gpf_content = f"""npts {npts[0]} {npts[1]} {npts[2]}
gridfld {receptor_stem}.maps.fld
spacing {spacing:.4f}
receptor_types A C HD N NA OA SA
ligand_types {" ".join(ligand_types)}
receptor {state.receptor_pdbqt.name}
gridcenter {cx:.4f} {cy:.4f} {cz:.4f}
smooth 0.5
{ligand_maps}
elecmap {receptor_stem}.e.map
dsolvmap {receptor_stem}.d.map
dielectric -0.1465
"""
    
    with open(gpf_path, 'w') as f:
        f.write(gpf_content)
    print(f'   GPF saved: {gpf_path.name}')
    
    # ============== STEP 2: Run AutoGrid4 ==============
    print('\n[2/4] Running AutoGrid4 (generating grid maps)...')
    print('   This may take several minutes for large receptors...')
    
    grid_cmd = f'"{state.autogrid4_exe}" -p "{gpf_path}" -l "{glg_path}"'
    grid_result = subprocess.run(grid_cmd, shell=True, capture_output=True, text=True, cwd=str(state.work_dir))
    
    if grid_result.returncode != 0:
        print(f'   [ERROR] AutoGrid4 failed!')
        if grid_result.stderr:
            print(f'   Error: {grid_result.stderr}')
        # Check log file for errors
        if glg_path.exists():
            with open(glg_path, 'r') as f:
                log_content = f.read()
            if 'error' in log_content.lower() or 'fatal' in log_content.lower():
                print('   Last 20 lines of log:')
                for line in log_content.split('\n')[-20:]:
                    print(f'   {line}')
        return None
    
    print('   AutoGrid4 completed successfully!')
    
    # ============== STEP 3: Create DPF (Docking Parameter File) ==============
    print('\n[3/4] Creating Docking Parameter File (DPF)...')
    
    # Build ligand type maps for DPF
    dpf_maps = '\n'.join([f'map {receptor_stem}.{atype}.map' for atype in ligand_types])
    
    dpf_content = f"""autodock_parameter_version 4.2
outlev 1
intelec
seed pid time
ligand_types {" ".join(ligand_types)}
fld {receptor_stem}.maps.fld
{dpf_maps}
elecmap {receptor_stem}.e.map
desolvmap {receptor_stem}.d.map
move {state.ligand_pdbqt.name}
about {cx:.4f} {cy:.4f} {cz:.4f}
tran0 random
quaternion0 random
dihe0 random
torsdof 6
rmstol 2.0
extnrg 1000.0
e0max 0.0 10000
ga_pop_size {pop_size}
ga_num_evals {num_evals}
ga_num_generations 27000
ga_elitism 1
ga_mutation_rate 0.02
ga_crossover_rate 0.8
ga_window_size 10
ga_cauchy_alpha 0.0
ga_cauchy_beta 1.0
set_ga
sw_max_its 300
sw_max_succ 4
sw_max_fail 4
sw_rho 1.0
sw_lb_rho 0.01
ls_search_freq 0.06
set_psw1
unbound_model bound
ga_run {ga_runs}
analysis
"""
    
    with open(dpf_path, 'w') as f:
        f.write(dpf_content)
    print(f'   DPF saved: {dpf_path.name}')
    
    # ============== STEP 4: Run AutoDock4 ==============
    print(f'\n[4/4] Running AutoDock4 ({ga_runs} GA runs)...')
    print('   This may take 10-30 minutes depending on settings...')
    
    dock_cmd = f'"{state.autodock4_exe}" -p "{dpf_path}" -l "{dlg_path}"'
    dock_result = subprocess.run(dock_cmd, shell=True, capture_output=True, text=True, cwd=str(state.work_dir))
    
    if dock_result.returncode != 0:
        print(f'   [ERROR] AutoDock4 failed!')
        if dock_result.stderr:
            print(f'   Error: {dock_result.stderr}')
        return None
    
    print('   AutoDock4 completed successfully!')
    
    # ============== Parse Results ==============
    print('\n' + '=' * 50)
    print('AUTODOCK4 RESULTS')
    print('=' * 50)
    
    output_pdbqt = None
    
    if dlg_path.exists():
        with open(dlg_path, 'r') as f:
            dlg_content = f.read()
        
        # Extract binding energies
        energies = []
        for line in dlg_content.split('\n'):
            if 'Estimated Free Energy of Binding' in line:
                try:
                    energy = float(line.split('=')[1].split('kcal')[0].strip())
                    energies.append(energy)
                except:
                    pass
        
        if energies:
            print(f'\nBinding Energies (kcal/mol):')
            for i, e in enumerate(energies[:10], 1):
                print(f'   Run {i}: {e:.2f}')
            print(f'\nBest binding energy: {min(energies):.2f} kcal/mol')
            print(f'Number of runs: {len(energies)}')
        
        # Convert DLG to PDBQT for visualization
        output_pdbqt = state.work_dir / f'{ligand_stem}_docked.pdbqt'
        models = []
        in_model = False
        current_model = []
        model_num = 1
        
        for line in dlg_content.split('\n'):
            if 'DOCKED: MODEL' in line:
                in_model = True
                current_model = [f'MODEL {model_num}\n']
                model_num += 1
            elif 'DOCKED: ENDMDL' in line:
                current_model.append('ENDMDL\n')
                models.append(''.join(current_model))
                in_model = False
            elif in_model and line.startswith('DOCKED:'):
                # Remove 'DOCKED: ' prefix
                pdbqt_line = line[8:]
                current_model.append(pdbqt_line + '\n')
        
        if models:
            with open(output_pdbqt, 'w') as f:
                f.writelines(models)
            print(f'\nDocked poses saved to: {output_pdbqt.name}')
        else:
            print('\nNo docked poses found in output.')
            output_pdbqt = dlg_path
    
    print('\n[COMPLETE] AutoDock4 docking finished!')
    return output_pdbqt
