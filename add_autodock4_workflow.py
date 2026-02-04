#!/usr/bin/env python3
"""
Script to add complete AutoDock4 workflow to the notebook.
This replaces the placeholder with a full implementation.
"""

import json
from pathlib import Path

def add_autodock4_workflow(notebook_path):
    """Replace the AutoDock4 placeholder with full implementation."""
    
    # Read the notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        notebook = json.load(f)
    
    # Find the docking cell (Cell 6)
    for i, cell in enumerate(notebook.get('cells', [])):
        if cell.get('cell_type') == 'code':
            source = ''.join(cell.get('source', []))
            if 'CELL 6: DOCKING' in source and 'Full AD4 integration coming' in source:
                print(f"Found docking cell at index {i}")
                
                # Build the new source with complete AutoDock4 implementation
                new_source = '''# =============================================================================
# CELL 6: DOCKING (COMPLETE - Both Vina and AutoDock4 fully implemented)
# =============================================================================

# ==================== COMPLETE DOCKING IMPLEMENTATION ====================

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
run_docking_btn = widgets.Button(description='Run Docking', button_style='success', icon='play')

def get_ligand_atom_types(pdbqt_path):
    """Extract unique atom types from ligand PDBQT file."""
    atom_types = set()
    with open(pdbqt_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                # Atom type is in the last column of PDBQT
                parts = line.split()
                if len(parts) >= 12:
                    atom_type = parts[-1]
                    # Only include valid AutoDock atom types
                    if atom_type in ['A', 'C', 'N', 'O', 'S', 'H', 'HD', 'HS', 'P', 'F', 'Cl', 'Br', 'I', 
                                      'NA', 'OA', 'SA', 'NS', 'OS', 'Mg', 'Mn', 'Zn', 'Ca', 'Fe']:
                        atom_types.add(atom_type)
    return sorted(list(atom_types))

def run_docking(btn):
    with docking_output:
        clear_output()
        
        # Validation checks
        if not state.receptor_pdbqt or not state.receptor_pdbqt.exists():
            print('Please prepare a receptor first!')
            return
        
        if not state.ligand_pdbqt or not state.ligand_pdbqt.exists():
            print('Please prepare a ligand first!')
            return
        
        engine = docking_engine.value
        print(f'Running {engine}...')
        
        if engine == 'AutoDock Vina':
            # Check Vina exists
            if not state.vina_exe.exists():
                print(f'Vina not found at: {state.vina_exe}')
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
            
            print(f'Config saved: {config_path}')
            
            cmd = f'"{state.vina_exe}" --config "{config_path}"'
            print(f'Executing: {cmd}')
            print('This may take a few minutes...\\n')
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=str(state.work_dir))
            
            if result.stdout:
                print('=' * 50)
                print('VINA OUTPUT:')
                print('=' * 50)
                
                in_table = False
                for line in result.stdout.split('\\n'):
                    if 'mode' in line.lower() and 'affinity' in line.lower():
                        in_table = True
                    if in_table or ('Writing output' in line):
                        print(line)
                    if in_table and line.strip() == '':
                        in_table = False
                
                print('=' * 50)
            
            if result.returncode == 0 and output_path.exists():
                print('Docking complete!')
                
                with open(output_path, 'r') as f:
                    num_poses = f.read().count('MODEL')
                print(f'Generated {num_poses} poses')
                
                state.docking_results = output_path
            else:
                if result.stderr:
                    print(f'Error: {result.stderr}')
        
        else:  # AutoDock4 - FULL IMPLEMENTATION
            print('=' * 50)
            print('AUTODOCK4 WORKFLOW')
            print('=' * 50)
            
            # Check executables exist
            if not state.autodock4_exe.exists():
                print(f'AutoDock4 not found at: {state.autodock4_exe}')
                return
            if not state.autogrid4_exe.exists():
                print(f'AutoGrid4 not found at: {state.autogrid4_exe}')
                return
            
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
            print('\\n[1/4] Creating Grid Parameter File (GPF)...')
            
            gpf_content = f"""npts {npts[0]} {npts[1]} {npts[2]}
gridfld {receptor_stem}.maps.fld
spacing {spacing:.4f}
receptor_types A C N O H HD HS OA NA SA
ligand_types {" ".join(ligand_types)}
receptor {state.receptor_pdbqt.name}
gridcenter {cx:.4f} {cy:.4f} {cz:.4f}
smooth 0.5
map {receptor_stem}.A.map
map {receptor_stem}.C.map
map {receptor_stem}.N.map
map {receptor_stem}.O.map
map {receptor_stem}.HD.map
map {receptor_stem}.OA.map
map {receptor_stem}.NA.map
elecmap {receptor_stem}.e.map
dsolvmap {receptor_stem}.d.map
dielectric -0.1465
"""
            # Add maps for each ligand atom type
            for atype in ligand_types:
                if atype not in ['A', 'C', 'N', 'O', 'HD', 'OA', 'NA']:
                    gpf_content += f"map {receptor_stem}.{atype}.map\\n"
            
            with open(gpf_path, 'w') as f:
                f.write(gpf_content)
            print(f'   GPF saved: {gpf_path.name}')
            
            # ============== STEP 2: Run AutoGrid4 ==============
            print('\\n[2/4] Running AutoGrid4 (generating grid maps)...')
            print('   This may take several minutes for large receptors...')
            
            grid_cmd = f'"{state.autogrid4_exe}" -p "{gpf_path}" -l "{glg_path}"'
            grid_result = subprocess.run(grid_cmd, shell=True, capture_output=True, text=True, cwd=str(state.work_dir))
            
            if grid_result.returncode != 0:
                print(f'   AutoGrid4 failed!')
                if grid_result.stderr:
                    print(f'   Error: {grid_result.stderr}')
                # Check log file for errors
                if glg_path.exists():
                    with open(glg_path, 'r') as f:
                        log_content = f.read()
                    if 'error' in log_content.lower() or 'fatal' in log_content.lower():
                        for line in log_content.split('\\n')[-20:]:
                            print(f'   {line}')
                return
            
            print('   AutoGrid4 completed successfully!')
            
            # ============== STEP 3: Create DPF (Docking Parameter File) ==============
            print('\\n[3/4] Creating Docking Parameter File (DPF)...')
            
            dpf_content = f"""autodock_parameter_version 4.2
outlev 1
intelec
seed pid time
ligand_types {" ".join(ligand_types)}
fld {receptor_stem}.maps.fld
map {receptor_stem}.A.map
map {receptor_stem}.C.map
map {receptor_stem}.N.map
map {receptor_stem}.O.map
map {receptor_stem}.HD.map
map {receptor_stem}.OA.map
map {receptor_stem}.NA.map
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
ga_pop_size {ad4_pop_size.value}
ga_num_evals {ad4_num_evals.value}
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
ga_run {ad4_ga_runs.value}
analysis
"""
            # Add maps for each ligand atom type
            for atype in ligand_types:
                if atype not in ['A', 'C', 'N', 'O', 'HD', 'OA', 'NA']:
                    at_idx = dpf_content.find('elecmap')
                    dpf_content = dpf_content[:at_idx] + f"map {receptor_stem}.{atype}.map\\n" + dpf_content[at_idx:]
            
            with open(dpf_path, 'w') as f:
                f.write(dpf_content)
            print(f'   DPF saved: {dpf_path.name}')
            
            # ============== STEP 4: Run AutoDock4 ==============
            print(f'\\n[4/4] Running AutoDock4 ({ad4_ga_runs.value} GA runs)...')
            print('   This may take 10-30 minutes depending on settings...')
            
            dock_cmd = f'"{state.autodock4_exe}" -p "{dpf_path}" -l "{dlg_path}"'
            dock_result = subprocess.run(dock_cmd, shell=True, capture_output=True, text=True, cwd=str(state.work_dir))
            
            if dock_result.returncode != 0:
                print(f'   AutoDock4 failed!')
                if dock_result.stderr:
                    print(f'   Error: {dock_result.stderr}')
                return
            
            print('   AutoDock4 completed successfully!')
            
            # ============== Parse Results ==============
            print('\\n' + '=' * 50)
            print('AUTODOCK4 RESULTS')
            print('=' * 50)
            
            if dlg_path.exists():
                with open(dlg_path, 'r') as f:
                    dlg_content = f.read()
                
                # Extract binding energies
                energies = []
                for line in dlg_content.split('\\n'):
                    if 'Estimated Free Energy of Binding' in line:
                        try:
                            energy = float(line.split('=')[1].split('kcal')[0].strip())
                            energies.append(energy)
                        except:
                            pass
                
                if energies:
                    print(f'\\nBinding Energies (kcal/mol):')
                    for i, e in enumerate(energies[:10], 1):
                        print(f'   Run {i}: {e:.2f}')
                    print(f'\\nBest binding energy: {min(energies):.2f} kcal/mol')
                    print(f'Number of runs: {len(energies)}')
                
                # Convert DLG to PDBQT for visualization
                output_pdbqt = state.work_dir / f'{ligand_stem}_docked.pdbqt'
                models = []
                in_model = False
                current_model = []
                model_num = 1
                
                for line in dlg_content.split('\\n'):
                    if 'DOCKED: MODEL' in line:
                        in_model = True
                        current_model = [f'MODEL {model_num}\\n']
                        model_num += 1
                    elif 'DOCKED: ENDMDL' in line:
                        current_model.append('ENDMDL\\n')
                        models.append(''.join(current_model))
                        in_model = False
                    elif in_model and line.startswith('DOCKED:'):
                        # Remove 'DOCKED: ' prefix
                        pdbqt_line = line[8:]
                        current_model.append(pdbqt_line + '\\n')
                
                if models:
                    with open(output_pdbqt, 'w') as f:
                        f.writelines(models)
                    print(f'\\nDocked poses saved to: {output_pdbqt.name}')
                    state.docking_results = output_pdbqt
                else:
                    print('\\nNo docked poses found in output.')
                    state.docking_results = dlg_path
            
            print('\\nDocking complete!')

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
    widgets.HTML('<h3>Docking Configuration</h3>'),
    docking_engine,
    HBox([vina_params, ad4_params]),
    run_docking_btn,
    docking_output
])

display(docking_panel)'''

                # Convert to list of lines for notebook format
                cell['source'] = [line + '\n' for line in new_source.split('\n')]
                # Remove trailing newline from last element
                if cell['source'] and cell['source'][-1].endswith('\n'):
                    cell['source'][-1] = cell['source'][-1].rstrip('\n')
                
                print("Replaced AutoDock4 placeholder with full implementation")
                break
    
    # Write back
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=1)
    
    print("[OK] AutoDock4 workflow has been fully implemented!")

if __name__ == "__main__":
    notebook_file = Path(__file__).parent / "autodock_workbench.ipynb"
    add_autodock4_workflow(notebook_file)
