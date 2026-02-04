# ==================== FIXED DOCKING ====================
# 
# This cell fixes the Vina --log option issue for Vina 1.2.x
# Replace the docking cell in your notebook with this code

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

def run_docking(btn):
    with docking_output:
        clear_output()
        
        if not state.receptor_pdbqt or not state.ligand_pdbqt:
            print('Please prepare both receptor and ligand first!')
            return
        
        engine = docking_engine.value
        print(f'Running {engine}...')
        
        if engine == 'AutoDock Vina':
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
            
            # Run Vina (FIXED: removed --log option for Vina 1.2.x)
            cmd = f'"{state.vina_exe}" --config "{config_path}"'
            print(f'Executing: {cmd}')
            print('This may take a few minutes...')
            
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=str(state.work_dir))
            
            # Vina outputs results to stdout
            if result.stdout:
                print('\n' + '='*50)
                print('VINA OUTPUT:')
                print('='*50)
                print(result.stdout)
            
            if result.returncode == 0 and output_path.exists():
                print('Docking complete!')
                
                # Count poses
                with open(output_path, 'r') as f:
                    num_poses = f.read().count('MODEL')
                print(f'Generated {num_poses} poses')
                
                state.docking_results = output_path
            else:
                if result.stderr:
                    print(f'Error: {result.stderr}')
        
        else:  # AutoDock4
            print('AutoDock4 workflow requires grid preparation.')
            print('Generating parameter files...')
            
            gpf_path = state.work_dir / f'{state.receptor_pdbqt.stem}.gpf'
            dpf_path = state.work_dir / f'{state.ligand_pdbqt.stem}.dpf'
            
            print(f'GPF will be saved to: {gpf_path}')
            print(f'DPF will be saved to: {dpf_path}')
            print('\nFull AD4 integration coming soon!')

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

display(docking_panel)
