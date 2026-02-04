# ==================== FIXED RECEPTOR PREPARATION ====================
# 
# This cell fixes issues with:
#   1. Duplicate atoms at same coordinates (causes ZeroDivisionError)
#   2. 3D preview not displaying
#   3. Better error handling

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
        ('No repairs (skip hydrogen addition)', '')
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

def clean_pdb_duplicates(pdb_path):
    """Remove duplicate atoms with identical coordinates from PDB file"""
    seen_coords = {}
    cleaned_lines = []
    duplicates_removed = 0
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                # Extract coordinates (columns 31-54 in PDB format)
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coord_key = (round(x, 3), round(y, 3), round(z, 3))
                    
                    if coord_key in seen_coords:
                        duplicates_removed += 1
                        continue  # Skip duplicate
                    seen_coords[coord_key] = True
                except:
                    pass
            cleaned_lines.append(line)
    
    # Write cleaned file
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
        
        # Clean duplicate atoms
        receptor_path, dups = clean_pdb_duplicates(receptor_path)
        if dups > 0:
            print(f'Removed {dups} duplicate atoms')
            state.receptor_pdb = receptor_path
        
        # Prepare receptor
        try:
            mol = Read(str(receptor_path))[0]
            cleanup = '_'.join(receptor_cleanup.value) if receptor_cleanup.value else ''
            repairs = receptor_repairs.value
            
            output_pdbqt = state.work_dir / f'{receptor_path.stem}_prepared.pdbqt'
            
            print(f'Atoms: {len(mol.allAtoms)}')
            if len(mol.allAtoms) > 3000:
                print('Note: Large protein - this may take a moment...')
            
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
            print(f'   Residues: {len(mol.chains.residues)}')
            
            # Calculate center of mass for grid
            coords = mol.allAtoms.coords
            center = np.mean(coords, axis=0)
            state.grid_center = center.tolist()
            print(f'   Center: ({center[0]:.2f}, {center[1]:.2f}, {center[2]:.2f})')
            
        except Exception as e:
            error_msg = str(e)
            print(f'Error: {error_msg}')
            
            # Suggest workaround for common issues
            if 'ZeroDivision' in error_msg or 'same coordinates' in error_msg:
                print('\nTry selecting "No repairs" to skip hydrogen addition,')
                print('or clean the PDB file to remove duplicate atoms.')
            elif 'babel_type' in error_msg:
                print('\nThe structure may have unusual atom types.')
                print('Try selecting "No repairs" option.')
            
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

# Warning for problematic PDBs
warning_text = widgets.HTML('''
<div style="background:#4a3d2d; padding:10px; border-radius:5px; margin:5px 0;">
<b>Tips:</b><br>
- For PDBs with errors, try "No repairs" option<br>
- Large proteins (>3000 atoms) may take a while<br>
- Consider using a single chain for multimeric proteins
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

display(receptor_panel)
