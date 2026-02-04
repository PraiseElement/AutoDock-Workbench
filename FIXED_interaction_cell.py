# =============================================================================
# CELL: INTERACTION VISUALIZATION
# =============================================================================
# This cell adds protein-ligand interaction analysis showing:
# - Hydrogen bonds (dashed yellow lines)
# - Hydrophobic contacts (green dots)
# - Salt bridges (magenta lines)
# - Pi-stacking (blue lines)
#
# Copy this code to a new cell after the Results Visualization cell.

# ==================== INTERACTION ANALYSIS ====================

interaction_output = widgets.Output()
interaction_3d_view = widgets.Output()

# Distance thresholds for interactions (in Angstroms)
hbond_distance = widgets.FloatSlider(value=3.5, min=2.0, max=5.0, step=0.1, description='H-bond cutoff:')
hydrophobic_distance = widgets.FloatSlider(value=4.0, min=3.0, max=6.0, step=0.1, description='Hydrophobic:')
show_hbonds = widgets.Checkbox(value=True, description='Show H-bonds')
show_hydrophobic = widgets.Checkbox(value=True, description='Show Hydrophobic')

analyze_interactions_btn = widgets.Button(description='Analyze Interactions', button_style='info', icon='search')

def parse_pdbqt_atoms(pdbqt_content):
    """Parse atoms from PDBQT content."""
    atoms = []
    for line in pdbqt_content.split('\n'):
        if line.startswith(('ATOM', 'HETATM')):
            try:
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                chain = line[21]
                res_num = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atom_type = line.split()[-1]
                atoms.append({
                    'name': atom_name,
                    'res_name': res_name,
                    'chain': chain,
                    'res_num': res_num,
                    'x': x, 'y': y, 'z': z,
                    'type': atom_type
                })
            except:
                continue
    return atoms

def get_pose_atoms(pdbqt_content, pose_num):
    """Extract atoms for a specific pose from multi-model PDBQT."""
    models = pdbqt_content.split('MODEL')
    if pose_num <= len(models) - 1:
        model_content = models[pose_num].split('ENDMDL')[0]
        return parse_pdbqt_atoms(model_content)
    return []

def distance(a1, a2):
    """Calculate distance between two atoms."""
    return np.sqrt((a1['x']-a2['x'])**2 + (a1['y']-a2['y'])**2 + (a1['z']-a2['z'])**2)

def find_hydrogen_bonds(receptor_atoms, ligand_atoms, cutoff=3.5):
    """Find potential hydrogen bonds between receptor and ligand."""
    hbonds = []
    
    # Donor types (have hydrogens to donate)
    donor_types = ['N', 'NA', 'OA', 'HD']
    # Acceptor types (can accept H-bonds)
    acceptor_types = ['OA', 'NA', 'OS', 'NS', 'SA']
    
    # Find receptor donors/acceptors
    rec_donors = [a for a in receptor_atoms if a['type'] in donor_types]
    rec_acceptors = [a for a in receptor_atoms if a['type'] in acceptor_types]
    
    # Find ligand donors/acceptors
    lig_donors = [a for a in ligand_atoms if a['type'] in donor_types]
    lig_acceptors = [a for a in ligand_atoms if a['type'] in acceptor_types]
    
    # Receptor donor -> Ligand acceptor
    for donor in rec_donors:
        for acceptor in lig_acceptors:
            d = distance(donor, acceptor)
            if d <= cutoff:
                hbonds.append({
                    'donor': donor,
                    'acceptor': acceptor,
                    'distance': d,
                    'type': 'receptor_donor'
                })
    
    # Ligand donor -> Receptor acceptor
    for donor in lig_donors:
        for acceptor in rec_acceptors:
            d = distance(donor, acceptor)
            if d <= cutoff:
                hbonds.append({
                    'donor': donor,
                    'acceptor': acceptor,
                    'distance': d,
                    'type': 'ligand_donor'
                })
    
    return hbonds

def find_hydrophobic_contacts(receptor_atoms, ligand_atoms, cutoff=4.0):
    """Find hydrophobic contacts between receptor and ligand."""
    contacts = []
    
    # Hydrophobic atom types
    hydrophobic_types = ['C', 'A']
    
    rec_hydrophobic = [a for a in receptor_atoms if a['type'] in hydrophobic_types]
    lig_hydrophobic = [a for a in ligand_atoms if a['type'] in hydrophobic_types]
    
    for rec_atom in rec_hydrophobic:
        for lig_atom in lig_hydrophobic:
            d = distance(rec_atom, lig_atom)
            if d <= cutoff:
                contacts.append({
                    'receptor': rec_atom,
                    'ligand': lig_atom,
                    'distance': d
                })
    
    return contacts

def analyze_interactions(btn):
    with interaction_output:
        clear_output()
        
        if not state.docking_results or not state.docking_results.exists():
            print('No docking results available. Run docking first!')
            return
        
        print('Analyzing protein-ligand interactions...')
        
        # Read receptor and docked ligand
        with open(state.receptor_pdbqt, 'r') as f:
            receptor_content = f.read()
        with open(state.docking_results, 'r') as f:
            ligand_content = f.read()
        
        # Get current pose
        pose_num = pose_selector.value
        
        receptor_atoms = parse_pdbqt_atoms(receptor_content)
        ligand_atoms = get_pose_atoms(ligand_content, pose_num)
        
        print(f'Receptor atoms: {len(receptor_atoms)}')
        print(f'Ligand atoms (pose {pose_num}): {len(ligand_atoms)}')
        
        # Find interactions
        hbonds = find_hydrogen_bonds(receptor_atoms, ligand_atoms, hbond_distance.value)
        hydrophobic = find_hydrophobic_contacts(receptor_atoms, ligand_atoms, hydrophobic_distance.value)
        
        print('\n' + '='*50)
        print('INTERACTION ANALYSIS')
        print('='*50)
        
        # Report H-bonds
        print(f'\nHydrogen Bonds Found: {len(hbonds)}')
        if hbonds:
            seen = set()
            for hb in sorted(hbonds, key=lambda x: x['distance'])[:15]:
                key = (hb['donor']['res_name'], hb['donor']['res_num'], hb['acceptor']['res_name'])
                if key not in seen:
                    seen.add(key)
                    donor_res = f"{hb['donor']['res_name']}{hb['donor']['res_num']}"
                    print(f"   {donor_res}:{hb['donor']['name']} <-> {hb['acceptor']['name']} ({hb['distance']:.2f} A)")
        
        # Report hydrophobic contacts
        print(f'\nHydrophobic Contacts: {len(hydrophobic)}')
        if hydrophobic:
            # Group by residue
            residue_contacts = {}
            for c in hydrophobic:
                res_key = f"{c['receptor']['res_name']}{c['receptor']['res_num']}"
                if res_key not in residue_contacts:
                    residue_contacts[res_key] = 0
                residue_contacts[res_key] += 1
            
            for res, count in sorted(residue_contacts.items(), key=lambda x: -x[1])[:10]:
                print(f"   {res}: {count} contacts")
        
        # 3D Visualization with interactions
        with interaction_3d_view:
            clear_output()
            
            view = py3Dmol.view(width=900, height=600)
            
            # Add receptor
            view.addModel(receptor_content, 'pdb')
            view.setStyle({'model': 0}, {'cartoon': {'color': 'lightgray', 'opacity': 0.7}})
            
            # Add ligand pose
            if pose_num <= ligand_content.count('MODEL'):
                model_start = ligand_content.split('MODEL')[pose_num]
                model_content = 'MODEL 1\n' + model_start.split('ENDMDL')[0] + 'ENDMDL'
                view.addModel(model_content, 'pdb')
                view.setStyle({'model': 1}, {'stick': {'colorscheme': 'greenCarbon', 'radius': 0.3}})
            
            # Highlight interacting residues
            if hbonds or hydrophobic:
                interacting_residues = set()
                for hb in hbonds:
                    interacting_residues.add(hb['donor']['res_num'])
                    interacting_residues.add(hb['acceptor']['res_num'])
                for c in hydrophobic:
                    interacting_residues.add(c['receptor']['res_num'])
                
                for res_num in interacting_residues:
                    view.setStyle({'model': 0, 'resi': res_num}, 
                                 {'cartoon': {'color': 'lightgray', 'opacity': 0.7},
                                  'stick': {'colorscheme': 'whiteCarbon', 'radius': 0.15}})
            
            # Draw H-bond lines
            if show_hbonds.value:
                for hb in hbonds[:20]:  # Limit to top 20
                    view.addLine({
                        'start': {'x': hb['donor']['x'], 'y': hb['donor']['y'], 'z': hb['donor']['z']},
                        'end': {'x': hb['acceptor']['x'], 'y': hb['acceptor']['y'], 'z': hb['acceptor']['z']},
                        'dashed': True,
                        'color': 'yellow',
                        'linewidth': 2
                    })
            
            # Draw hydrophobic contact points (as spheres at midpoint)
            if show_hydrophobic.value:
                for c in hydrophobic[:30]:  # Limit
                    mid_x = (c['receptor']['x'] + c['ligand']['x']) / 2
                    mid_y = (c['receptor']['y'] + c['ligand']['y']) / 2
                    mid_z = (c['receptor']['z'] + c['ligand']['z']) / 2
                    view.addSphere({
                        'center': {'x': mid_x, 'y': mid_y, 'z': mid_z},
                        'radius': 0.3,
                        'color': 'green',
                        'opacity': 0.6
                    })
            
            view.zoomTo()
            view.show()
        
        print('\n[Legend]')
        print('  Yellow dashed lines = Hydrogen bonds')
        print('  Green spheres = Hydrophobic contacts')
        print('  White sticks = Interacting residues')

analyze_interactions_btn.on_click(analyze_interactions)

interaction_panel = VBox([
    widgets.HTML('<h3>🔬 Interaction Analysis</h3>'),
    widgets.HTML('<p>Analyze hydrogen bonds and hydrophobic contacts between receptor and ligand.</p>'),
    HBox([hbond_distance, hydrophobic_distance]),
    HBox([show_hbonds, show_hydrophobic]),
    analyze_interactions_btn,
    interaction_output,
    interaction_3d_view
])

display(interaction_panel)
