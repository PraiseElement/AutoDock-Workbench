#!/usr/bin/env python3
"""
Script to make Jupyter notebook code cells collapsible (hidden by default).
This adds the 'source_hidden' metadata to each code cell.
"""

import json
from pathlib import Path

def make_cells_collapsible(notebook_path):
    """Add source_hidden metadata to all code cells in a notebook."""
    
    # Read the notebook
    with open(notebook_path, 'r', encoding='utf-8') as f:
        notebook = json.load(f)
    
    # Count modified cells
    modified = 0
    
    # Iterate through cells
    for cell in notebook.get('cells', []):
        if cell.get('cell_type') == 'code':
            # Initialize metadata if not present
            if 'metadata' not in cell:
                cell['metadata'] = {}
            
            # Add jupyter source_hidden metadata
            if 'jupyter' not in cell['metadata']:
                cell['metadata']['jupyter'] = {}
            
            cell['metadata']['jupyter']['source_hidden'] = True
            modified += 1
    
    # Write back
    with open(notebook_path, 'w', encoding='utf-8') as f:
        json.dump(notebook, f, indent=1)
    
    print(f"[OK] Modified {modified} code cells to be collapsible (source hidden)")
    return modified

if __name__ == "__main__":
    notebook_file = Path(__file__).parent / "autodock_workbench.ipynb"
    make_cells_collapsible(notebook_file)
