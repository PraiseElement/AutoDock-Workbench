# 🧬 AutoDock Workbench

**Interactive Molecular Docking Interface for AutoDock4 and AutoDock Vina**

A complete, user-friendly Jupyter Notebook workbench for performing molecular docking studies. This interface provides an interactive GUI for receptor preparation, ligand preparation, docking configuration, and results visualization.

---

## 📋 Table of Contents

- [Features](#-features)
- [Requirements](#-requirements)
- [Installation](#-installation)
- [Quick Start Guide](#-quick-start-guide)
- [Detailed User Guide](#-detailed-user-guide)
- [Offline Usage](#-offline-usage)
- [Uploading to GitHub](#-uploading-to-github)
- [Troubleshooting](#-troubleshooting)
- [Credits](#-credits)
- [Contact](#-contact)

---

## ✨ Features

- **🧬 Receptor Preparation**: Upload PDB files, clean structures, remove waters/heteroatoms, and prepare PDBQT files
- **💊 Ligand Preparation**: Upload SDF/MOL/PDB ligands, add hydrogens, compute Gasteiger charges, and prepare for docking
- **📦 Grid Box Configuration**: Visual grid box setup with interactive 3D visualization
- **⚙️ Docking Engines**: Support for both **AutoDock Vina** and **AutoDock4**
- **📊 Results Visualization**: 3D molecular viewer with docking poses and binding affinity analysis
- **🎨 Interactive UI**: User-friendly widgets for all operations

---

## 💻 Requirements

### System Requirements

- Windows 10/11 (64-bit), Linux, or macOS
- Python 3.8 or higher
- Jupyter Notebook or JupyterLab
- At least 4GB RAM recommended

### Python Dependencies

```
ipywidgets>=8.0.0
py3Dmol>=2.0.0
rdkit>=2022.09.1
numpy>=1.21.0
pandas>=1.3.0
```

### Docking Executables (Included)

The following executables are included in this package:

- `autodock4.exe` - AutoDock4 docking engine
- `autogrid4.exe` - AutoGrid4 for grid map generation
- `vina.exe` - AutoDock Vina
- `vina_split.exe` - Vina output splitter

---

## 🔧 Installation

### Step 1: Create a Conda Environment (Recommended)

```bash
# Create new environment
conda create -n autodock_workbench python=3.10 -y

# Activate environment
conda activate autodock_workbench
```

### Step 2: Install Dependencies

```bash
# Install Jupyter
pip install jupyter jupyterlab ipywidgets

# Install visualization tools
pip install py3Dmol

# Install RDKit (via conda for best compatibility)
conda install -c conda-forge rdkit -y

# Install other dependencies
pip install numpy pandas
```

### Step 3: Install AutoDockTools_py3

```bash
# Install from GitHub
pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3
```

**OR** install locally from this repository:

```bash
# Navigate to this directory and install
pip install -e .
```

### Step 4: Launch the Notebook

```bash
# Start Jupyter
jupyter notebook autodock_workbench.ipynb
```

---

## 🚀 Quick Start Guide

1. **Launch the notebook** in Jupyter
2. **Run Cell 1** (Setup & Imports) - Initializes all libraries
3. **Run Cell 2** (Global State) - Sets up working directories
4. **Upload your receptor** (PDB file) and click "Prepare Receptor"
5. **Upload your ligand** (SDF/MOL file) and click "Prepare Ligand"
6. **Configure the grid box** around your binding site
7. **Select docking engine** (Vina or AutoDock4)
8. **Run docking** and visualize results

---

## 📖 Detailed User Guide

### Cell 1: Setup & Imports

This cell loads all required Python libraries:

- **RDKit**: For molecular manipulation and 2D structure handling
- **py3Dmol**: For interactive 3D molecular visualization
- **ipywidgets**: For GUI elements
- **AutoDockTools**: For receptor/ligand preparation

**Expected Output**: `✅ All imports successful!`

---

### Cell 2: Global State

Creates a `WorkbenchState` object that tracks:

- Working directory paths
- Receptor and ligand file paths
- Grid box parameters
- Docking results

**Expected Output**: Shows paths to docking executables

---

### Cell 3: Receptor Preparation

**Steps:**

1. Click "Upload" button and select your PDB file
2. Select cleanup options:
   - `nphs` - Merge non-polar hydrogens
   - `lps` - Merge lone pairs
   - `waters` - Remove water molecules
   - `nonstdres` - Remove non-standard residues
3. Click **"Prepare Receptor"**
4. Wait for processing (may take 1-5 minutes for large structures)

**Output Files:**

- `receptor_prepared.pdbqt` in `workbench_files/`

**Tips:**

- For multi-chain proteins, you may need to select a specific chain
- Ensure the PDB has proper atom naming
- Remove any ligands already in the structure if doing blind docking

---

### Cell 4: Ligand Preparation

**Supported Formats:** SDF, MOL, MOL2, PDB

**Steps:**

1. Click "Upload" and select your ligand file
2. Review the 2D structure preview
3. Click **"Prepare Ligand"**

**What happens:**

- Adds hydrogens if missing
- Computes Gasteiger charges
- Identifies rotatable bonds
- Generates PDBQT format

**Output Files:**

- `ligand_prepared.pdbqt` in `workbench_files/`

---

### Cell 5: Grid Box Configuration

**Steps:**

1. Set center coordinates (X, Y, Z) for the binding site
2. Adjust box dimensions (default: 40×40×40 Å)
3. Set grid spacing (default: 0.375 Å for AutoDock4)
4. Visualize the grid box in 3D

**Tips:**

- Center the box on known binding site residues
- Ensure the box is large enough to encompass the entire binding pocket
- For Vina, larger boxes increase search time

---

### Cell 6: Docking Configuration & Execution

**Engine Options:**

| Parameter      | Vina | AutoDock4 |
| -------------- | ---- | --------- |
| Exhaustiveness | 8-32 | N/A       |
| Num Modes      | 1-20 | N/A       |
| GA Runs        | N/A  | 10-100    |
| Speed          | Fast | Slower    |
| Accuracy       | Good | Excellent |

**Steps:**

1. Select docking engine (Vina recommended for speed)
2. Adjust parameters as needed
3. Click **"Run Docking"**
4. Wait for completion

---

### Cell 7: Results Visualization

**Features:**

- 3D visualization of docking poses
- Binding affinity scores (kcal/mol)
- Pose ranking and selection
- Export options

---

## 📴 Offline Usage

**Yes, this workbench works completely offline!**

Once you have:

1. ✅ Installed all Python dependencies
2. ✅ Downloaded the repository with executables
3. ✅ Prepared your input PDB/SDF files

No internet connection is required for:

- Receptor preparation
- Ligand preparation
- Docking calculations
- Results visualization

The only features that may require internet:

- Installing packages for the first time
- Fetching PDB structures from online databases (not implemented in this version)

---

## 📤 Uploading to GitHub

### Step 1: Create a GitHub Account

If you don't have one, go to [github.com](https://github.com) and sign up.

### Step 2: Install Git

Download and install Git from [git-scm.com](https://git-scm.com/downloads)

### Step 3: Initialize Your Repository

```bash
# Navigate to this folder
cd "C:\Users\Administrator\Downloads\Compressed\AutoDockTools-master\AutoDockTools-master"

# Initialize git repository
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: AutoDock Workbench"
```

### Step 4: Create Repository on GitHub

1. Go to [github.com/new](https://github.com/new)
2. Name your repository (e.g., `AutoDock-Workbench`)
3. Add a description
4. Keep it **Public** or select **Private**
5. **Don't** initialize with README (you already have one)
6. Click **"Create repository"**

### Step 5: Push to GitHub

Replace `YOUR_USERNAME` with your GitHub username:

```bash
# Add remote origin
git remote add origin https://github.com/YOUR_USERNAME/AutoDock-Workbench.git

# Push to GitHub
git branch -M main
git push -u origin main
```

### Step 6: Verify Upload

Visit `https://github.com/YOUR_USERNAME/AutoDock-Workbench` to see your repository!

### ⚠️ Important Notes:

- The `.exe` files are large (~3.8 MB total). GitHub allows files up to 100 MB.
- Consider using [Git LFS](https://git-lfs.github.com/) for large files if needed.
- Make sure to update `.gitignore` to exclude any sensitive or temporary files.

---

## 🔧 Troubleshooting

### Common Issues

| Problem                      | Solution                                                                               |
| ---------------------------- | -------------------------------------------------------------------------------------- |
| `ModuleNotFoundError: rdkit` | Install via: `conda install -c conda-forge rdkit`                                      |
| `ImportError: AutoDockTools` | Install via: `pip install git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3` |
| Widgets not displaying       | Run: `jupyter nbextension enable --py widgetsnbextension`                              |
| Docking fails                | Check that all `.exe` files are present in the main directory                          |
| Receptor preparation slow    | Normal for large proteins (>5000 atoms). Be patient.                                   |

### Getting Help

If you encounter issues, please:

1. Check the error message carefully
2. Ensure all dependencies are installed
3. Contact us using the information below

---

## 👏 Credits

### Original AutoDockTools_py3 Authors

This workbench is built upon the excellent **AutoDockTools_py3** library:

- **MS. Valdes-Tresanco** - [GitHub](https://github.com/Valdes-Tresanco-MS)
- **ME. Valdes-Tresanco**

Original repository: [github.com/Valdes-Tresanco-MS/AutoDockTools_py3](https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3)

### Interactive Workbench Development

The Jupyter Notebook interface and workbench enhancements were developed by:

**Praize Kene Ugwu**

- 📧 Email: praizekene1@gmail.com
- 📱 Phone/WhatsApp: +2348144680669

---

## 📄 License

This project uses components from AutoDockTools_py3 which is distributed under its original license. Please refer to the original repository for licensing details.

---

## 🤝 Contributing

Contributions are welcome! Feel free to:

- Report bugs
- Suggest new features
- Submit pull requests

---

_Last updated: February 2026_
