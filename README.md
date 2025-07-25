# 🧬 Glycan-Masking Rosetta Analysis Pipeline

This is a modular pipeline for simulating, analyzing, and visualizing glycan-masked antigens and their binding profiles with antibodies using Rosetta and Python-based tools. The workflow consists of **four sequential modules**:

---

## 🧭 Pipeline Overview

1. **[Glycan Masking Simulation (Rosetta)]**
2. **[Glycan Score Plotting]**
3. **[Antibody Binding Simulation (Rosetta)]**
4. **[Antibody Binding Score Plotting]**

---

## 🔁 Step-by-Step Usage

### 1️⃣ Glycan Masking Simulation
📄 Scripts: `single_glycan.py`, `local_parallel_glycan_run.py`

Simulates glycan attachment at N-linked sites using Rosetta `FastDesign` and `GlycanTreeModeler`.

- Inputs: FASTA, PDB, glycosylation positions
- Outputs: glycosylated PDB structures + score files
- Parallelized across positions

📘 Detailed guide: [glycan_masking_readme](https://github.com/sjiang29/glycan_masking_pipeline/blob/main/glycan_masking/README_local_parallel_glycan_run.md)

---

### 2️⃣ Glycan Score Plotting
📄 Script: `improved_plot_script.py`  
⚙️ Config: `plot_config.json`

Plots Rosetta energy scores (REU) of glycosylated variants.

- Box plots or violin plots per position
- Supports filtering by score threshold
- Highlight wild-type baseline

📘 Detailed guide: [README_plotting.md](README_plotting.md)

---

### 3️⃣ Antibody Binding Simulation
📄 Scripts: `single_ab_binding.py`, `local_parallel_pos_binding_run.py`

Performs antibody chain addition and binding energy calculations using Rosetta on selected glycosylated antigens.

- Inputs: top-scoring glycosylated models from step 2
- Outputs: antibody-bound complex PDBs + binding scores
- Parallelized across positions and antibodies

📘 Detailed guide: [README_binding.md](README_binding.md)

---

### 4️⃣ Antibody Binding Score Plotting
📄 Script: `refactored_ab_plot_script.py`  
⚙️ Config: `ab_plot_config.json`

Plots antibody binding scores across antigen positions using violin plots.

- One violin plot per antibody
- Highlights wild-type (position 200) as baseline
- Configurable columns (e.g., `dG_separated`, `total_score`)

📘 Detailed guide: [README_ab_plotting.md](README_ab_plotting.md)

---

## 🗂️ Organizational Note

Although the pipeline is conceptually organized into four stages (masking → glycan plotting → binding → binding plotting), **all Python scripts and config files are expected to be placed in a single flat directory** to ensure compatibility with relative file paths in the code.

You can still organize input and output data (e.g., `template_files/`, `antigens/`, `binding_results/`) in subdirectories.

---

## 🔄 Suggested Workflow

```text
1. Run glycan simulations across positions          👉 Step 1
2. Analyze scores and visualize top variants        👉 Step 2
3. Perform antibody binding simulations on hits     👉 Step 3
4. Compare antibody binding profiles across sites   👉 Step 4
```

---

## 🧠 Tip

To repeat the process for a new antigen:
- Update your input FASTA/PDB files
- Adjust `antigens.csv` and `plot_config.json`
- Run all four modules in sequence

---

## 📫 Support

For help, issues, or contributions, feel free to reach out or open a GitHub issue.
