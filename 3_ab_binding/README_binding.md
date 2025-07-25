
# 🤝 Antibody-Antigen Binding Energy Analysis Pipeline

This pipeline automates the process of analyzing antibody-antigen binding by:
1. Adding antibody chains to antigen structures
2. Running Rosetta FastRelax and docking
3. Calculating binding energy scores
4. Selecting top-ranked antigen structures

It supports parallel execution across multiple glycosylation positions and antibodies.

---

## 📁 Project Structure

```
.
├── binding_template_files/       # XML, options, and shell templates for binding
├── Abs_all.csv                   # Metadata for antibodies
├── Glyc_score.sc                 # Score file for all glycosylated structures
├── single_ab_binding.py          # Class to run antibody chain addition and binding
├── local_parallel_pos_binding_run.py # Driver for parallel binding energy calculations
```

---

## 🧪 Prerequisites

- Python 3.6+
- Rosetta suite with FastRelax and scoring apps
- Shell scripts:
  - `addChain.sh`
  - `cal_binding_energy.sh`
  - `clean.sh`

- Template files:
  - `template_AddChain_FastRelax_Multi.xml`
  - `template_docking_analysis.xml`

---

## 📄 Input: `Abs_all.csv`

A CSV with the following format:

```csv
antibody,antibody_heavy_chain_pdb,antibody_light_chain_pdb
Ab1,Ab1_H.pdb,Ab1_L.pdb
Ab2,Ab2_H.pdb,Ab2_L.pdb
...
```

---

## ⚙️ How to Run

```bash
python local_parallel_pos_binding_run.py <n_cores> <positions> <top_n> <abs_csv>
```

### Arguments:

- `n_cores`: Number of parallel processes to use
- `positions`: Positions as range or comma-separated list (e.g., `1-10,15,20`)
- `top_n`: Top N glycosylated structures (by score) to analyze per position
- `abs_csv`: Path to antibody CSV file

### Example:

```bash
python local_parallel_pos_binding_run.py 4 1-5,10 50 Abs_all.csv
```

---

## 🧬 Workflow Summary

For each glycosylation site (pos_X):

1. **Select Top N Models**:
   - Based on total score in `Glyc_score.sc`.

2. **Add Antibody Chains**:
   - Use `AddChain` XML and `addChain.sh` to attach antibody H/L chains.

3. **Run Rosetta Relaxation**:
   - FastRelax is run on each antigen-antibody complex.

4. **Clean and Prepare for Docking**:
   - Files are cleaned using `clean.sh`.

5. **Binding Energy Calculation**:
   - Docking analysis via `cal_binding_energy.sh` and `template_docking_analysis.xml`.

---

## 📦 Output

- Cleaned and relaxed antigen-antibody complexes named:
  ```
  Add_Rlx_<Ab><Antigen>.pdb
  ```

- Binding scores printed to output logs

---

## 🧵 Parallel Execution

- One process per position using Python `multiprocessing`
- Antibodies are looped within each process
- Each process:
  - Copies binding templates
  - Runs chain addition and binding score calculation

---

## 📌 Notes

- Make sure your Rosetta applications and scripts are executable (`chmod +x *.sh`)
- Template files and antibodies must be provided before running the pipeline
- Input antigen structures should already be generated and scored

---

## 📫 Contact

For issues, bug reports, or collaboration inquiries, please open an issue or contact the maintainer.
