
# 🧬 Glycan Attachment Pipeline using Rosetta

This project automates the glycan attachment process to desired asparagine (N) residues on antigen proteins using Rosetta's `FastDesign` and `GlycanTreeModeler`. The pipeline supports batch processing of multiple antigens and parallel execution across multiple CPU cores.

---

## 📁 Project Structure

```
.
├── glycan_template_files/        # Contains example XMLs, options, and shell scripts
├── antigens/                     # Contains input PDB and FASTA files for antigens
├── antigens.csv                  # Metadata for antigen inputs
├── single_glycan.py              # Main class to run glycan addition on a single site
├── local_parallel_glycan_run.py # Parallel processing script for batch execution
├── FastDesign.sh                 # Shell wrapper for Rosetta FastDesign
├── GlycanTreeModeler.sh         # Shell wrapper for Rosetta GlycanTreeModeler
```

---

## 🧪 Prerequisites

- Python 3.6+
- Rosetta suite installed and `FastDesign` and `GlycanTreeModeler` available via command-line
- Linux or macOS environment
- Required Rosetta `.xml` and `.options` templates in `glycan_template_files/`
- PDB and FASTA input files for antigens in `antigens/`

---

## 📄 Input: `antigens.csv`

A CSV file specifying antigen inputs with the following header:

```csv
antigen,pdb,fasta,native_glycan_pos
Spike,spike.pdb,spike.fasta,154+234+678
...
```

- `antigen`: folder name for the antigen
- `pdb`: structure file in PDB format
- `fasta`: sequence file in FASTA format
- `native_glycan_pos`: positions (1-based) of native glycans separated by `+`, or empty if none

---

## ⚙️ How to Run

```bash
python local_parallel_glycan_run.py <n_cores> <positions> <n_structs> <antigens.csv>
```

### Arguments:

- `n_cores`: Number of parallel processes to use
- `positions`: Comma-separated positions (1-based) to add glycans, e.g., `12,45,78`
- `n_structs`: Number of structures to generate per position
- `antigens.csv`: Path to the CSV input file

### Example:

```bash
python local_parallel_glycan_run.py 4 12,34 50 antigens.csv
```

---

## 🚀 Pipeline Steps

For each specified antigen and glycosylation site:

1. **Prepare working directory** per position.
2. **Create Rosetta resfile** to enforce N-X-S/T motif.
3. **Modify FastDesign XML** to include the resfile.
4. **Run FastDesign** to stabilize the new motif.
5. **Modify GlycanTreeModeler XML** to specify glycosylation positions.
6. **Run GlycanTreeModeler** to attach glycan trees.

---

## 🧠 Class: `Single_Glycan`

Core functionality is handled by `single_glycan.py`, which:
- Parses FASTA and infers target sequence
- Checks and mutates sequence to ensure glycosylation motif
- Generates Rosetta XML and resfile templates
- Calls Rosetta FastDesign and GlycanTreeModeler using subprocess

---

## 🛠 Template Files

Ensure the following files exist in `glycan_template_files/`:

- `template_FastDesign.xml`
- `FastDesign.options`
- `template_GlycanTreeModeler.xml`
- `GlycanTreeModeler.options`
- `FastDesign.sh`
- `GlycanTreeModeler.sh`

These are copied per run into antigen/position folders.

---

## 📦 Output

Outputs are saved in a folder hierarchy like:

```
Spike/
├── pos_12/
│   ├── Glyc_Des_spike_0001.pdb
│   └── ...
├── pos_34/
│   └── Glyc_Des_spike_0001.pdb
```

Each `Glyc_Des_*.pdb` contains the glycosylated model.

---

## 🧵 Parallelism Notes

- Uses Python's `multiprocessing` to distribute glycan addition jobs across CPU cores.
- Each glycosylation position is processed independently.
- Antigens are processed sequentially; sites within each antigen are parallelized.

---

## 📌 Known Assumptions

- Glycan is added only to N-linked motifs (N-X-S/T).
- The motif is checked and enforced via mutation in the resfile.
- Positions in input are **1-based**.

---

## 📫 Contact

For issues or contributions, feel free to open a GitHub issue or contact the maintainer.
