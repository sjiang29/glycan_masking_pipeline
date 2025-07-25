
# 📊 Rosetta Score Plotting Utility

This utility generates customizable plots (box plots, violin plots, etc.) from Rosetta score files (e.g., `Glyc_score.sc`) across multiple glycosylation positions for antigens.

It uses a flexible configuration file so you can reuse the script easily for different antigens, thresholds, wild types, and plot types — without modifying the code.

---

## 🧪 Requirements

- Python 3.6+
- Libraries:
  - `matplotlib`
  - `numpy`

Install them using:

```bash
pip install matplotlib numpy
```

---

## 📁 Files

- `improved_plot_script.py` — Main plotting script
- `plot_config.json` — Configuration file specifying how to generate the plots

---

## ⚙️ How to Use

1. **Edit the Configuration File**

Modify the parameters in `plot_config.json`:

```json
{
  "working_dir": ".",                // Folder with pos_XXX subdirectories
  "score_file": "Glyc_score.sc",     // Score file name in each folder
  "top_n": 75,                       // Top N structures to use for plotting
  "start_pos": 1,                    // Starting AA position
  "end_pos": 268,                    // Ending AA position
  "threshold": -180,                 // Score threshold for filtering
  "wild_index": 200,                 // WT position for comparison
  "plot_type": "box_grid",           // Supported: box_grid
  "n_splits": 1,                     // Number of rows in the plot grid
  "title": "Influenza Top Scores"   // Plot title
}
```

2. **Run the Script**

```bash
python improved_plot_script.py
```

The script reads `plot_config.json`, processes the score files in the specified directory, and shows the generated plot interactively.

---

## 📦 Output

- Interactive matplotlib plots (shown in window)
- Can be saved manually or extended to auto-save to PNG/PDF

---

## ✏️ Example Use Cases

- Compare Rosetta scores across mutant vs. wild-type antigen positions
- Visualize structure selection quality per position
- Filter glycosylation positions based on energy thresholds

---

## 🧠 Tip

To switch to a new antigen or experiment:
- Just update the `working_dir`, `wild_index`, or `title` in `plot_config.json`
- No need to change the Python code!

---

## 📫 Contact

For improvements, issues, or feature requests, open an issue or contact the author.
