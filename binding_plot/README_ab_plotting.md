
# 🧬 Antibody-Antigen Binding Score Plotting Utility

This utility allows you to visualize antibody-antigen interaction scores across glycosylation positions using violin plots. It's designed to work with Rosetta-generated `.fasc` or `.sc` files and supports multiple antibodies, score columns, and binding poses — all customizable via a configuration file.

---

## 📦 Files

- `refactored_ab_plot_script.py`: The main script for plotting scores.
- `ab_plot_config.json`: Configuration file to specify antibodies, positions, and score types.

---

## 🧪 Requirements

- Python 3.6+
- Libraries:
  - `matplotlib`
  - `numpy`

Install them with:

```bash
pip install matplotlib numpy
```

---

## ⚙️ Configuration: `ab_plot_config.json`

Example:

```json
{
  "score_type": "binding",              // "binding" or "docking"
  "column_name": "dG_separated",        // Column from .fasc/.sc file to extract
  "antibodies": ["H7-200", "m826", "FluA-20", "H7point5"],  // List of antibodies
  "binding_positions": [4, 5, 7, 8, 13, 16, 22, 34, 37, 40], // Positions to analyze
  "wild_index": 200,                   // Position for wild-type reference
  "working_dir": "."                   // Directory containing pos_XXX subfolders
}
```

You can modify this file to switch antibodies, adjust plot positions, or change the metric to analyze (e.g., `total_score`, `fa_atr`, `dG_cross`, etc.).

---

## 🚀 How to Run

```bash
python refactored_ab_plot_script.py
```

This will:
1. Parse the configuration from `ab_plot_config.json`
2. Extract scores for each antibody across defined positions
3. Generate violin plots with wild-type (position 200) as a reference
4. Print score summary (mean, median) per position

---

## 📈 Output

- Interactive violin plots (one per antibody)
- Each plot compares binding scores across antigen positions
- Red dotted line marks average score for wild-type (if available)

---

## 🧠 Use Cases

- Visual comparison of binding energy changes due to glycosylation
- Antibody performance analysis across antigen mutants
- Identification of escape or stabilizing mutations

---

## 🧰 Notes

- Each `pos_XXX` folder must contain one `.fasc` or `.sc` file per antibody
- Score files must include Rosetta `SCORE:` lines with the specified column
- Set `"score_type"` to `"binding"` or `"docking"` depending on the file format

---

## 📫 Contact

Feel free to suggest improvements, request support for new formats, or contribute enhancements!

