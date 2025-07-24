
import os
import sys
import json
import matplotlib.pyplot as plt
import numpy as np
import statistics

def load_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

def safe_open_file(filename, mode='r'):
    try:
        return open(filename, mode)
    except FileNotFoundError:
        print(f"File not found: {filename}")
        return None
    except IOError as e:
        print(f"I/O error({e.errno}): {e.strerror}")
        return None

def get_column_values(sc_file, column_name):
    values = []
    f = safe_open_file(sc_file)
    if f is not None:
        with f:
            for line in f:
                if line.startswith("SCORE:") and 'total_score' in line:
                    header = line.strip().split()
                    break
            column_index = header.index(column_name)
            for line in f:
                if line.startswith("SCORE:"):
                    parts = line.strip().split()
                    try:
                        values.append(float(parts[column_index]))
                    except (IndexError, ValueError):
                        continue
    return values

def get_sc_file_path(sc_type, ab_name, pos_dir):
    if sc_type == "binding":
        return os.path.join(pos_dir, f"Add_Rlx_{ab_name}score.fasc")
    elif sc_type == "docking":
        return os.path.join(pos_dir, f"InAn_cleaned_{ab_name}.sc")
    else:
        raise ValueError("Unsupported score type")

def get_single_pos_score(pos, sc_type, column_name, ab_names, working_dir):
    pos_dir = os.path.join(working_dir, f"pos_{pos}")
    scores = {}
    for ab_name in ab_names:
        sc_file_path = get_sc_file_path(sc_type, ab_name, pos_dir)
        try:
            score = get_column_values(sc_file_path, column_name)
            if len(score) == 75:
                scores[ab_name] = score
        except Exception as e:
            print(f"Error reading {sc_file_path}: {e}")
    return scores

def get_score_for_ab(ab_name, sc_type, column_name, binding_poses, working_dir):
    scores = {}
    for pos in binding_poses:
        score_for_pos = get_single_pos_score(pos, sc_type, column_name, [ab_name], working_dir)
        if ab_name in score_for_pos:
            scores[pos] = score_for_pos[ab_name]
    return scores

def make_violin_plot_for_all(all_ab_scores_dict, sc_type, column_name, wild_index=200):
    n_abs = len(all_ab_scores_dict)
    fig, axes = plt.subplots(n_abs, 1, figsize=(16, 3.5 * n_abs), squeeze=False)

    for i, (ab_name, pos_and_scores) in enumerate(all_ab_scores_dict.items()):
        keys = list(pos_and_scores.keys())
        values = list(pos_and_scores.values())

        if wild_index in pos_and_scores:
            wild_idx = keys.index(wild_index)
            wild_value = values[wild_idx]
            del keys[wild_idx]
            del values[wild_idx]
            keys = [wild_index] + keys
            values = [wild_value] + values
            avg_value_of_wild = np.mean(wild_value)
        else:
            avg_value_of_wild = None

        pos = np.arange(1, len(keys) + 1)
        ax = axes[i, 0]
        ax.violinplot(values, positions=pos, showmeans=False, showmedians=True, widths=0.7)
        ax.set_xticks(pos)
        ax.set_xticklabels(keys, rotation=45, fontsize=8)
        ax.set_ylim(np.min(values), np.max(values))

        if avg_value_of_wild is not None:
            ax.axhline(avg_value_of_wild, color='r', linestyle='dotted', linewidth=1, label=f"WT mean: {avg_value_of_wild:.2f}")
            ax.legend(fontsize=8)

        if i < n_abs - 1:
            ax.set_xticklabels([])
            ax.set_xlabel('')

    fig.text(0.04, 0.5, column_name, va='center', rotation='vertical', fontsize=12)
    fig.text(0.5, 0.01, 'Amino Acid Position', ha='center', fontsize=12)
    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])
    plt.show()

def print_score_summary(ab_name, pos_and_score, column_name):
    print(f"----- {len(pos_and_score)} positions for {ab_name} -----")
    for k, v in pos_and_score.items():
        print(f"Pos {k}: mean = {statistics.mean(v):.2f}, median = {statistics.median(v):.2f}")

def main():
    config = load_config("ab_plot_config.json")
    sc_type = config["score_type"]
    column_name = config["column_name"]
    ab_names = config["antibodies"]
    binding_poses = config["binding_positions"]
    wild_index = config.get("wild_index", 200)
    working_dir = config.get("working_dir", os.getcwd())

    all_scores = {}
    for ab_name in ab_names:
        scores = get_score_for_ab(ab_name, sc_type, column_name, binding_poses, working_dir)
        print_score_summary(ab_name, scores, column_name)
        all_scores[ab_name] = scores

    make_violin_plot_for_all(all_scores, sc_type, column_name, wild_index=wild_index)

if __name__ == "__main__":
    main()
