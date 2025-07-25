
import matplotlib.pyplot as plt
import os
import sys
import statistics
import numpy as np
import json

def load_config(config_path):
    with open(config_path, 'r') as f:
        return json.load(f)

def get_top_n(top_n, sc_file):
    with open(sc_file) as file:
        lines = file.readlines()[2:102]
    scores = []
    for line in lines:
        try:
            score = float(line.split(":")[1].split()[0])
            scores.append(score)
        except Exception:
            continue
    scores_100 = sorted(scores[:100])
    scores = sorted(scores)
    return scores[:top_n], scores_100[:50], scores

def get_top_n_for_all(working_dir, sc_file, top_n, start, end, wild_index):
    pos_and_scores, pos_and_scores_50, pos_and_avg_score = {}, {}, {}
    for filename in os.listdir(working_dir):
        f = os.path.join(working_dir, filename)
        if os.path.isdir(f) and "pos_" in f:
            score_file_path = os.path.join(f, sc_file)
            if os.path.exists(score_file_path):
                pos = int(f.split("_")[-1])
                top_n_scores, scores_50, scores = get_top_n(top_n, score_file_path)
                pos_and_avg_score[pos] = statistics.mean(scores)
                if start <= pos <= end or pos == wild_index:
                    pos_and_scores[pos] = top_n_scores
                    pos_and_scores_50[pos] = scores_50
    return dict(sorted(pos_and_scores.items())), dict(sorted(pos_and_scores_50.items())), pos_and_avg_score

def make_box_plot_grid(pos_and_scores, title, threshold, n_splits, wild_type=200):
    pos_and_scores = {k: v for k, v in pos_and_scores.items() if statistics.mean(v) < threshold}
    keys = list(pos_and_scores.keys())
    values = list(pos_and_scores.values())
    avg_wild = np.mean(pos_and_scores[wild_type]) if wild_type in pos_and_scores else None
    split_size = len(keys) // n_splits
    remainder = len(keys) % n_splits
    splits, start = [], 0
    for i in range(n_splits):
        end = start + split_size + (1 if i < remainder else 0)
        splits.append((keys[start:end], values[start:end]))
        start = end
    fig, axes = plt.subplots(n_splits, 1, figsize=(10, 3.5 * n_splits))
    if n_splits == 1:
        axes = [axes]
    for i, (ks, vs) in enumerate(splits):
        box = axes[i].boxplot(vs, labels=ks, patch_artist=True)
        for patch in box['boxes']:
            patch.set(facecolor='grey')
        axes[i].tick_params(labelsize=12)
        axes[i].set_ylim(np.min(values) - 1000, threshold)
        axes[i].set_ylabel('Rosetta Score', fontsize=14)
        if avg_wild:
            axes[i].axhline(avg_wild, color='r', linestyle='dotted', label=f"WT@{wild_type}: {avg_wild:.2f}")
            axes[i].legend(fontsize=8)
    plt.suptitle(title, fontsize=12)
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.show()

def main():
    config = load_config("glycan_plot_config.json")
    working_dir = config["working_dir"]
    sc_file = config.get("score_file", "Glyc_score.sc")
    top_n = config.get("top_n", 75)
    start = config.get("start_pos", 1)
    end = config.get("end_pos", 268)
    threshold = config.get("threshold", -180)
    wild_index = config.get("wild_index", 200)
    plot_type = config.get("plot_type", "box_grid")
    n_splits = config.get("n_splits", 1)
    plot_title = config.get("title", "Antigen")

    pos_scores, pos_scores_50, _ = get_top_n_for_all(working_dir, sc_file, top_n, start, end, wild_index)

    if plot_type == "box_grid":
        make_box_plot_grid(pos_scores, plot_title, threshold, n_splits, wild_index)
    else:
        print("Unsupported plot type in config.")

if __name__ == "__main__":
    main()
