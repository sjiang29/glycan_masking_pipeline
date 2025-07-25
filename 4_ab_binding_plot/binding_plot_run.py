import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import sys
                   
import statistics                

binding_poses = [4,5,7,8,13,16,22,34,37,40,53,60,62,65,83,85,87,92,94,96,103,115,116,118,120,122,125,127,133,146,152,153,157,159,161,165,167,168,170,174,176,178,179,180,181,184,185,187,196,200,202,221,222,223,228,230,231,232,233,234,235,236,240,245,249,251,253,254,256,260,262,267]
wild_index = 200
ab_names = ["H7-200", "m826", "FluA-20", "H7point5", "H7-167","L3A-44", "L4A-14","HNIgGA6","H7-235_update"]



def get_score(sc_file):
    #print(sc_file)
    
    file = open(sc_file) 
    line = file.readline()
    
    cnt = 1
    scores = []
    
    while line:
        if cnt >= 3 and cnt<=102:
            parts = line.split(":")[1].split()
            score = parts[0].strip()
            
            try:
                score = float(score)
                
            except ValueError:
                #print(parts[4])
                pass
                       
            scores.append(score)
            
        line = file.readline()
        cnt = cnt +1
    file.close()

    return scores

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
                    except (IndexError, ValueError) as e:
                        continue  # Skip malformed linespython parallel_pos_binding_run.py 5,7 2 updated_9bt5.csv 2

    return values

def box_plot_dict(scores_for_abs, pos):
     
    vals = list(scores_for_abs.values())

    # Create the box plot
    labels = scores_for_abs.keys()
    plt.figure(figsize = (15, 6))
    plt.boxplot(vals, labels=labels)
    plt.title(f"ab binding for pos_{pos}", fontsize = 16)
    #plt.figure(figsize = (15, 6))
    plt.xlabel('Ab', fontsize = 16)
    #plt.xticks(np.arange(101, 150, 10))
    plt.xticks(rotation=90)
    plt.tick_params(labelsize = 12)
    plt.ylim(np.min(vals), np.max(vals))
    plt.ylabel('Rosetta_score', fontsize = 16)
    plt.show()

def box_pot_with_wild(scores_for_pos, pos, scores_for_wild):
     
    labels = sorted(scores_for_pos.keys())
    
    data1 = [scores_for_wild[label] for label in labels]
    data2 = [scores_for_pos[label] for label in labels]
    
    
    x = np.arange(len(labels))
    width = 0.3

    fig, ax = plt.subplots(figsize=(10, 5))

    # Plot boxplots side by side by shifting positions
    box1 = ax.boxplot(data1, positions=x - width/2, widths=0.25, patch_artist=True)
    box2 = ax.boxplot(data2, positions=x + width/2, widths=0.25, patch_artist=True)

    # Color boxes for clarity
    for patch in box1['boxes']:
        patch.set(facecolor='lightblue')
    for patch in box2['boxes']:
        patch.set(facecolor='lightgreen')

    # Labels and formatting
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_xlabel('Antibody')
    ax.set_ylabel('REU')
    ax.set_title(f'Box Plot Comparison: wildtype vs pos_{pos} ')
    ax.legend([box1["boxes"][0], box2["boxes"][0]], ['wildtype',f'pos_{pos}'], loc='upper left')

    plt.tight_layout()
    plt.show()

def get_sc_file_path(sc_type, ab_name, pos_dir):
    if sc_type == "binding":
        prefix = f"Add_Rlx_{ab_name}"

        binding_sc_file = f"{prefix}score.fasc"
        binding_sc_file_path = os.path.join(pos_dir, binding_sc_file)
        return binding_sc_file_path
    if sc_type == "docking":
        docking_sc_file = f"InAn_cleaned_{ab_name}.sc"
        docking_sc_file_path = os.path.join(pos_dir, docking_sc_file)
        return docking_sc_file_path
      

def get_single_pos_score(pos, sc_type, column_name):
    working_dir = os.getcwd()
    pos_dir = os.path.join(working_dir, f"pos_{pos}")
    #abs_file_path = os.path.join(pos_dir, abs_csv)

    scores = {}

    for ab_name in ab_names:
        
        sc_file_path = get_sc_file_path(sc_type, ab_name, pos_dir)
        try:
            f = open(sc_file_path)
            score = get_column_values(sc_file_path, column_name)
            #print(score)
            if len(score) == 75:
                scores[ab_name] = score
            else:
                return -1

        except Exception:
            return -1
        
    return scores

def make_box_plot_grid(pos_and_scores, ab_name,  n_splits, sc_type, column_name):
    
    keys = list(pos_and_scores.keys())  # Extract amino acid positions
    values = list(pos_and_scores.values())

    # Compute the mean of values when key == 0(wild type)
    #wild_type = 200
    avg_value_of_wild = np.mean(pos_and_scores[wild_index]) if wild_index in pos_and_scores else None
    
    # Split data into n_splits parts
    split_size = len(keys) // n_splits  # Determine approximate split size
    remainder = len(keys) % n_splits  # Handle uneven splits
    
    data_splits = []
    start = 0
    for i in range(n_splits):
        end = start + split_size + (1 if i < remainder else 0)  # Distribute remainder across first few splits
        data_splits.append((keys[start:end], values[start:end], f"Segment {i+1}"))
        start = end
    
    fig, axes = plt.subplots(n_splits, 1, figsize=(16, 3.5 * n_splits))  # Adjust figure size based on n_splits
    
    if n_splits == 1:
        axes = [axes]  # Ensure axes is iterable for a single subplot
    
    for i, (keys_split, values_split, title) in enumerate(data_splits):
        box = axes[i].boxplot(values_split, tick_labels=keys_split, patch_artist=True)
        
        # Set grey color for each box
        for patch in box['boxes']:
            patch.set(facecolor='grey')
        
        # Formatting
        axes[i].set_xticklabels(keys_split, rotation=0, fontsize=8)
        axes[i].tick_params(labelsize=8)
        axes[i].set_ylim(np.min(values), 0)
        axes[i].set_ylabel(f'{column_name}', fontsize=14)
        

        # Add horizontal line at avg_value (if available)
        if avg_value_of_wild is not None:
            axes[i].axhline(avg_value_of_wild, color='r', linestyle='dotted', linewidth=1, label=f"Mean (Pos 0) = {avg_value_of_wild:.2f}")
            axes[i].legend(fontsize=8)

    plt.xlabel('Amino Acid Position', fontsize=14)
    plt.suptitle(f'{sc_type} for {ab_name} across all pos', fontsize=10)
    plt.subplots_adjust(bottom=0.3)  # Ensure x-axis labels are fully visible
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])  # Adjust layout to prevent overlap
    plt.show()

def make_violin_plot_grid(pos_and_scores, ab_name, n_splits, sc_type, column_name):
    keys = list(pos_and_scores.keys())  # Extract amino acid positions
    values = list(pos_and_scores.values())

    # Separate wild type (index 200) and move it to the beginning
    if wild_index in pos_and_scores:
        wild_idx = keys.index(wild_index)
        wild_value = values[wild_idx]

    # Remove wild type entry from original lists
        del keys[wild_idx]
        del values[wild_idx]

    # Prepend wild type to both
        keys = [wild_index] + keys
        values = [wild_value] + values
        avg_value_of_wild = np.mean(wild_value)
    else:
        avg_value_of_wild = None

    # Split data into n_splits parts
    split_size = len(keys) // n_splits
    remainder = len(keys) % n_splits

    data_splits = []
    start = 0
    for i in range(n_splits):
        end = start + split_size + (1 if i < remainder else 0)
        data_splits.append((keys[start:end], values[start:end], f"Segment {i+1}"))
        start = end

    fig, axes = plt.subplots(n_splits, 1, figsize=(16, 3.5 * n_splits))

    if n_splits == 1:
        axes = [axes]

    for i, (keys_split, values_split, title) in enumerate(data_splits):
        pos = np.arange(1, len(keys_split) + 1)

        # Create violin plot
        axes[i].violinplot(values_split, positions=pos, showmeans=False, showmedians=True, widths=0.7)

        axes[i].set_xticks(pos)
        axes[i].set_xticklabels(keys_split, rotation=0, fontsize=6)
        axes[i].tick_params(labelsize=8)
        axes[i].set_ylim(np.min(values), 0)
        axes[i].set_ylabel(f'{column_name}', fontsize=10)

        # Add horizontal line for wild-type mean
        if avg_value_of_wild is not None:
            axes[i].axhline(avg_value_of_wild, color='r', linestyle='dotted', linewidth=1, label=f"Mean (WT) = {avg_value_of_wild:.2f}")
            axes[i].legend(fontsize=6)

    plt.xlabel('Amino Acid Position', fontsize=14)
    plt.suptitle(f'{sc_type} for {ab_name} across all pos', fontsize=8)
    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.show()

def get_score_for_ab(ab_name,sc_type, column_name):
     binding_for_ab = {}
     for pos in binding_poses:
          score_for_pos = get_single_pos_score(pos, sc_type, column_name)
          if score_for_pos != -1:
               binding_for_ab[pos] = score_for_pos[ab_name]

     return binding_for_ab
               
def print_dict(ab_name, pos_and_score, column):
    print(f"----- There are {len(pos_and_score)} amino acid in the dict for antibody:{ab_name}")
    for k,v in pos_and_score.items():
        print(f"pos {k}: avg(mean) {column} = {statistics.mean(v):.2f}, median {column} = {statistics.median(v):.2f}")

def make_violin_plot_for_all(all_ab_scores_dict, sc_type, column_name, wild_index=200):
    n_abs = len(all_ab_scores_dict)
    fig, axes = plt.subplots(n_abs, 1, figsize=(16, 3.5 * n_abs), squeeze=False)
    
    for i, (ab_name, pos_and_scores) in enumerate(all_ab_scores_dict.items()):
        keys = list(pos_and_scores.keys())
        values = list(pos_and_scores.values())

        # Move wild type to the beginning
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
        #ax.set_ylabel(f'{column_name}', fontsize=8)
        #ax.set_title(f'{ab_name}', fontsize=8)

        if avg_value_of_wild is not None:
            ax.axhline(avg_value_of_wild, color='r', linestyle='dotted', linewidth=1, label=f"Mean (WT): {avg_value_of_wild:.2f}, Ab:{ab_name}")
            ax.legend(fontsize=8)

        # Hide x-axis labels for all but last subplot
        if i < n_abs - 1:
            ax.set_xticklabels([])
            ax.set_xlabel('')

    # One shared y-axis label
    fig.text(0.04, 0.5, column_name, va='center', rotation='vertical', fontsize=12)

    # Shared x-axis label
    fig.text(0.5, 0.01, 'Amino Acid Position', ha='center', fontsize=12)

    plt.tight_layout(rect=[0.05, 0.05, 1, 0.95])  # Leave space for shared labels
    plt.show()
    
                   
def main(argv):

    cur_dir = os.getcwd()
    sc_type = sys.argv[1]
    pos = int(sys.argv[2])
    column_name = sys.argv[3]
    #abs_csv = sys.argv[2]
    #ab_name = sys.argv[3]
    
    scores_for_pos = get_single_pos_score(pos, sc_type, column_name)
    scores_for_wild = get_single_pos_score(wild_index, sc_type, column_name)
    #print_dict(scores_for_pos)
    #box_pot_with_wild(scores_for_pos, pos, scores_for_wild)
    all_scores = {}
    for ab_name in ab_names:
        binding_for_ab = get_score_for_ab(ab_name, sc_type, column_name)
        print_dict(ab_name, binding_for_ab, column_name)
        all_scores[ab_name] = binding_for_ab
        #
    make_violin_plot_for_all(all_scores, sc_type, column_name,wild_index=200)

    
    

if __name__ == "__main__":
        main(sys.argv[1:])

    