import matplotlib.pyplot as plt
import os
import sys
import statistics
import numpy as np

wild_index = 200
def get_top_n(top_n, sc_file):
    '''
    Extracts the top n scores from a given Rosetta score file.

    Parameters:

        top_n (int): The number of top scores to extract.

        sc_file (str): The path to the score file.

    Returns:

        list: The top n scores from the file.

        list: The top 50 scores from the first 100 entries.

        list: All extracted scores.
    '''
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
    
    scores_100 = scores[0:100]
    scores_100.sort()
    scores.sort()
    #print(scores)   
    return scores[0:top_n], scores_100[0:50],scores

def get_top_n_for_all(working_dir, sc_file, top_n, start, end):
    '''
    Extracts and organizes top scores from multiple subdirectories in a given working directory.

    Parameters:

        working_dir (str): The root directory containing position-based subdirectories.

        sc_file (str): The filename of the score file.

        top_n (int): The number of top scores to extract.

        start (int): The starting amino acid position.

        end (int): The ending amino acid position.

    Returns:

        dict: Ordered dictionary of top n scores for each position.

        dict: Ordered dictionary of top 50 scores from 100 selections.

        dict: Dictionary mapping positions to average scores.
    '''
    pos_and_scores = {}
    pos_and_scores_50 = {}
    pos_and_avg_score = {}
    for filename in os.listdir(working_dir):
        f = os.path.join(working_dir, filename)
        # checking if it is a file
        print(f)
        if os.path.isdir(f) and f.find("pos_")!= -1:
            score_file_path = os.path.join(f, sc_file)
            if os.path.exists(score_file_path):
                
                #top_n_scores = get_top_n(top_n, score_file_path)
                
                pos = f.split("/")[-1]
                pos = pos.split("_")[-1]

                top_n_scores, scores_50, scores = get_top_n(top_n, score_file_path)
                avg_score = statistics.mean(scores)
                pos_and_avg_score[int(pos)] = avg_score
               
                if (int(pos) >= start and int(pos) <= end):
                    #top_n_scores, scores_50 = get_top_n(top_n, score_file_path)
                    
                    pos_and_scores[int(pos)] = top_n_scores
                    pos_and_scores_50[int(pos)] = scores_50

                if int(pos) == wild_index:
                    pos_and_scores[int(pos)] = top_n_scores
                    pos_and_scores_50[int(pos)] = scores_50

    # sort the dict by the keys in ascending order               
    ordered_all = dict(sorted(pos_and_scores.items()))
    ordered_50 = dict(sorted(pos_and_scores_50.items()))
                    
    return ordered_all, ordered_50, pos_and_avg_score

def compare_top75_and_top50from100(working_dir, sc_file):

    ordered_all, ordered_50, pos_and_avg_score = get_top_n_for_all(working_dir, sc_file, 75, 1, 268)

    for k, v in ordered_all.items():
        v = statistics.mean(v)
        ordered_all[k] = v

    for k,v in ordered_50.items():
        v = statistics.mean(v)
        ordered_50[k] = v
    x1,y1 = ordered_all.keys(), ordered_all.values()
    x2,y2 = ordered_50.keys(), ordered_50.values()

    fig = plt.figure()

    ax = fig.add_subplot(111)
    #ax.scatter(x1, y1, c = 'b', marker = 's', label = "top 75 from 150")
    #ax.scatter(x2, y2, c = 'r', marker = 'o', label = 'top 50 from 100')
    plt.scatter(x=y1, y=y2)
    plt.xlabel("Average_Rosetta_score \n of top 75 from 150", fontsize = 16)
    plt.ylabel("Average_Rosetta_score \n of top 50 from 100", fontsize = 16)
    
    #plt.legend(loc='upper left')
    plt.figure(figsize=(10,10))
    plt.show()

'''
grid of box plot, can be adjusted accordingly
'''
def make_box_plot_grid(pos_and_scores, dict_name, threshold, n_splits):
    pos_and_scores = pick_pos_using_threshold(pos_and_scores, threshold)
    print_dict(pos_and_scores)
    keys = list(pos_and_scores.keys())  # Extract amino acid positions
    values = list(pos_and_scores.values())

    # Compute the mean of values when key == 0(wild type)
    wild_type = 200
    avg_value_of_wild = np.mean(pos_and_scores[wild_type]) if wild_type in pos_and_scores else None
    
    # Split data into n_splits parts
    split_size = len(keys) // n_splits  # Determine approximate split size
    remainder = len(keys) % n_splits  # Handle uneven splits
    
    data_splits = []
    start = 0
    for i in range(n_splits):
        end = start + split_size + (1 if i < remainder else 0)  # Distribute remainder across first few splits
        data_splits.append((keys[start:end], values[start:end], f"Segment {i+1}"))
        start = end
    
    fig, axes = plt.subplots(n_splits, 1, figsize=(10, 3.5 * n_splits))  # Adjust figure size based on n_splits
    
    if n_splits == 1:
        axes = [axes]  # Ensure axes is iterable for a single subplot
    
    for i, (keys_split, values_split, title) in enumerate(data_splits):
        box = axes[i].boxplot(values_split, labels=keys_split, patch_artist=True)
        
        # Set grey color for each box
        for patch in box['boxes']:
            patch.set(facecolor='grey')
        
        # Formatting
        axes[i].set_xticklabels(keys_split, rotation=0, fontsize=12)
        axes[i].tick_params(labelsize=12)
        axes[i].set_ylim(np.min(values) - 1000, threshold)
        axes[i].set_ylabel('Rosetta Score', fontsize=14)

        # Add horizontal line at avg_value (if available)
        if avg_value_of_wild is not None:
            axes[i].axhline(avg_value_of_wild, color='r', linestyle='dotted', linewidth=1, label=f"Mean(WT:pos200) = {avg_value_of_wild:.2f}")
            axes[i].legend(fontsize=8)

        
    # Global figure formatting
    plt.xlabel('Amino Acid Position', fontsize=14)
    plt.suptitle(f"Influenza_{dict_name}", fontsize=10)
    plt.subplots_adjust(bottom=0.3)  # Ensure x-axis labels are fully visible
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])  # Adjust layout to prevent overlap
    plt.show()


def make_box_plot(pos_and_scores, dict_name, threshold):

    pos_and_scores = pick_pos_using_threshold(pos_and_scores, threshold)
    vals = list(pos_and_scores.values())

    # Create the box plot
    labels = pos_and_scores.keys()
    plt.figure(figsize = (15, 6))
    plt.boxplot(vals, labels=labels)
    plt.title(dict_name, fontsize = 16)
    #plt.figure(figsize = (15, 6))
    plt.xlabel('Amino_acid_position', fontsize = 16)
    #plt.xticks(np.arange(101, 150, 10))
    plt.xticks(rotation=90)
    plt.tick_params(labelsize = 12)
    plt.ylim(np.min(vals), np.max(vals))
    plt.ylabel('Rosetta_score', fontsize = 16)
    plt.show()

def make_violin_plot_grid(pos_and_scores, dict_name, threshold, n_splits):
    pos_and_scores = pick_pos_using_threshold(pos_and_scores, threshold)
    print_dict(pos_and_scores)
    keys = list(pos_and_scores.keys())  # Extract amino acid positions
    values = list(pos_and_scores.values())

    wild_type = 200
    avg_value_of_wild = np.mean(pos_and_scores[wild_type]) if wild_type in pos_and_scores else None

    # Split data into n_splits
    split_size = len(keys) // n_splits
    remainder = len(keys) % n_splits

    data_splits = []
    start = 0
    for i in range(n_splits):
        end = start + split_size + (1 if i < remainder else 0)
        data_splits.append((keys[start:end], values[start:end]))
        start = end

    fig, axes = plt.subplots(n_splits, 1, figsize=(16, 3.5 * n_splits))

    if n_splits == 1:
        axes = [axes]

    for i, (keys_split, values_split) in enumerate(data_splits):
        parts = axes[i].violinplot(values_split, showmeans=False, showextrema=True, showmedians=True)
        
        # Format violin parts (e.g., change color)
        for pc in parts['bodies']:
            pc.set_facecolor('lightgrey')
            pc.set_edgecolor('black')
            pc.set_alpha(0.8)

        if 'cmedians' in parts:
            parts['cmedians'].set_color('black')

        # Formatting
        axes[i].set_xticks(np.arange(1, len(keys_split) + 1))
        axes[i].set_xticklabels(keys_split, rotation=0, fontsize=8)
        axes[i].tick_params(labelsize=8)
        axes[i].set_ylim(np.min(values) -200, threshold)
        axes[i].set_ylabel('Rosetta Score', fontsize=14)
        

        # Add horizontal line for wild type average
        if avg_value_of_wild is not None:
            axes[i].axhline(avg_value_of_wild, color='red', linestyle='dotted', linewidth=1,
                            label=f"Mean (WT) = {avg_value_of_wild:.2f}")
            axes[i].legend(fontsize=8)

    plt.xlabel('Amino Acid Position', fontsize=14)
    plt.suptitle(f"Influenza_{dict_name}", fontsize=12)
    plt.subplots_adjust(bottom=0.3)
    plt.tight_layout(rect=[0, 0.03, 1, 0.97])
    plt.show()

def pick_pos_using_threshold(pos_and_scores, threshold):

    updated_pos_and_scores = {}
    for k, v in pos_and_scores.items():
        if statistics.mean(v) < threshold:
            updated_pos_and_scores[k] = v

    return updated_pos_and_scores

def print_dict(pos_and_score):
    print(f"----- There are {len(pos_and_score)} amino acid in the dict")
    for k,v in pos_and_score.items():
        print(f"pos {k}: avg(mean) score = {statistics.mean(v):.2f}, median score = {statistics.median(v):.2f}")
    
def plot_all_for_current_dir(cur_dir, top_n, start, end, threshold):
    pos_and_scores, pos_and_scores_50, pos_and_avg_score = get_top_n_for_all(cur_dir, "Glyc_score.sc", top_n, start, end)
   
    #pos_and_sasa = read_sasa(sasa_file)
    #plot_sasa_vs_score(pos_and_sasa, pos_and_avg_score)
    #make_box_plot(pos_and_scores, f"Top {top_n} from 100")
    num_of_splits = 1
    make_box_plot_grid(pos_and_scores, f"Top {top_n} from 100 with RSU threshold of {threshold}", threshold, num_of_splits)
    #make_violin_plot_grid(pos_and_scores, f"Top {top_n} from 100 with RSU threshold of {threshold}", threshold, num_of_splits)
    # make_box_plot(pos_and_scores_50, "Top 50 from 100 out of 150")

def read_sasa(sasa_file):
    file = open(sasa_file)
    line = file.readline()
    pos_and_sasa = {}
    while line:
        parts = line.split(":")

        pos = parts[0]
        pos = int(pos.split(" ")[0][3:])
        score = float(parts[1])
        pos_and_sasa[pos] = score

        line = file.readline()
    file.close()

    return pos_and_sasa

def plot_sasa_vs_score(pos_and_sasa, pos_and_avg_score):

    sasa_vs_score = {}

    for key,value in pos_and_avg_score.items():
        sasa_vs_score[pos_and_sasa[key]] = value

    x,y = sasa_vs_score.keys(), sasa_vs_score.values()

    #plt.figure(figsize=())
    plt.scatter(x, y)
    plt.xlabel("Solvent_Accessible_Surface_Area", fontsize = 16 )
    plt.ylabel("Average_Rosetta_score", fontsize = 16)
    plt.tick_params(labelsize = 14)
    plt.ylim(min(y), max(y))
    plt.show()

   

def main(argv):

    cur_dir = os.getcwd()
    
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    top_n = int(sys.argv[3])
    threshold = int(sys.argv[4])
    #sasa_file = sys.argv[5]
    plot_all_for_current_dir(cur_dir, top_n, start, end, threshold)

    #compare_top75_and_top50from100(cur_dir, "Glyc_score.sc")
    
    

if __name__ == "__main__":
        main(sys.argv[1:])