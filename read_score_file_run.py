import matplotlib.pyplot as plt
import os
import sys
import statistics
import numpy as np

def get_top_n(top_n, sc_file):
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

def main(argv):

    cur_dir = os.getcwd()

    pos = sys.argv[1]
    top_n = int(sys.argv[2])
    target_dir = f = os.path.join(cur_dir, "pos_"+pos)
    score_file_path = os.path.join(target_dir, "Glyc_score.sc")
    
    _,_,scores = get_top_n(top_n, score_file_path)
    
    l=len(scores)
    print(f"the length of the score file is {l}")
    for s in scores:
        print(f"{s},")
    print(f"the average score of pos {pos} is {statistics.mean(scores)} ")




if __name__ == "__main__":
    main(sys.argv[1:])
