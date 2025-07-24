import os
import csv
import matplotlib.pyplot as plt
import numpy as np
import sys


ab_names = ["H7-200", "m826", "FluA-20", "H7point5", "H7-167","L3A-44", "L4A-14","HNIgGA6","H7-235", "H7-235_update"]
binding_poses = [4,5,7,8,13,16,22,34,37,40,53,60,62,65,83,85,87,92,94,96,103,115,116,118,120,122,125,127,133,146,152,153,157,159,161,165,167,168,170,174,176,178,179,180,181,184,185,187,196,200,202,221,222,223,228,230,231,232,233,234,235,236,240,245,249,251,253,254,256,260,262,267]
wild_index = 200

def safe_open_file(filename, mode='r'):
    try:
        return open(filename, mode)
    except FileNotFoundError:
        print(f"File not found: {filename}")
        return None
    except IOError as e:
        print(f"I/O error({e.errno}): {e.strerror}")
        return None

def count_ab_files():
    working_dir = os.getcwd()
    for pos in binding_poses:
        pos_dir = os.path.join(working_dir, f"pos_{pos}")
        pos_count = {}
        ab_count = {}
        for ab_name in ab_names:
            prefix = f"Add_Rlx_{ab_name}"
            suffix = ".pdb"
            count = sum(1 for fname in os.listdir(pos_dir) if fname.startswith(prefix) and fname.endswith(suffix))
            ab_count[ab_name] = count
        pos_count[pos] = ab_count

    return pos_count

def print_count(pos_count):
    results = []
    total = 0
    poses = []
    for pos in pos_count.keys():
        ab_count = pos_count[pos]
        for ab_name in ab_count.keys():
            cnt = ab_count[ab_name]
            if cnt < 75:
                result = f"pos_{pos}: {ab_count}"
                print(result)
                results.append(result)
                total = total + 1
                poses.append(pos)
                break

    print(f"Thre are total of {total} that have less han 75")
    print(poses)
    return results

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
                        values.append(parts[column_index])
                    except IndexError:
                        continue  # Skip malformed lines

    return values

def count_column_values(sc_file):
    column_counts = {}

    f = safe_open_file(sc_file)

    if f is not None:
        with f:  
            for line in f:
                if line.startswith("SCORE:") and 'total_score' in line:
                    header = line.strip().split()
                    for column in header:
                        column_counts[column] = 0
                    break

            for line in f:
                if line.startswith("SCORE:"):
                    parts = line.strip().split()
                    for i, value in enumerate(parts):
                        if i < len(header):
                            column = header[i]
                            column_counts[column] += 1

    return column_counts

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
      
def count_in_scFiles(sc_type):
    working_dir = os.getcwd()
    pos_count = {}
    for pos in binding_poses:
        pos_dir = os.path.join(working_dir, f"pos_{pos}")

        ab_count = {}
        for ab_name in ab_names:
            sc_file = get_sc_file_path(sc_type, ab_name, pos_dir)
            count = len(get_column_values(sc_file, "total_score"))
            ab_count[ab_name] = count
        pos_count[pos] = ab_count

    return pos_count
    

def main(argv):

    cur_dir = os.getcwd()
    sc_type = sys.argv[1]
    count = count_in_scFiles(sc_type)
    print_count(count)
    
    #binding_file_count = count_ab_files()
    #print_count(binding_file_count)
    
    

if __name__ == "__main__":
        main(sys.argv[1:])
