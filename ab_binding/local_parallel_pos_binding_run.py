import single_ab_binding as sab
import sys
import os
import csv
import shutil
from multiprocessing import Pool
import multiprocessing

class SinglePosBinding:
    sc_file = "Glyc_score.sc"
    template_folder = "binding_template_files"
    #Abs_folder = "Abs_pdb"

    def __init__(self, pos: int, top_n: int, abs_csv: str):
        self.pos = pos
        self.abs_csv = abs_csv
        self.top_n = top_n

    def run_binding_for_all_abs(self):
        """Runs binding for all antibodies in the CSV file."""
        
        with open(self.abs_csv, mode='r') as file:
            csv_reader = csv.DictReader(file)
            for row in csv_reader:
                ab = row['antibody']
                ab_heavy_pdb = row['antibody_heavy_chain_pdb']
                ab_light_pdb = row['antibody_light_chain_pdb']
                antigen_list = self.get_list_file()
                
                single_ab_binding = sab.Single_ab_binding(
                    antigen_list=antigen_list,
                    ab=ab,
                    ab_heavy_pdb=ab_heavy_pdb,
                    ab_light_pdb=ab_light_pdb
                )
                single_ab_binding.run_addChain_cmd()
                

                #create_dir_for_each_abBinding(ab)

    def sort_files_on_score(self):
            
            file = open(SinglePosBinding.sc_file)

            line = file.readline()

            cnt = 1
            file_to_score = {}

            while line:

                if cnt >= 3 and cnt :
                    parts = line.split(" ")
                    parts_s = line.split(":")[1].split()
                    score = parts_s[0].strip()

                    try:
                        score = float(score)
                
                    except ValueError:
                        #print(parts[4])
                        pass
                       
                    file_name = parts[-1][0:-1]
                    file_to_score[file_name] = score
                
                line = file.readline()
                cnt = cnt + 1
            file.close()
            
            # sorted_files is tuple
            sorted_files = sorted(file_to_score.items(), key=lambda item: item[1])
            
            return sorted_files

    def get_top_sc_file(self):
            
        sorted_files = self.sort_files_on_score()

        top_sc_files = []

        if self.top_n > len(sorted_files):
             for i in range(len(top_sc_files)):
                 top_sc_files.append(sorted_files[i][0] + ".pdb")

        for i in range(self.top_n):
            top_sc_files.append(sorted_files[i][0] + ".pdb")

        return top_sc_files
    
    def get_list_file(self):
        """Creates a list file with the top scoring structures."""
        file_name = "antigen_list_file"

        top_sc_files = self.get_top_sc_file()
        with open(file_name, "w") as f:
            for file in top_sc_files:
                f.write(file + "\n")

        return file_name


def parse_positions(pos_str):
    """Parses position ranges like '1-10,15,20-30' into a list of integers."""
    positions = set()
    for part in pos_str.split(","):
        if "-" in part:
            start, end = map(int, part.split("-"))
            positions.update(range(start, end + 1))
        else:
            positions.add(int(part))
    return sorted(positions)


def worker(task_queue):
    """Runs the workflow for a single position."""
    while True:

        task = task_queue.get()
        
        if task is None:  # Exit signal
            break
        
    
        pos, top_n, abs_csv, working_dir = task
        print(f"---------------------Thread for pos_{pos} starts working")

        pos_dir = os.path.join(working_dir, f"pos_{pos}")
        
        if not os.path.exists(pos_dir):
            raise FileNotFoundError(f"Error: The directory '{pos_dir}' does not exist!")

        shutil.copy(abs_csv, pos_dir)
        os.chdir(pos_dir)

        spb = SinglePosBinding(pos, top_n, abs_csv)

        # Copy template folder
        template_dir = os.path.join(working_dir, SinglePosBinding.template_folder)
        if os.path.exists(template_dir):
            shutil.copytree(template_dir, pos_dir, dirs_exist_ok=True)  # Overwrite if exists

        # Run antibody binding workflow
        spb.run_binding_for_all_abs()
        print(f"---------------------Thread for pos_{pos} ends working-----------------")
        
def create_dir_for_each_abBinding(ab):
    prefix = f"Add_Rlx_{ab}"
    destination = f"{ab}"

    # Create the destination folder if it doesn't exist
    os.makedirs(destination, exist_ok=True)

    # Loop through files in the current directory
    for filename in os.listdir():
        if filename.startswith(prefix) and os.path.isfile(filename):
            #shutil.move() will overwrite the destination file if a file with the same name already exists
            shutil.move(filename, os.path.join(destination, filename))
    


def main(argv):
    """Main function to process input arguments and run jobs in parallel."""
    if len(argv) != 4:
        print("Usage: python local_parallel_pos_binding.py <num_cores> <positions> <top_n> <abs_csv>")
        print("Example: python local_parallel_pos_binding.py 4 1-10,15,20-30 75 Abs_all.csv 4")
        sys.exit(2)

    num_cores, positions_str, top_n, abs_csv = argv
    top_n = int(top_n)
    num_workers = int(num_cores)

    positions = parse_positions(positions_str)

    working_dir = os.getcwd()

    task_queue = multiprocessing.Queue()
                    
    # Start multiple worker processes
    processes = []
    for _ in range(num_workers):
        process = multiprocessing.Process(target=worker, args=(task_queue,))
        processes.append(process)
        process.start()

                    # Add tasks to the queue (each task is a tuple)
    for i in positions:  # Change to add more tasks if needed
        task_queue.put((i, top_n, abs_csv, working_dir))  # Task ID, Arg1, Arg2

    # Signal the workers to exit after all tasks are done
    for _ in range(num_workers):
        task_queue.put(None)  # Send a sentinel value for each worker

    # Wait for all processes to finish
    for process in processes:
        process.join()

    return


if __name__ == "__main__":
    main(sys.argv[1:])
