import single_glycan as sg
import sys
import csv
import os
import shutil
import multiprocessing
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
import time

'''
    Helper function to copy target antigen related files from source folder to destination folder   
'''
def copy_antigen_files(source_dir, dest_dir, antigen_pdf, antigen_fasta):
        files = [antigen_fasta, antigen_pdf]
        for item in files:
            source_path = os.path.join(source_dir, item)
            dest_path = os.path.join(dest_dir, item)

            if os.path.isfile(source_path):
                shutil.copy2(source_path, dest_path)



def worker(task_queue):

    while True:
        task = task_queue.get()
        
        if task is None:  # Exit signal
            break
        
        # Unpack the task arguments
        ith_pos, n_struct, antigen, antigen_pdb, antigen_fasta, native_glycan_pos, template_folder, antigen_folder = task

        #for i in range(start_pos, end_pos + 1):
        print(f"thread {ith_pos} starts to work")
                  
        glycan_pos = ith_pos
                        
        folder_for_pos = "pos_" + str(glycan_pos)
        if not os.path.exists(folder_for_pos):
            #shutil.rmtree(folder_for_pos) 
            shutil.copytree(template_folder, folder_for_pos)

            copy_antigen_files(source_dir=antigen_folder, dest_dir= folder_for_pos, antigen_pdf=antigen_pdb, antigen_fasta=antigen_fasta)
                        
        # cd to the folder to a spefic pos
        os.chdir(folder_for_pos)
        if native_glycan_pos == "":
            native_glycan_pos = -1
        else:
            native_glycan_pos = native_glycan_pos
        single_Glycan = sg.Single_Glycan(glycan_pos, n_struct, ag_fasta=antigen_fasta, ag_pdb=antigen_pdb, native_glycan_pos=native_glycan_pos)
        single_Glycan.run_single_glycan()

        # go back to the antigen-level folder
        os.chdir("..")

        #print(f'Process ID: {pid} is now set to run on Core {current_core}')
        return f"---------------Glycan was added to {ith_pos} of {antigen}-----------------------"
        
def convert_string_to_int_list(s):
    return [int(num.strip()) for num in s.split(",") if num.strip().isdigit()]
        

def main(argv):
        start_time = time.time()

        if len(argv) != 4:
            print("please provide number of cores you want to use, starting and ending aa index for glycan(1-based), number of designed structure you want, csv_file")
            print("Example usage: python local_parallel_glycan_run.py 10 12,15 100 antigens.csv")
            sys.exit(2)
        else:
            working_dir = os.getcwd()
            n_cores = int(sys.argv[1])
            positions = convert_string_to_int_list(sys.argv[2])
            n_struct = int(sys.argv[3])
            antigens_csv = sys.argv[4]
            template_folder = os.path.join(os.getcwd(), 'glycan_template_files')
            antigen_folder = os.path.join(os.getcwd(), 'antigens')

            with open(antigens_csv, mode ='r') as file:    
                csvFile = csv.DictReader(file)
                for lines in csvFile:
                    antigen_pdb = lines['pdb']
                    antigen_fasta = lines['fasta']
                    antigen = lines['antigen']
                    native_glycan_pos = lines['native_glycan_pos']

                    folder_for_antigen = antigen
                    path_for_antigen_folder = os.path.join(os.getcwd(), folder_for_antigen)

                    if not os.path.isdir(path_for_antigen_folder):
                         os.makedirs(folder_for_antigen)
                    # cd to the new folder for the antigen
                    os.chdir(folder_for_antigen)
                    print(">>>>>>>Current folder: {folder}".format(folder = folder_for_antigen))

                    num_workers = n_cores
                    print(f"there are {n_cores} processes")
                    task_queue = multiprocessing.Queue()
                    #num_workers = multiprocessing.cpu_count()  # Get the number of available CPU cores

                    # Start multiple worker processes
                    processes = []
                    for _ in range(num_workers):
                        process = multiprocessing.Process(target=worker, args=(task_queue,))
                        processes.append(process)
                        process.start()

                    # Add tasks to the queue (each task is a tuple)
                    for i in positions:  # Change to add more tasks if needed
                        task_queue.put((i, n_struct, antigen, antigen_pdb, antigen_fasta, native_glycan_pos, template_folder, antigen_folder))  # Task ID, Arg1, Arg2

                    # Signal the workers to exit after all tasks are done
                    for _ in range(num_workers):
                        task_queue.put(None)  # Send a sentinel value for each worker

                    # Wait for all processes to finish
                    for process in processes:
                        process.join()
    
                    print(f"All tasks completed for {antigen}")

                    # go back to the working dir
                    os.chdir("..")

                print("ALL JOBS ARE DONE. CURRENT WOKING DIR IS : {folder}".format(folder = os.getcwd()))
                os.chdir(working_dir)

                end_time = time.time()

                execution_time = (end_time - start_time) /60
                print(f"Execution time: {execution_time} minutes.")
                          


if __name__ == "__main__":
        main(sys.argv[1:])
