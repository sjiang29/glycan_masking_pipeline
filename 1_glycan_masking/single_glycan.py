import xml.etree.ElementTree as ET
import shutil
import sys
import subprocess
from subprocess import run
from shlex import split
import os
import stat



class Single_Glycan:

    #origin_ag = "TVERTNIPRICSKGKRTVDLGQCGLLGTITGPPQCDQFLEFSADLIIERREGSDVCYPGKFVNEEALRQILRESGGIDKEAMGFTYSGIRTNGATSACRRSGSSFYAEMKWLLSNTDNAAFPQMTKSYKNTRKSPALIVWGIHHSVSTAEQTKLYGSGNKLVTVGSSNYQQSFVPSPGARPQVNGLSGRIDFHWLMLNPNDTVTFSFNGAFIAPDRASFLRGKSMGIQSGVQVDANCEGDCYHSGGTIISNLPFQNIDSRAVGKCPRYV"
    #origin_ag_pdb = "head_6uig_cut_ABC_A.pdb"

    fast_design_options = "FastDesign.options"
    #example_resfile = "example.resfile"
    example_FastDesign_xml = "template_FastDesign.xml"
    fast_design_sh = "fast_design.sh"

    GlycanTreeModeler_options = "GlycanTreeModeler.options"
    example_GlycanTreeModeler_xml = "template_GlycanTreeModeler.xml"
    GlycanTreeModeler_sh = "GlycanTreeModeler.sh"

    '''
    Constructor for the Sinle_Glycan class

    Args:

        glycan_pos : an int to define the index of aa(N) for glycan adding
        n_struct : # of desired output designed structures
        parent_dir: path of parent dir which contains this python file
        all indexed are 1-based
    '''
    def __init__(self, glycan_pos: int, n_struct: int, ag_pdb:str, ag_fasta:str, native_glycan_pos:str):
        self.glycan_pos = glycan_pos
        self.n_struct = n_struct
        self.ag_pdb = ag_pdb
        self.ag_fasta = ag_fasta
    
        self.origin_ag_seq = self.read_ag_fasta()
        self.native_glycan_pos = native_glycan_pos.replace("+",",")
        
    
    
    def read_ag_fasta(self):
        
        sequences = {}
        with open(self.ag_fasta, 'r') as file:
            header = ""
            sequence = ""
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if header:
                        sequences[header] = sequence
                    header = line[1:]
                    sequence = ""
                else:
                    sequence += line
            if header:
                sequences[header] = sequence
        # partation the fasta file name so to get the header name
        parts = self.ag_fasta.split('.fasta')
        ag_header = parts[0]
        
        return sequences[ag_header]

    '''
    Helper function to confirm how to modify corresponding amino acids so to make plycan_pos desirable for glycan-attachment

    Returns:

        A dictionary: key of which is the amino acid index(1-based), value of which is the wanted amino acid in such position
    '''

    def confirm_mod_pos(self):
        ag_len = len(self.origin_ag_seq)
        mod_pos = {}
        s_pos = self.glycan_pos
        
        if s_pos >= 1 and s_pos <= (ag_len - 2):
            if self.origin_ag_seq[s_pos - 1] != "N":
                mod_pos[s_pos] = "N"
            if self.origin_ag_seq[s_pos] == "P":
                mod_pos[s_pos + 1] = "A"
            if self.origin_ag_seq[s_pos + 1] != "T":
                mod_pos[s_pos + 2] = "T"
        else:
            assert("starting position is not correct")

        return mod_pos


    '''
    Helper function to create a resfile needed for clarifying which amino acid is needed to be modified to what, such resfile is needed for running fast design protocol

    Returns:

        Name of the resfile: designed to be s{glycan_pos}.resfile, e.q. s1.resfile
    '''
    def create_resfile_mod_seq(self):

        mod_pos = self.confirm_mod_pos()
        file_content_head = "NATRO\nEX 1 EX 2\nstart\n"
        file_content_body = ""
        for key in mod_pos:
            file_content_body = file_content_body + str(key) + " " + "A" + " " + "PIKAA" + " "+str(mod_pos[key]).upper() + '\n'

        file_content = file_content_head + file_content_body

        resfile_name = "s"+ str(self.glycan_pos)+ ".resfile"

        f = open(resfile_name, "w")
        f.write(file_content)
        f.close()

        return resfile_name

    '''
    Helper function to create FastDesign.xml needed for running fast design protocol

    Returns:

        Name of the fast design xml: FastDesign.xml
    '''
    def create_FastDesign_xml(self):
        # Parse the XML file
        
        tree = ET.parse(Single_Glycan.example_FastDesign_xml)
        print(tree)
        root = tree.getroot()
        print(root)

        # Modify an element
        for resfile in root.iter('ReadResfile'):
            resfile.set("filename", "s"+ str(self.glycan_pos)+ ".resfile")
                
        
        # Write the changes back to the file
               
        fastdesign_file_name = "FastDesign.xml"

        dst = fastdesign_file_name
        
        tree.write(dst)

        return fastdesign_file_name
    
    
    '''
    Helper function to run fast design protocol

    '''

    def run_fast_design_cmd(self):
        print("****start of fast design")
        fast_design_xml = self.create_FastDesign_xml()
        
        st = os.stat(Single_Glycan.fast_design_sh)
        os.chmod(Single_Glycan.fast_design_sh, st.st_mode | stat.S_IEXEC)

        #fast_design_cmd = "./" + Single_Glycan.fast_design_sh

        subprocess.run( f"bash fast_design.sh '{self.ag_pdb}'", shell = True)

        #subprocess.call(fast_design_cmd, shell = True)

        parts = self.ag_pdb.split(".pdb")

        designed_pdb = "Des_" + parts[0] + "_0001.pdb"

        print("*****end of fast design")

        return designed_pdb
    

    '''
    Helper function to create GlycanTreeModeler.xml needed for running  glycan tree modeler protocol

    Returns:

        Name of the fast design xml: FastDesign.xml
    '''

    def create_GlycanTreeModeler_xml(self):
        # Parse the XML file
        
        tree = ET.parse(Single_Glycan.example_GlycanTreeModeler_xml)
        root = tree.getroot()
        postions_for_glycan = ""
        # Modify an element
        if self.native_glycan_pos == "-1":
            postions_for_glycan = str(self.glycan_pos)
        else:
            if self.glycan_pos == 0:
                postions_for_glycan = str(self.native_glycan_pos)
            else:             
                postions_for_glycan = str(self.glycan_pos) + "," + str(self.native_glycan_pos)

        for element in root.iter('SimpleGlycosylateMover'):
                element.set ("positions", postions_for_glycan)

        
        # Write the changes back to the file
    
        GlycanTreeModeler_xml = "GlycanTreeModeler.xml"

        dst = GlycanTreeModeler_xml
        
        tree.write(dst)
        
        return GlycanTreeModeler_xml
    

    '''
    Helper function to run glycan tree modeler protocol

    '''
    
    def run_add_glycan_cmd(self):

        GlycanTreeModeler_xml = self.create_GlycanTreeModeler_xml()
        
        if self.glycan_pos == 0:
            designed_pdb = self.ag_pdb
        else:
            designed_pdb = self.run_fast_design_cmd()

        print("*****start of adding glycan")

        st = os.stat(Single_Glycan.GlycanTreeModeler_sh)
        os.chmod(Single_Glycan.GlycanTreeModeler_sh, st.st_mode | stat.S_IEXEC)

        GlycanTreeModeler_cmd = "./" + Single_Glycan.GlycanTreeModeler_sh
        #subprocess.run(f"bash my_script.sh '{my_variable}'", shell=True)
        cmd = f"bash GlycanTreeModeler.sh {designed_pdb} {str(self.n_struct)}"
        subprocess.run(cmd, shell = True)

        print("*****end of adding glycan")

        return "Glyc_" + designed_pdb
    

    def copy_files_only(source_dir, dest_dir):
        for item in os.listdir(source_dir):
            source_path = os.path.join(source_dir, item)
            dest_path = os.path.join(dest_dir, item)

            if os.path.isfile(source_path):
                shutil.copy2(source_path, dest_path)


         
    '''
    Function to run a complete process to add a single glcan to the wanted position

    Returns:
        0: process is completed successfully, when done, it goes back the parent_dir
        -1: failed

    '''

    def run_single_glycan(self):

        '''
        parent_dir = os.getcwd()
        
        print(parent_dir)
        template_files_dir = parent_dir + "/template_files"

        new_folder_name = "pos_" + str(self.glycan_pos)
    
        src_dir = template_files_dir
        dst_dir = new_folder_name

        # copy all example xml, options, resifle and protein structure related files(.pdb and fasta) from example_files dir to new dir
        # check if dst_dir exists
        if os.path.exists(dst_dir):
            shutil.rmtree(dst_dir) 
        shutil.copytree(src_dir, dst_dir)
        
        print(os.getcwd())
        # cd to the new dir
        os.chdir(new_folder_name)
        print(os.getcwd())
        '''
        print("-----------Start of adding glycan at position {pos} for {antigen}".format(pos = self.glycan_pos, antigen = self.ag_pdb))
        self.create_resfile_mod_seq()
        # run single glycan
        self.run_add_glycan_cmd()
        print("-----------End of adding glycan at position {pos} for {antigen}".format(pos = self.glycan_pos, antigen = self.ag_pdb))

        # go back to the parent dir
        #os.chdir(parent_dir)
        #print(os.getcwd())

        #cur_dir = os.getcwd()
        #double check the current location
        #if parent_dir == cur_dir:
            #return 0
        #else:
            #return -1


        


