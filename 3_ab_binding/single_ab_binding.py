import xml.etree.ElementTree as ET
import shutil
import sys
import subprocess
from subprocess import run
from shlex import split
import os
import stat

class Single_ab_binding:
    example_AddChain_xml = "template_AddChain_FastRelax_Multi.xml"
    addChain_sh = "addChain.sh"
    docking_analysis_xml = "template_docking_analysis.xml"
    cal_binding_sh = "cal_binding_energy.sh"
    clean_sh = "clean.sh"

    def __init__(self, antigen_list: str, ab:str, ab_heavy_pdb:str, ab_light_pdb:str):

        self.antigen_list = antigen_list
        self.ab = ab
        self.ab_heavy_pdb = ab_heavy_pdb
        self.ab_light_pdb= ab_light_pdb



    def create_AddChain_xml(self):
        tree = ET.parse(Single_ab_binding.example_AddChain_xml)
        print(tree)
        root = tree.getroot()
        print(root)

        #heavy_chain_pdb = self.ab_heavy[self.ab]
        #light_chain_pdb = self.ab_light[self.ab]
        # Modify an element
        for addChain in root.iter('AddChain'):
            if addChain.get('name') == 'Add_H':
                addChain.set("file_name", self.ab_heavy_pdb)
            if addChain.get('name') == 'Add_L':
                addChain.set("file_name", self.ab_light_pdb)
                  
        # Write the changes back to the file
               
        addChain_xml_name = "AddChain_FastRelax_Multi_{ab}.xml".format(ab = self.ab)

        dst = addChain_xml_name
        
        tree.write(dst)

        return addChain_xml_name
    
    
    def run_addChain_cmd(self):

        addChain_xml_name = self.create_AddChain_xml()

        print("*****start of adding antibody chains")

        st = os.stat(Single_ab_binding.addChain_sh)
        os.chmod(Single_ab_binding.addChain_sh, st.st_mode | stat.S_IEXEC)

        addChain_cmd = "./" + Single_ab_binding.addChain_sh
        #subprocess.call(addChain_cmd, shell = True)
        #subprocess.run(f"bash my_script.sh '{my_variable}'", shell=True)

        cmd = f"bash addChain.sh {self.antigen_list} {addChain_xml_name} {self.ab}"
        print(cmd)
        subprocess.run(cmd, shell = True)

        print("*****end of adding antibody chains")

        #return 0
        self.run_docking_analysis_cmd()
        return 0
    
    def creat_docking_list(self):
        out_file_name = f"{self.ab}_docking_list"

        with open(self.antigen_list, 'r') as infile, open(out_file_name, 'w') as outfile:
            for line in infile:
                pdb_suffix_for_clean = line.strip().replace(".pdb", "_0001.pdb")
                # clean the add_rlx pdb
                pdb_for_clean = f"Add_Rlx_{self.ab}{pdb_suffix_for_clean}"
                self.run_clean_cmd(pdb_for_clean)
                
                pdb_suffix_for_docking = line.strip().replace(".pdb", "_0001_ABC.pdb")
                new_line = f"Add_Rlx_{self.ab}{pdb_suffix_for_docking}"
                outfile.write(new_line + "\n")

        infile.close()
        outfile.close()
        return out_file_name

    def run_clean_cmd(self, pdb_for_clean):
        print(f"**************clean {pdb_for_clean}")

        st = os.stat(Single_ab_binding.clean_sh)
        os.chmod(Single_ab_binding.clean_sh, st.st_mode | stat.S_IEXEC)

        cmd = f"bash clean.sh {pdb_for_clean}"
        print(cmd)
        subprocess.run(cmd, shell = True)

        print(f"*****end of cleaning {pdb_for_clean}")
    
    def run_docking_analysis_cmd(self):
        print("*********start calculating binding enengy")

        st = os.stat(Single_ab_binding.cal_binding_sh)
        os.chmod(Single_ab_binding.cal_binding_sh, st.st_mode | stat.S_IEXEC)

        docking_list = self.creat_docking_list()

        cmd = f"bash cal_binding_energy.sh {docking_list} {self.ab}"
        print(cmd)
        subprocess.run(cmd, shell = True)

        print("*****end of calculating binding energy")

        return 0
    


    
