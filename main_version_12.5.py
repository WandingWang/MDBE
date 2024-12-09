import os
import datetime
import logging
import shutil
import pandas as pd
import glob
import csv
#from concurrent.futures import ThreadPoolExecutor
import concurrent.futures
import time
import threading
import os
import shutil
import time
import multiprocessing

from FUNCTION import make_top_protein, fill_water_ions, energy_min, make_new_minim_nvt_npt, make_new_minim_config_samd, run_md
from FUNCTION import files_gmxmmpbsa, gmx_mmpbsa, Data_Analysis_Pre, Data_Analysis_Cal, clean_for_each_cycle

#PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.getcwd()

DATA_DIR = os.path.join(PROJECT_ROOT, "DATA")
VMD_DIR = os.path.join(PROJECT_ROOT, "VMD_FUNCTION")
FUNCTION_DIR = os.path.join(PROJECT_ROOT, "FUNCTION")
FORCE_FIELD_PATH = os.path.join(PROJECT_ROOT, "FORCE_FIELD")
MMPBSA_INFILE_PATH = os.path.join(PROJECT_ROOT, "gmx_mmpbsa_in")
# pdb file
protein_infile = "HLA_BiAB_protein_50ns" 
#protein_infile = "mtbind"
protein_file_path = os.path.join(DATA_DIR, f"{protein_infile}.pdb")

# MDP files
ions_mdp_file = "ions"
minim_mdp_file = "minim"
nvt_mdp_file = "NVT"
npt_mdp_file = "NPT"
samd_mdp_file = "SAMD"
md_mdp_file = "EngComp_ff14sb_custom"
only_protein_md_mdp_file = "Protein_EngComp_ff14sb_custom"

ions_mdp_path = os.path.join(DATA_DIR, f"{ions_mdp_file}.mdp")
minim_mdp_path = os.path.join(DATA_DIR, f"{minim_mdp_file}.mdp")
nvt_mdp_path = os.path.join(DATA_DIR, f"{nvt_mdp_file}.mdp")
npt_mdp_path = os.path.join(DATA_DIR, f"{npt_mdp_file}.mdp")
samd_mdp_path = os.path.join(DATA_DIR, f"{samd_mdp_file}.mdp")
md_mdp_path = os.path.join(DATA_DIR, f"{md_mdp_file}.mdp")
only_protein_md_mdp_path = os.path.join(DATA_DIR, f"{only_protein_md_mdp_file}.mdp")



def create_output_directory():
    
    current_dir = os.getcwd()
    
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir_name = f"output_{timestamp}"
    output_dir_path = os.path.join(current_dir, output_dir_name)
    
    os.mkdir(output_dir_path)
    print(f"Created directory: {output_dir_path}")

    os.chdir(output_dir_path)

    
    
    return output_dir_path

ROOT_OUTPUT = create_output_directory()

logging.basicConfig(
    filename = "OUTPUT.out",
    level = logging.INFO,
    format="%(asctime)s - %(levelname)s -%(message)s"
)

logging.info(f"PATH: {ROOT_OUTPUT}")

# create configuration folder
configuration_path = os.path.join(os.getcwd(),"configuration")
os.mkdir(configuration_path)
print(f"Create directory: {configuration_path}")
os.chdir(configuration_path)
logging.info("####Begin with configuration####")

# build folders
# cycle_num, MAKE IT 1 NOW, NEED TO CHANGE
cycle_num = 2
def build_folders(current_dir, cycle_num):
    # Create folder for each cycle
    folders  ={}
    
    for cycle_n in range (1,cycle_num + 1):
        folder_name = f"cycle{cycle_n}_MD"
        folder_path = os.path.join(current_dir, folder_name)
        os.makedirs(folder_path, exist_ok = True)
        folders[f"cycle{cycle_n}_MD"] = folder_path

    folders["repository"] = os.path.join(current_dir,"REPOSITORY")
    folders["TEMP_FILES_FOLDER"] = os.path.join(current_dir,"TEMP_FILES_FOLDER")
    folders["REMOVED_FILES_FOLDER"] = os.path.join(current_dir,"REMOVED_FILES_FOLDER")
    folders["results"] = os.path.join(current_dir,"RESULTS")

    for folder in folders.values():
        os.makedirs(folder,exist_ok = True)
  
    header = [
    "#RUNnumber", "DeltaG(kJ/mol)", "Coul(kJ/mol)", "vdW(kJ/mol)",
    "PolSol(kJ/mol)", "NpoSol(kJ/mol)", "ScoreFunct", "ScoreFunct2",
    "Canonica_AVG", "MedianDG", "DeltaG_2s", "dG_PotEn"]

    df = pd.DataFrame(columns=header)
    results_file_path = os.path.join(folders["results"], "MoleculesResults.dat")
    df.to_csv(results_file_path, sep='\t', index=False, header=True)

    return folders
    

current_dir = os.getcwd()
folders = build_folders(current_dir,cycle_num)

# generating a topology and build box
make_top_protein(protein_file_path, "amber99sb-ildn", "tip3p", "system", "topol")

# cp system.pdb {protein_infile}.pdb in current folder
source = os.path.join(current_dir, "system.pdb")
destination = os.path.join(current_dir, f"{protein_infile }.pdb")
try:
    shutil.copy(source,destination)
except Exception:
    print("Copy system.pdb failed.")

## OC2 DOESN'T NORMAL ONE, NEED TO CHANGE
def add_ter_to_pdb(pdb_file_name):
    
    temp_file_name = f"{pdb_file_name}_temp"  
    with open(pdb_file_name, 'r') as f:
        lines = f.readlines()

    new_lines = []  
    for i, line in enumerate(lines):
        new_lines.append(line)
        if "OC2" in line:
            if i + 1 >= len(lines) or not lines[i + 1].startswith("TER"):
                new_lines.append("TER\n")

    with open(temp_file_name, 'w') as f:
        f.writelines(new_lines)

    os.rename(temp_file_name, pdb_file_name)

add_ter_to_pdb(f"{protein_infile }.pdb")

def replace_his_(input_pdb, output_pdb):
    # input_pdb: the pdb in configuration/mutant x
    # output_pdb: make it in ROOT_OUTPUT folder
    with open(input_pdb, 'r') as infile:
        data = infile.read()
    data = data.replace("HISD", "HIS").replace("HISE", "HIS").replace("HISP", "HIS")

    with open(output_pdb, 'w') as outfile:
        outfile.write(data)

output_pdb = os.path.join(ROOT_OUTPUT, f"{protein_infile}.pdb")
replace_his_(f"{protein_infile}.pdb",output_pdb)

# Adding water and ions
fill_water_ions("system", "topol", ions_mdp_path)
# Energy Minimiization
energy_min(minim_mdp_path, "system_ions", "topol", "system_compl")

# Nvt and Npt
sequence = 0
make_new_minim_nvt_npt("system_compl_minim.gro", nvt_mdp_path, npt_mdp_path, "system_equil", 0)

# Move .cpt, .top, and .itp files to repository folder
for file_pattern in [f"{current_dir}/*.cpt", f"{current_dir}/*.top", f"{current_dir}/*.itp"]:
    for file in glob.glob(file_pattern):
        shutil.move(file, folders["repository"])

# Move specific files to repository folder
shutil.move(f"{current_dir}/{protein_infile}.pdb", folders["repository"])
shutil.move(f"{current_dir}/system_compl_minim.gro", folders["repository"])
shutil.move(f"{current_dir}/system_equil.gro", folders["repository"])


# Move temp* and *out files to removed files folder
for file in glob.glob("./*temp*.*") + glob.glob("./*.temp") + glob.glob("./*out"):
    shutil.move(file, folders["REMOVED_FILES_FOLDER"])

# Remove files with # in their name
for file in glob.glob("./#*"):
    os.remove(file)

'''
def MD_for_each_cycle(cycle_num,input_structure_file,samd_mdp_path,output_gro, sequence, md_mdp_path, tpr_file, trj_name):
    print("start MD in MD function")
    cycle_number = 1
    while cycle_number <= cycle_num:
        
        cycle_MD_path = folders[f"cycle{cycle_number}_MD"] 
        os.chdir(cycle_MD_path)
        shutil.copy(os.path.join(folders["repository"], "system_equil.gro"), "./")
        shutil.copy(os.path.join(folders["repository"], "topol.top"), "./")

        for itp_file in glob.glob(os.path.join(folders["repository"], "*rotein_chain_*.itp")):
            shutil.copy(itp_file, "./")

        for itp_file in glob.glob(os.path.join(folders["repository"], "posres_*.itp")):
            shutil.copy(itp_file, "./")

        for cpt_file in glob.glob(os.path.join(folders["repository"], "*NPT*.cpt")):
            shutil.copy(cpt_file, "./")

        fifo_path = os.path.join(cycle_MD_path, "mypipe")
        if not os.path.exists(fifo_path):
            os.mkfifo(fifo_path)

        #make_new_minim_config_samd("system_equil.gro", samd_mdp_path, "system_Compl_MDstart", 0)
        make_new_minim_config_samd(input_structure_file, samd_mdp_path, output_gro, sequence)
        #run_md(md_mdp_path,"system_Compl_MD", "traj_MD", 0, 1)
        run_md(md_mdp_path, tpr_file, trj_name, sequence, cycle_number)
        shutil.copy("system_Compl_MD.gro", f"LastFrame_cycle{cycle_number}.gro")
       
        with open (fifo_path, "w") as fifo:
            fifo.write("done")
        
        cycle_number += 1

def gmx_mmpbsa_for_each_cycle(cycle_num,only_protein_md_mdp_path,VMD_DIR,temp_files_folder, FORCE_FIELD_PATH, MMPBSA_INFILE_PATH, REMOVED_FILES_FOLDER, results_folder, repository_folder, current_conf_path):
    cycle_number = 1
    ConfName = f"cycle{cycle_number}"
    RootName = f"cycle{cycle_number}_BE"
    cycle_number_MD_FOLDER = folders[f"cycle{cycle_number}_MD"]
    print(f"Cycle Number: {cycle_number}")
    print(f"Configuration Name: {ConfName}")
    print(f"Root Name: {RootName}")
    print(f"MD Folder Path: {cycle_number_MD_FOLDER}")
    while cycle_number <= cycle_num:

        os.chdir(cycle_number_MD_FOLDER)
        waiting_time = 0
        max_waiting_time = 15
        fifo_path = os.path.join(cycle_number_MD_FOLDER, "mypipe")
        while not os.path.exists(fifo_path):
            time.sleep(60)
            waiting_time +=1
            if waiting_time == max_waiting_time:
                raise Exception(f"Unable to find {fifo_path} Waiting time={waiting_time}m")
        
        with open(fifo_path, "r") as fifo:
            fifo_output = fifo.read().strip()
            if fifo_output == "done":
                print("MD suscessuly")

                repository_pdb_file = os.path.join(repository_folder, f"{protein_infile}.pdb")
                #startingFrameGMXPBSA="2000"
                # make files for gmx_mmpbsa
                # files_gmxmmpbsa(starting_gro_file, repository_pdb_file, trj_file, tpr_file, top_file, mdp_name, root_name, conf_name, vmd_function_folder, temp_files_folder)

                files_gmxmmpbsa("system_Compl_MD", repository_pdb_file, "traj_MD", "system_Compl_MD", "topol", only_protein_md_mdp_path, RootName, ConfName, VMD_DIR, temp_files_folder, cycle_number)
                # get number of frames
                try:
                    with open("trj_check.out", "r") as file:
                        number_of_frames = next(
                            (line.split()[1] for line in file if line.startswith("Step")), None
                        )
                except FileNotFoundError:
                    print(f"Error: File trj_check.out not found.")
                    number_of_frames = None
                #conda_activate_path="/home/bio/ls/bin"
                conda_path = shutil.which("conda")
                conda_activate_path = os.path.dirname(conda_path)
                conda_gmxmmpbsa_name="gmxMMPBSA"
                forcefield="amber99sb-ildn"
                #FORCE_FIELD_PATH = "/home/bio/Desktop/jupyter_test/antibody_test/FORCE_FIELD"
                mmpbsa_inFILE="mmpbsa_LinearPB_amber99SB_ILDN.in"
                #MMPBSA_INFILE_PATH = "/home/bio/Desktop/jupyter_test/antibody_test/gmx_mmpbsa_in"
                np_value = 32
                #gmx_mmpbsa(1, conda_activate_path, conda_gmxmmpbsa_name, cycle_number_MD_FOLDER, ConfName, RootName, forcefield, FORCE_FIELD_PATH, 
                #             mmpbsa_inFILE, MMPBSA_INFILE_PATH , np_value, number_of_frames)
                gmx_mmpbsa(cycle_number, conda_activate_path, conda_gmxmmpbsa_name, cycle_number_MD_FOLDER, ConfName, RootName, forcefield, FORCE_FIELD_PATH, mmpbsa_inFILE, MMPBSA_INFILE_PATH, np_value, number_of_frames)
                # data analysis
                NUMframe = "all"
                Data_Analysis_Pre(cycle_number_MD_FOLDER, REMOVED_FILES_FOLDER, NUMframe)
                Data_Analysis_Cal(cycle_number, results_folder)
                # clean and move files
                clean_for_each_cycle(cycle_number, repository_folder, cycle_number_MD_FOLDER, RootName, REMOVED_FILES_FOLDER, current_conf_path)
                    # 更新周期号及相关变量
                #cycle_number += 1
                #conf_name = f"cycle{cycle_number}"
                #root_name = f"cycle{cycle_number}_BE"
                #cycle_number_md_folder = os.path.join(current_conf_path, f"cycle{cycle_number}_MD")
            else:
                raise Exception ('MD did not completed')
        cycle_number += 1



cycle_num =2
sequence = 0
MD_for_each_cycle(cycle_num,"system_equil.gro",samd_mdp_path,"system_Compl_MDstart", sequence, md_mdp_path, "system_Compl_MD", "traj_MD")
gmx_mmpbsa_for_each_cycle(cycle_num, only_protein_md_mdp_path,VMD_DIR,folders["TEMP_FILES_FOLDER"], FORCE_FIELD_PATH, MMPBSA_INFILE_PATH, folders["REMOVED_FILES_FOLDER"], folders["results"], folders["repository"], configuration_path)
'''

def MD_for_each_cycle(work_dir, cycle_number,input_structure_file,samd_mdp_path,output_gro, sequence, md_mdp_path, tpr_file, trj_name):
    print("start MD in MD function")
    #cycle_number = 1
    #while cycle_number <= cycle_num:
        
    cycle_MD_path = folders[f"cycle{cycle_number}_MD"] 
    os.chdir(work_dir)
    shutil.copy(os.path.join(folders["repository"], "system_equil.gro"), "./")
    shutil.copy(os.path.join(folders["repository"], "topol.top"), "./")

    for itp_file in glob.glob(os.path.join(folders["repository"], "*rotein_chain_*.itp")):
        shutil.copy(itp_file, "./")

    for itp_file in glob.glob(os.path.join(folders["repository"], "posres_*.itp")):
        shutil.copy(itp_file, "./")

    for cpt_file in glob.glob(os.path.join(folders["repository"], "*NPT*.cpt")):
        shutil.copy(cpt_file, "./")

    #make_new_minim_config_samd("system_equil.gro", samd_mdp_path, "system_Compl_MDstart", 0)
    make_new_minim_config_samd(input_structure_file, samd_mdp_path, output_gro, sequence)
    #run_md(md_mdp_path,"system_Compl_MD", "traj_MD", 0, 1)
    run_md(md_mdp_path, tpr_file, trj_name, sequence, cycle_number)
    shutil.copy("system_Compl_MD.gro", f"LastFrame_cycle{cycle_number}.gro")
    #cycle_number += 1


def gmx_mmpbsa_for_each_cycle(work_dir, cycle_number,only_protein_md_mdp_path,VMD_DIR,temp_files_folder, FORCE_FIELD_PATH, MMPBSA_INFILE_PATH, REMOVED_FILES_FOLDER, results_folder, repository_folder, current_conf_path):
    #cycle_number = 1
    #while cycle_number <= cycle_num:
    ConfName = f"cycle{cycle_number}"
    RootName = f"cycle{cycle_number}_BE"
    cycle_number_MD_FOLDER = folders[f"cycle{cycle_number}_MD"]
    # 输出相关信息
    print(f"Cycle Number: {cycle_number}")
    print(f"Configuration Name: {ConfName}")
    print(f"Root Name: {RootName}")
    print(f"MD Folder Path: {cycle_number_MD_FOLDER}")
    #os.chdir(cycle_number_MD_FOLDER )
    os.chdir(work_dir)
        
    repository_pdb_file = os.path.join(repository_folder, f"{protein_infile}.pdb")
    #startingFrameGMXPBSA="2000"
    # make files for gmx_mmpbsa
    # files_gmxmmpbsa(starting_gro_file, repository_pdb_file, trj_file, tpr_file, top_file, mdp_name, root_name, conf_name, vmd_function_folder, temp_files_folder)

    files_gmxmmpbsa("system_Compl_MD", repository_pdb_file, "traj_MD", "system_Compl_MD", "topol", only_protein_md_mdp_path, RootName, ConfName, VMD_DIR, temp_files_folder, cycle_number)
    # get number of frames
    try:
        with open("trj_check.out", "r") as file:
            number_of_frames = next(
                (line.split()[1] for line in file if line.startswith("Step")), None
            )
    except FileNotFoundError:
        print(f"Error: File trj_check.out not found.")
        number_of_frames = None
    #conda_activate_path="/home/bio/ls/bin"
    conda_path = shutil.which("conda")
    conda_activate_path = os.path.dirname(conda_path)
    conda_gmxmmpbsa_name="gmxMMPBSA"
    forcefield="amber99sb-ildn"
    #FORCE_FIELD_PATH = "/home/bio/Desktop/jupyter_test/antibody_test/FORCE_FIELD"
    mmpbsa_inFILE="mmpbsa_LinearPB_amber99SB_ILDN.in"
    #MMPBSA_INFILE_PATH = "/home/bio/Desktop/jupyter_test/antibody_test/gmx_mmpbsa_in"
    np_value = 32
    #gmx_mmpbsa(1, conda_activate_path, conda_gmxmmpbsa_name, cycle_number_MD_FOLDER, ConfName, RootName, forcefield, FORCE_FIELD_PATH, 
    #             mmpbsa_inFILE, MMPBSA_INFILE_PATH , np_value, number_of_frames)
    gmx_mmpbsa(cycle_number, conda_activate_path, conda_gmxmmpbsa_name, cycle_number_MD_FOLDER, ConfName, RootName, forcefield, FORCE_FIELD_PATH, mmpbsa_inFILE, MMPBSA_INFILE_PATH, np_value, number_of_frames)
    # data analysis
    NUMframe = "all"
    Data_Analysis_Pre(cycle_number_MD_FOLDER, REMOVED_FILES_FOLDER, NUMframe)
    Data_Analysis_Cal(cycle_number, results_folder)
    # clean and move files
    clean_for_each_cycle(cycle_number, repository_folder, cycle_number_MD_FOLDER, RootName, REMOVED_FILES_FOLDER, current_conf_path)
        # 更新周期号及相关变量
    #cycle_number += 1
    #conf_name = f"cycle{cycle_number}"
    #root_name = f"cycle{cycle_number}_BE"
    #cycle_number_md_folder = os.path.join(current_conf_path, f"cycle{cycle_number}_MD")
    
 


cycle_num =2
sequence = 0
md_args = ("system_equil.gro",samd_mdp_path,"system_Compl_MDstart", sequence, md_mdp_path, "system_Compl_MD", "traj_MD")
gmx_args = (only_protein_md_mdp_path,VMD_DIR,folders["TEMP_FILES_FOLDER"], FORCE_FIELD_PATH, MMPBSA_INFILE_PATH, folders["REMOVED_FILES_FOLDER"], folders["results"], folders["repository"], configuration_path)




def run_cycle(cycle_number, cycle_num, md_args, gmx_args):
    """
    处理当前 cycle 的 gmx_mmpbsa 和下一个 cycle 的 MD。
    """
    process = []

    # 启动当前 cycle 的 gmx_mmpbsa
    gmx_process = multiprocessing.Process(target=gmx_mmpbsa_for_each_cycle, args=(folders[f"cycle{cycle_number}_MD"],cycle_number, *gmx_args))
    process.append(gmx_process)
    gmx_process.start()

    # 如果还有下一个 cycle，启动下一个 cycle 的 MD
    if cycle_number < cycle_num:
        next_md_process = multiprocessing.Process(target=MD_for_each_cycle, args=(folders[f"cycle{cycle_number+1}_MD"], cycle_number + 1, *md_args))
        process.append(next_md_process)
        next_md_process.start()

    # 等待所有线程完成
    for p in process:
        p.join()



# 首先启动第一个 cycle 的 MD
MD_for_each_cycle(folders["cycle1_MD"],1, *md_args)

# 然后依次处理每个 cycle 的 gmx_mmpbsa 和下一个 MD
for cycle_number in range(1, cycle_num + 1):
    run_cycle(cycle_number, cycle_num, md_args, gmx_args)




