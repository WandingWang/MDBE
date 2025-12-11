########################################## IMPORT MODULES ######################################
import os
import datetime
import logging
import shutil
import pandas as pd
import glob
import multiprocessing
import re
import subprocess
import numpy as np
#import math
#import random
import yaml

from FUNCTION import make_top_protein, fill_water_ions, energy_min, run_md, make_new_minim_nvt_npt
#  make_new_minim_nvt_npt
from FUNCTION import files_gmxmmpbsa, gmx_mmpbsa, Data_Analysis_Pre, Data_Analysis_Cal, clean_for_each_cycle, GRO_to_PDB
from FUNCTION import Data_Analysis_Cal_child

############################################### FUNCTION DEFINATION #####################################
def load_config(config_file):
    with open(config_file, "r") as file:
        config = yaml.safe_load(file)
    return config

def get_version(command):
    try:
        results = subprocess.run([command,"--version"],capture_output=True, text=True, check=True)
        return results.stdout.strip()
    except subprocess.CalledProcessError:
        return None

def which_program(command):
    try:
        results = subprocess.run(["which",command],capture_output=True, text=True, check=True)
        return results.stdout.strip()
    except subprocess.CalledProcessError:
        return None

def check_conda_env(env):
    try:
        check_env = subprocess.run(['conda', 'info', '--envs'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
        if env in check_env.stdout:
            return True
        else:
            return False
    except subprocess.CalledProcessError:
        return False
'''    
def create_output_directory():
    
    current_dir = os.getcwd()
    
    timestamp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')
    output_dir_name = f"output_{timestamp}"
    output_dir_path = os.path.join(current_dir, output_dir_name)
    
    os.mkdir(output_dir_path)
    print(f"Created directory: {output_dir_path}")

    os.chdir(output_dir_path)

    
    
    return output_dir_path
'''

def create_output_directory():
    current_dir = os.getcwd()
    
    # Get the current date formatted as MMDD
    date_str = datetime.datetime.now().strftime('%m%d')
    base_name = f"output_{date_str}_"
    index = 1
    
    # Generate a new directory name that doesn't already exist
    while True:
        output_dir_name = f"{base_name}{index}"
        output_dir_path = os.path.join(current_dir, output_dir_name)
        
        if not os.path.exists(output_dir_path):
            os.mkdir(output_dir_path)
            print(f"Created directory: {output_dir_path}")
            break
        
        index += 1

    os.chdir(output_dir_path)
    return output_dir_path

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
    '''
    header = [
    "#RUNnumber", "DeltaG(kcal/mol)", "Coul(kcal/mol)", "vdW(kcal/mol)",
    "PolSol(kcal/mol)", "NpoSol(kcal/mol)", "ScoreFunct", "ScoreFunct2",
    "Canonica_AVG", "MedianDG", "DeltaG_2s", "dG_PotEn"]
    '''
    
    header = [
    "#RUNnumber", "DeltaG(kcal/mol)", "Coul(kcal/mol)", "vdW(kcal/mol)",
    "PolSol(kcal/mol)", "NpoSol(kcal/mol)", "ScoreFunct", "ScoreFunct2",
    "Canonica_AVG", "MedianDG", "DeltaG_2s"]

    df = pd.DataFrame(columns=header)
    results_file_path = os.path.join(folders["results"], "MoleculesResults.dat")
    df.to_csv(results_file_path, sep='\t', index=False, header=True)

    return folders

def add_ter_to_pdb(pdb_file_name):
    temp_file_name = f"{pdb_file_name}_temp"  

    with open(pdb_file_name, 'r') as f:
        lines = f.readlines()

    new_lines = []
    prev_chain_id = None  
    num_lines = len(lines)

    for i, line in enumerate(lines):
  
        if not line.startswith("ATOM"):
            new_lines.append(line)
            continue


        current_chain_id = line[21]  # get chain ID (the 22rd column)

        # add
        new_lines.append(line)


        if prev_chain_id is not None:

            if (i == num_lines - 1 or 
                (lines[i + 1].startswith("ATOM") and lines[i + 1][21] != current_chain_id) or 
                not lines[i + 1].startswith("ATOM")):
                new_lines.append("TER\n")

 
        prev_chain_id = current_chain_id


    with open(temp_file_name, 'w') as f:
        f.writelines(new_lines)

  
    os.rename(temp_file_name, pdb_file_name)

def replace_his_residues_flexible(input_pdb, output_pdb):
    with open(input_pdb, "r", encoding="utf-8") as infile, open(output_pdb, "w", encoding="utf-8") as outfile:
        for line in infile:
            if line.startswith("ATOM"):
                match = re.match(r"(.{17})(HISE|HISD|HISP)(.*)", line)
                if match:
                    line = f"{match.group(1)}{'HIS':<4}{match.group(3)}"
            
            outfile.write(line.rstrip('\n') + '\n')

def MD_for_each_cycle(work_dir, cycle_number,sequence, md_mdp_path, tpr_file, trj_name, gmx_path):
    print("start MD in MD function")
    #cycle_number = 1
    #while cycle_number <= cycle_num:
        
    #cycle_MD_path = folders[f"cycle{cycle_number}_MD"] 
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
    #make_new_minim_config_samd(input_structure_file, samd_mdp_path, output_gro, sequence)
    #run_md(md_mdp_path,"system_Compl_MD", "traj_MD", 0, 1)
    run_md(md_mdp_path, tpr_file, trj_name, sequence, cycle_number, gmx_path)
    shutil.copy("system_Compl_MD.gro", f"LastFrame_cycle{cycle_number}.gro")
    #cycle_number += 1


def gmx_mmpbsa_for_each_cycle(work_dir, cycle_number,only_protein_md_mdp_path,VMD_DIR,temp_files_folder, FORCE_FIELD_PATH, MMPBSA_INFILE_PATH, REMOVED_FILES_FOLDER, results_folder, repository_folder, current_conf_path):
    #cycle_number = 1
    #while cycle_number <= cycle_num:
    ConfName = f"cycle{cycle_number}"
    RootName = f"cycle{cycle_number}_BE"
    cycle_number_MD_FOLDER = folders[f"cycle{cycle_number}_MD"]

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

    files_gmxmmpbsa("system_Compl_MD", repository_pdb_file, "traj_MD", "system_Compl_MD", "topol", only_protein_md_mdp_path, RootName, ConfName, VMD_DIR, temp_files_folder, cycle_number, startingFrameGMXPBSA, receptorFRAG, ABchains,gmx_path, VMD_path )
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

    #conda_gmxmmpbsa_name="gmxMMPBSA"
    #forcefield="amber99sb-ildn"
    #FORCE_FIELD_PATH = "/home/bio/Desktop/jupyter_test/antibody_test/FORCE_FIELD"
    mmpbsa_inFILE="mmpbsa_LinearPB_amber99SB_ILDN.in"
    #MMPBSA_INFILE_PATH = "/home/bio/Desktop/jupyter_test/antibody_test/gmx_mmpbsa_in"
    np_value = config['run']['num_processors']
    #gmx_mmpbsa(1, conda_activate_path, conda_gmxmmpbsa_name, cycle_number_MD_FOLDER, ConfName, RootName, forcefield, FORCE_FIELD_PATH, 
    #             mmpbsa_inFILE, MMPBSA_INFILE_PATH , np_value, number_of_frames)
    gmx_mmpbsa(cycle_number, conda_actiavte_path, conda_gmxmmpbsa_name, cycle_number_MD_FOLDER, ConfName, RootName, mmpbsa_inFILE, MMPBSA_INFILE_PATH, np_value, number_of_frames)
    # data analysis
    NUMframe = "all"
    Data_Analysis_Pre(cycle_number_MD_FOLDER, REMOVED_FILES_FOLDER, NUMframe)
    Data_Analysis_Cal(cycle_number, results_folder)
    # clean and move files
    clean_for_each_cycle(cycle_number, repository_folder, cycle_number_MD_FOLDER, RootName, REMOVED_FILES_FOLDER, current_conf_path)
    #cycle_number += 1
    #conf_name = f"cycle{cycle_number}"
    #root_name = f"cycle{cycle_number}_BE"
    #cycle_number_md_folder = os.path.join(current_conf_path, f"cycle{cycle_number}_MD")

def run_cycle(cycle_number, cycle_num, md_args, gmx_args):
    """
    deal with gmx_MMPBSA for current cycle and MD for the next cycle
    """
    process = []

    # current cycle gmx_mmpbsa
    gmx_process = multiprocessing.Process(target=gmx_mmpbsa_for_each_cycle, args=(folders[f"cycle{cycle_number}_MD"],cycle_number, *gmx_args))
    process.append(gmx_process)
    gmx_process.start()

    # if we have next cycle，run MD for the next cycle 
    if cycle_number < cycle_num:
        next_md_process = multiprocessing.Process(target=MD_for_each_cycle, args=(folders[f"cycle{cycle_number+1}_MD"], cycle_number + 1, *md_args))
        process.append(next_md_process)
        next_md_process.start()

    # all processes finished
    for p in process:
        p.join()


############################################# LOAD FILES ###########################################

# Command-line argument parsing
parser = argparse.ArgumentParser(description='MDBE Pipeline')
parser.add_argument('-i', '--infile', type=str, required=True,
                   help='Input YAML configuration file')

args = parser.parse_args()

# Load configuration from the provided YAML file
config_file = args.infile
config = load_config(config_file )


#PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.getcwd()

DATA_DIR = os.path.join(PROJECT_ROOT, "DATA")
VMD_DIR = os.path.join(PROJECT_ROOT, "VMD_FUNCTION")
FUNCTION_DIR = os.path.join(PROJECT_ROOT, "FUNCTION")
FORCE_FIELD_PATH = os.path.join(PROJECT_ROOT, "FORCE_FIELD")
MMPBSA_INFILE_PATH = os.path.join(PROJECT_ROOT, "gmx_mmpbsa_in")
# pdb file
#protein_infile = "HLA_BiAB_protein_50ns" 
#protein_infile = "mtbind"
#protein_infile = "antibody_zixuan"

#protein_file_path = os.path.join(DATA_DIR, f"{protein_infile}.pdb")

#make_mutation_modeller_py = os.path.join(FUNCTION_DIR,"MakeNewMutant_Modeller.py") 
# MDP files
ions_mdp_file = "ions"
minim_mdp_file = "minim"
md_mdp_file = "EngComp_ff14sb_custom"
only_protein_md_mdp_file = "Protein_EngComp_ff14sb_custom"
nvt_mdp_file = "NVT"
npt_mdp_file = "NPT"
ions_mdp_path = os.path.join(DATA_DIR, f"{ions_mdp_file}.mdp")
minim_mdp_path = os.path.join(DATA_DIR, f"{minim_mdp_file}.mdp")
nvt_mdp_path = os.path.join(DATA_DIR, f"{nvt_mdp_file}.mdp")
npt_mdp_path = os.path.join(DATA_DIR, f"{npt_mdp_file}.mdp")
md_mdp_path = os.path.join(DATA_DIR, f"{md_mdp_file}.mdp")
only_protein_md_mdp_path = os.path.join(DATA_DIR, f"{only_protein_md_mdp_file}.mdp")

ROOT_OUTPUT = create_output_directory()

logging.basicConfig(
    filename = "OUTPUT.out",
    level = logging.INFO,
    format="%(asctime)s - %(levelname)s -%(message)s"
)

logging.info(f"ROOT FOLDER PATH: {ROOT_OUTPUT}")

####################################### CHECK SYSTEM #############################################
#config_file = os.path.join(PROJECT_ROOT, 'infile.yaml')
#config = load_config(config_file )
conda_actiavte_path = config['Basic_setting']['conda_activate_script_path']
VMD_path = config['Basic_setting']['VMD_path']
gmx_path = config['Basic_setting']['GROMACS_executable_path']

conda_gmxmmpbsa_name = config['Basic_setting']['conda_gmx_MMPBSA_name']
#conda_modeller_name = config['Basic_setting']['conda_Modeller_name']


# check conda
if not os.path.isfile(conda_actiavte_path):
    logging.error(f"ERROR: cannot find conda activate scrpt path as {conda_actiavte_path}")
else:
    logging.info(f"conda activate path --> {conda_actiavte_path} version: {get_version('conda')}")
    
# check python
python_version = get_version("python")
logging.info(f"Python --> {which_program('python')} version: {python_version}")

# check VMD
if VMD_path:
    if os.path.isfile(VMD_path) and os.access(VMD_path, os.X_OK):
        logging.info(f"VMD path --> {VMD_path} ")
    else:
        logging.error(f"ERROR: cannot find VMD path as {VMD_path}")
else:
    logging.info(f"VMD path --> {which_program('vmd')} ")

# check gromacs
if gmx_path:
    if os.path.isfile(gmx_path) and os.access(gmx_path, os.X_OK):
        logging.info(f"gmx path --> {gmx_path} ")
    else:
        logging.error(f"ERROR: cannot find gmx path as {gmx_path}")
else:
    logging.info(f"gmx path --> {which_program('gmx')} ")

# check gmx_mmpbsa

# check parameters
receptorFRAG = str(config['gmx_mmpbsa']['receptorFRAG'])
ABchains = str(config['gmx_mmpbsa']['ABchains'])
startingFrameGMXPBSA = config['gmx_mmpbsa']['startingFrameGMXPBSA']
#protein_infile = config['input_files']['structure_infile_name']
protein_file_path = config['input_files']['structure_file_path']
if not os.path.exists(protein_file_path):
    raise ValueError(f"No path {protein_file_path}")
else:
    logging.info(f"The structure file is {protein_file_path}")
protein_infile = os.path.basename(protein_file_path)
protein_infile, _ =os.path.splitext(protein_infile)


#max_mutant = config['modeller']['max_mutant']
#cycle_num = config['run']['cycle_num'] # the run cycle numbers for each configuration  Default:10

cycle_num = 2



######################################### MAIN PROCESS ######################################



try:
    os.chdir(ROOT_OUTPUT)
except OSError:
    logging.error(f"Cannot enter {ROOT_OUTPUT} folder")
    exit()


sequence = 0
# create configuration folder
configuration_path = os.path.join(ROOT_OUTPUT,"configuration")
os.mkdir(configuration_path)
current_path_store = configuration_path
        
print(f"Create directory: {configuration_path}")
os.chdir(configuration_path)
logging.info(f"#### Begin with configuration ####")
logging.info(f"PATH : {configuration_path}")


current_dir = os.getcwd()
folders = build_folders(current_dir,cycle_num)

# generating a topology and build box
make_top_protein(protein_file_path, "amber99sb-ildn", "tip3p", "system", "topol", gmx_path)

# cp system.pdb {protein_infile}.pdb in current folder
source = os.path.join(current_dir, "system.pdb")
destination = os.path.join(current_dir, f"{protein_infile }.pdb")
try:
    shutil.copy(source,destination)
except Exception:
    print("Copy system.pdb failed.")

add_ter_to_pdb(f"{protein_infile }.pdb")
output_pdb = os.path.join(ROOT_OUTPUT, f"{protein_infile}.pdb")
replace_his_residues_flexible(f"{protein_infile}.pdb",output_pdb)
# Adding water and ions
fill_water_ions("system", "topol", ions_mdp_path, gmx_path)
# Energy Minimiization
energy_min(minim_mdp_path, "system_ions", "topol", "system_compl",gmx_path)

    
# Nvt and Npt
make_new_minim_nvt_npt("system_compl_minim.gro", nvt_mdp_path, npt_mdp_path, "system_equil", 0, gmx_path)
#shutil.copy("system_compl_minim.gro",os.path.join(folders["repository"],"system_equil.gro"))

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

md_args = (sequence, md_mdp_path, "system_Compl_MD", "traj_MD", f"{gmx_path}")
gmx_args = (only_protein_md_mdp_path,VMD_DIR,folders["TEMP_FILES_FOLDER"], FORCE_FIELD_PATH, MMPBSA_INFILE_PATH, folders["REMOVED_FILES_FOLDER"], folders["results"], folders["repository"], current_path_store)
# 1st cycle MD
MD_for_each_cycle(folders["cycle1_MD"],1, *md_args)

# each cycle: gmx_mmpbsa and next MD
for cycle_number in range(1, cycle_num + 1):
    run_cycle(cycle_number, cycle_num, md_args, gmx_args)

    
last_cycle_MD_FOLDER = os.path.join(folders["repository"],f"cycle{cycle_number}_MD")
last_cycle_gro = os.path.join(last_cycle_MD_FOLDER,f"LastFrame_cycle{cycle_number}.gro")
shutil.copy(last_cycle_gro, os.path.join(folders["repository"],f"LastFrame_cycle{cycle_number}.gro"))

    
os.chdir(current_path_store)
    
repository_pdb_file = os.path.join(folders["repository"], f"{protein_infile}.pdb")

pathGRO = folders["repository"]
fileNameGRO = f"LastFrame_cycle{cycle_number}"
pathPDB = os.path.dirname(repository_pdb_file)
pdb_name_with_extension = os.path.basename(repository_pdb_file) #xxxx.pdb
pdb_name_without_extension = os.path.splitext(pdb_name_with_extension)[0] #xxxx
fileNamePDB = pdb_name_without_extension
FileNamePDB_OUT = f"LastFrame_cycle{cycle_number}"
GRO_to_PDB(pathGRO, fileNameGRO, pathPDB, fileNamePDB, FileNamePDB_OUT, VMD_DIR, folders["TEMP_FILES_FOLDER"], VMD_path)
last_cycle_pdb = os.path.join(folders["repository"], f"LastFrame_cycle{cycle_number}.pdb")
add_ter_to_pdb(last_cycle_pdb)        
output_last_cycle_pdb = os.path.join(ROOT_OUTPUT, f"Mutant{sequence}_cycle{cycle_number}_LastFrameMD.pdb")
replace_his_residues_flexible(last_cycle_pdb,output_last_cycle_pdb)

protein_infile = f"Mutant{sequence}_cycle{cycle_number}_LastFrameMD"
logging.info("Making the average of the cycles results.")
    
all_cycle_data = "All_cycle_data.out"
MoleculesResults_data = os.path.join(folders["results"], "MoleculesResults.dat")
data_analysis_temp = "DataAnalysis_temp.csv"

if os.path.exists(all_cycle_data):
    os.remove(all_cycle_data)

flag_header = True

with open(MoleculesResults_data, 'r') as infile, open(all_cycle_data, 'w') as outfile:
    for line in infile:
        if flag_header == True:
            # head row
            outfile.write(f"#{'configNum':<10} \t{line}")
            flag_header = False
        else:
            # data row
            outfile.write(f"{'avg':<10} \t{line}")
# get the data from the second line
df = pd.read_csv(all_cycle_data, sep='\t')  
df_filtered = df.iloc[:, 1:]  # get the data from the second column

# save to csv
df_filtered.to_csv(data_analysis_temp, sep = '\t', index=False,header = False)
#Data_Analysis_Signal = False
Data_Analysis_Cal_child(data_analysis_temp, "AllData.temp", False)

frame_count = 0
# get AVG and STD to AllData.out
with open("AllData.temp", 'r') as temp_file, open(all_cycle_data, 'a') as outfile:
    for line in temp_file:

        if line.startswith("#frame"):
            frame_count +=1
            if frame_count ==2:
                outfile.write(line)
        if line.startswith("#AVG") or line.startswith("#STD"):
            outfile.write(line)

frame = []
avg = []
std = []
with open(all_cycle_data, 'r') as infile:
    for line in infile:
        if line.startswith("#frame"):
            frame = line.strip().split()[1:]
        elif line.startswith("#AVG"):
            avg = line.strip().split()[1:]
        elif line.startswith("#STD"):
            std = line.strip().split()[1:]
if frame and avg and std:
    output_lines = ["Results for Configuation"]
    output_lines += [f"{frame[i]}: {avg[i]} +- {std[i]} kcal/mol"
                    for i in range(len(frame))
                    ]
    logging.info("\n".join(output_lines))
shutil.move("AllData.temp", folders["REMOVED_FILES_FOLDER"])
shutil.move(data_analysis_temp, folders["REMOVED_FILES_FOLDER"])
    

logging.info(f"Finished with Configuration")
        
os.chdir(ROOT_OUTPUT)

logging.info("ALL DONE.")
    
