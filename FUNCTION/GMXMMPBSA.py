import os
import subprocess
import sys
import shutil

def check_file(file_path):
    """检查指定的文件是否存在"""
    if not os.path.isfile(file_path):
        print(f"ERROR: File {file_path} not found!")
        return False
    return True

#GRO_to_PDB(pathGRO, fileNameGRO, pathPDB, fileNamePDB, FileNamePDB_OUT, vmd_function_folder, temp_files_folder)
def GRO_to_PDB(pathGRO, fileNameGRO, pathPDB, fileNamePDB, FileNamePDB_OUT, vmd_function_folder, temp_files_folder):
    """
    Converts a GRO file to a PDB file using VMD.

    Parameters:
        pathGRO (str): Path to the input GRO file.
        fileNameGRO (str): Name of the GRO file.
        pathPDB (str): Path to the output PDB file. 
        fileNamePDB (str): Name of the output PDB file.
        temp_files_folder (str): Temporary files folder.
        vmd_executable (str): Path to VMD executable.
    """
    
    vmd = shutil.which("vmd")
    if not vmd:
        raise FileNotFoundError("VMD not found. please ensure VMD is installed and added to PATH")
    # Paths and filenames
    vmd_script = os.path.join(temp_files_folder, "vmd_GRO_to_PDB.tcl")
    log_file = os.path.join(temp_files_folder, "vmd_GRO_to_PDB.out")
    function_tcl = os.path.join(vmd_function_folder, "VMD_function_GRO_to_PDB.tcl")

    # Ensure paths and required files exist
    if not os.path.exists(temp_files_folder):
        os.makedirs(temp_files_folder)
    if not os.path.isfile(function_tcl):
        raise FileNotFoundError(f"Function script missing: {function_tcl}")

    if not FileNamePDB_OUT:
        FileNamePDB_OUT = fileNameGRO

    # Write VMD TCL script
    tcl_content = f"""
    variable _pathGRO "{pathGRO}"
    variable _FileNameGRO "{fileNameGRO}"
    variable _pathPDB "{pathPDB}"
    variable _FileNamePDB "{fileNamePDB}"
    variable _FileNamePDB_OUT "{FileNamePDB_OUT}"
    """
    with open(vmd_script, 'w') as f:
        f.write(tcl_content)
        f.write(open(function_tcl).read())

    # Execute VMD
    result = subprocess.run(
        [vmd, "-dispdev", "none", "-e", vmd_script],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Log results
    with open(log_file, 'w') as log:
        log.write(result.stdout + result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"VMD execution failed. Check log: {log_file}")

    print(f"GRO to PDB completed successfully: {fileNamePDB}")



def make_ndx_string_gmxpbsa(receptor_frag, ab_chains):
    """
    Generate the make_ndx string for GROMACS PBSA analysis.
    
    receptor_frag: Number of fragments in the receptor (protein part).
    ab_chains: Number of antibody chains (ligand part).
    return: A formatted make_ndx string for GROMACS.
    """
    # Initialize the make_ndx string with the initial commands:
    # 1. `keep 1`: Keep the first group (usually the whole structure).
    # 2. `splitch 0`: Split chains in the structure (0 is the default group for all atoms).
    make_ndx_string = "keep 1\nsplitch 0\n"
    
    # Handle receptor fragments
    if int(receptor_frag) > 1:
        # If there are multiple fragments in the receptor, create a group by combining them.
        # Start with the first fragment.
        flag = "1"
        # Delete the initial combined group after creating a new group for each fragment.
        flag2 = "del 1\n"
        
        # Loop through the receptor fragments, starting from the second fragment.
        for i in range(2, int(receptor_frag) + 1):
            # Combine fragments using the `|` operator in GROMACS make_ndx syntax.
            flag += f"|{i}"
            # Add a deletion command for each fragment.
            flag2 += "del 1\n"
        
        # Append the constructed commands to the make_ndx string.
        make_ndx_string += f"{flag}\n{flag2}"
    else:
        # If there's only one receptor fragment, directly select the combined group 0 and 1.
        make_ndx_string += "0&1\ndel 1\n"
    
    # Handle antibody chains
    if int(ab_chains) > 1:
        # If there are multiple antibody chains, create a group by combining them.
        flag = "1"
        # Delete the initial group for antibodies after creating new groups for each chain.
        flag2 = "del 1\n"
        
        # Loop through the antibody chains, starting from the second chain.
        for i in range(2, int(ab_chains) + 1):
            # Combine chains using the `|` operator.
            flag += f"|{i}"
            # Add a deletion command for each chain.
            flag2 += "del 1\n"
        
        # Append the constructed commands to the make_ndx string.
        # Finally, create a combined group for receptor and antibody chains.
        make_ndx_string += f"{flag}\n{flag2}del 0\n0|1\n"
    else:
        # If there's only one antibody chain, directly select and combine group 0 and 1.
        make_ndx_string += "0&1\ndel 0\ndel 0\n0|1\n"
    
    # Add naming commands to label the created groups:
    # `name 0 receptor`: Name group 0 as "receptor".
    # `name 1 ligand`: Name group 1 as "ligand".
    # `name 2 complex`: Name the combined group as "complex".
    # The final `q` command exits the make_ndx tool.
    make_ndx_string += "name 0 receptor\nname 1 ligand\nname 2 complex\n\nq\n"
    
    # Return the complete make_ndx string.
    return make_ndx_string


def run_command(cmd, input_data):
    """runnning command with inteaction"""
    try:
        process = subprocess.Popen(
            cmd,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        stdout, stderr = process.communicate(input=input_data)
        if process.returncode != 0:
            raise RuntimeError(f"Command failed: {stderr}")
        return stdout
    except Exception as e:
        print(f"Error executing command: {cmd}\n{e}")
        raise

def remove_pbc(trj_file, tpr_file, startingFrameGMXPBSA, root_name, conf_name):
    """remove PBC"""
    print("\t\t--running TRJCONV to remove the pbc from the trajectory..")
    
    # first TRJCONV
    cmd1 = [
        "gmx", "trjconv", "-n", "index.ndx", "-f", f"{trj_file}.xtc",
        "-o", "./nptMD_nojump_temp.xtc", "-s", tpr_file, "-pbc", "nojump", "-b", startingFrameGMXPBSA
    ]
    run_command(cmd1, input_data="0\n")
    
    # second TRJCONV
    cmd2 = [
        "gmx", "trjconv", "-n", "index.ndx", "-f", "./nptMD_nojump_temp.xtc",
        "-o", f"./{root_name}/{conf_name}_noPBC.xtc", "-s", tpr_file, "-pbc", "mol", "-center"
    ]
    run_command(cmd2, input_data="1\n0\n")

    # check
    #cmd_check = ["gmx", "check", "-f", f"./{root_name}/{conf_name}_noPBC.xtc"]
    #run_command(cmd_check,input_data=None)
    #print("\t\t--TRJCONV completed successfully!")
    # check
    output_file = "trj_check.out"
    cmd_check = ["gmx", "check", "-f", f"./cycle1_BE/cycle1_noPBC.xtc", ">", output_file, "2>&1"]
    run_command(["bash", "-c", " ".join(cmd_check)],input_data=None)
    print("\t\t--TRJCONV completed successfully!")
    

def make_index(conf_name, root_name, receptor_frag, ab_chains):
    """make index file"""
    print("\t\t--running MAKE_NDX to make index with only receptor, ligand and complex..")
        
    make_ndx_cmd = [
        "gmx", "make_ndx", "-f", f"{conf_name}_starting_protein.pdb",
        "-o", f"{root_name}/index.ndx"
    ]
    #ndx_string = f"{receptor_frag}\n{ab_chains}\nq\n"
    ndx_string = make_ndx_string_gmxpbsa(receptor_frag, ab_chains)
    
    run_command(make_ndx_cmd, input_data=ndx_string)
    print("\t\t--MAKE_NDX completed successfully!")

def create_protein_top(top_file):
    """create top with protein only"""
    print("\t\t--Creating a protein-only topology file...")
    excluded_keywords = ["SOL", "K", "CL"]  # keyword

    with open(f"{top_file}.top", "r") as infile:
        lines = infile.readlines()

    with open(f"{top_file}_protein.top", "w") as outfile:
        for line in lines:
            # check if have keywords
            if not any(keyword in line for keyword in excluded_keywords):
                outfile.write(line)

    print("\t\t--Protein topology file created successfully!")

def run_grompp(mdp_name, conf_name, top_file, root_name):
    """run GROMPP to get tpr file"""
    print("\t\t--running GROMPP to make a protein tpr..")
    
    # first GROMPP
    cmd1 = [
        "gmx", "grompp", "-v", "-f", mdp_name, "-c", f"{conf_name}_starting_protein.pdb",
        "-p", f"{top_file}_protein.top", "-o", f"{root_name}/{conf_name}.tpr", "-maxwarn", "1"
    ]
    run_command(cmd1, input_data = None)
    
    # second GROMPP
    cmd2 = [
        "gmx", "grompp", "-v", "-f", mdp_name, "-c", f"{conf_name}_starting_protein.pdb",
        "-p", f"{top_file}_protein.top", "-o", f"{root_name}/{conf_name}_newGRO.tpr", "-maxwarn", "1"
    ]
    run_command(cmd2, input_data = None)
    print("\t\t--GROMPP completed successfully!")

def count_his_residues(pdb_file):
    """calcuate the number of  HIS residuals and get the residual string"""
    print("\t\t--Counting HIS residues in the PDB file..")
    his_patterns = ["HIS CA", "HID CA", "HIE CA", "HIP CA", "HSD CA", "HSE CA", "HSP CA"]
    his_string = ""
    with open(pdb_file, "r") as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if any(pattern in line for pattern in his_patterns):
                his_string += "1\n"
    print(f"\t\t--Found HIS residues: {his_string.count('1')}")
    return his_string

def files_gmxmmpbsa(starting_gro_file, repository_pdb_file, trj_file, tpr_file, top_file, mdp_name, root_name, conf_name, vmd_function_folder, temp_files_folder, startingFrameGMXPBSA = "2000", receptor_frag = "2", ab_chains = "2"):
    
    if not check_file(f"{starting_gro_file}.gro") or not check_file(f"{trj_file}.xtc") or not check_file(f"{tpr_file}.tpr") or not check_file(f"{top_file}.top"):
        return  "Files for doing the GMXMMPBSA are not found!"
    
    # create cycle{n}_BE
    if not os.path.exists(root_name):
        os.makedirs(root_name)
    for filename in os.listdir(root_name):
        file_path = os.path.join(root_name, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

    # use make_ndx to get index.ndx
    print("\t\t--running MAKE_NDX to create index.ndx..")
    make_ndx_command = f"echo -e 'keep 1\n\nq\n' | gmx make_ndx -f {starting_gro_file}.gro -o index.ndx"
    try:
        subprocess.run(make_ndx_command, shell=True, check=True,executable = "/bin/bash")
    except subprocess.CalledProcessError:
        print("Something wrong during running MAKE_NDX!")

    # file for GRO_TO_PDB
    pathGRO = os.getcwd()
    fileNameGRO = starting_gro_file
    pathPDB = os.path.dirname(repository_pdb_file)
    pdb_name_with_extension = os.path.basename(repository_pdb_file) #xxxx.pdb
    pdb_name_without_extension = os.path.splitext(pdb_name_with_extension)[0] #xxxx
    fileNamePDB = pdb_name_without_extension
    FileNamePDB_OUT = f"{conf_name}_starting_protein"
    # RUN GRO_TO_PDB
    GRO_to_PDB(pathGRO, fileNameGRO, pathPDB, fileNamePDB, FileNamePDB_OUT, vmd_function_folder, temp_files_folder)
    
    remove_pbc(trj_file, tpr_file, startingFrameGMXPBSA, root_name, conf_name)
    make_index(conf_name, root_name, receptor_frag, ab_chains)
    create_protein_top(top_file)
    run_grompp(mdp_name, conf_name, top_file, root_name)
    his_string = count_his_residues(f"{conf_name}_starting_protein.pdb")

    print("Files created successfully!")

# 示例运行
# make_files_gmxt_pbsa("system_Compl_MD", "traj", "input.tpr", "topology.top", "energy_comp.mdp")
