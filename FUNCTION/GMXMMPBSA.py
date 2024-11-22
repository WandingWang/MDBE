import os
import subprocess
import sys

def check_file(file_path):
    """检查指定的文件是否存在"""
    if not os.path.isfile(file_path):
        print(f"ERROR: File {file_path} not found!")
        return False
    return True

def GRO_to_PDB(pathPDB, fileNameGRO, fileNamePDB, vmd_function_folder, temp_files_folder, pathGRO = "."):
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
    # Paths and filenames
    vmd_script = os.path.join(temp_files_folder, "vmd_GRO_to_PDB.tcl")
    log_file = os.path.join(temp_files_folder, "vmd_GRO_to_PDB.out")
    function_tcl = os.path.join(vmd_function_folder, "VMD_function_GRO_to_PDB.tcl")

    # Ensure paths and required files exist
    if not os.path.exists(temp_files_folder):
        os.makedirs(temp_files_folder)
    if not os.path.isfile(function_tcl):
        raise FileNotFoundError(f"Function script missing: {function_tcl}")

    # Write VMD TCL script
    tcl_content = f"""
    variable _pathGRO "{pathGRO}"
    variable _FileNameGRO "{fileNameGRO}"
    variable _pathPDB "{pathPDB}"
    variable _FileNamePDB "{fileNamePDB}"
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

#GRO_to_PDB(pathGRO, fileNameGRO, pathPDB, fileNamePDB, vmd_function_folder, vmd_executable, temp_files_folder)

def run_command(cmd, input_data=None):
    try:
        result = subprocess.run(
            cmd, input=input_data, text=True, capture_output=True, check=True
        )
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error: {e.stderr}")
        raise

def remove_pbc(trj_file, tpr_file, startingFrameGMXPBSA, root_name, conf_name):
    """移除周期性边界条件 (PBC)"""
    print("\t\t--running TRJCONV to remove the pbc from the trajectory..")
    
    # 第一次 TRJCONV
    cmd1 = [
        "gmx", "trjconv", "-n", "index.ndx", "-f", f"{trj_file}.xtc",
        "-o", "./nptMD_nojump_temp.xtc", "-s", tpr_file, "-pbc", "nojump", "-b", str(starting_frame)
    ]
    run_command(cmd1, input="0\n")
    
    # 第二次 TRJCONV
    cmd2 = [
        "gmx", "trjconv", "-n", "index.ndx", "-f", "./nptMD_nojump_temp.xtc",
        "-o", f"./{root_name}/{conf_name}_noPBC.xtc", "-s", tpr_file, "-pbc", "mol", "-center"
    ]
    run_command(cmd2, input="1\n0\n")
    
    # 检查生成的轨迹文件
    cmd_check = ["gmx", "check", "-f", f"./{root_name}/{conf_name}_noPBC.xtc"]
    run_command(cmd_check)
    print("\t\t--TRJCONV completed successfully!")

def make_index(conf_name, root_name, receptor_frag, ab_chains):
    """生成索引文件"""
    print("\t\t--running MAKE_NDX to make index with only receptor, ligand and complex..")
    make_ndx_cmd = [
        "gmx", "make_ndx", "-f", f"{conf_name}_starting_protein.pdb",
        "-o", f"{root_name}/index.ndx"
    ]
    ndx_string = f"{receptor_frag}\n{ab_chains}\nq\n"
    run_command(make_ndx_cmd, input=ndx_string)
    print("\t\t--MAKE_NDX completed successfully!")

def create_protein_top(top_file):
    """创建仅包含蛋白质的拓扑文件"""
    print("\t\t--using HEAD to make a only-protein top..")
    with open(f"{top_file}.top", "r") as infile:
        lines = infile.readlines()
    with open(f"{top_file}_protein.top", "w") as outfile:
        outfile.writelines(lines[:-3])  # 去掉最后3行
    print("\t\t--Protein topology file created successfully!")

def run_grompp(mdp_name, conf_name, top_file, root_name):
    """运行 GROMPP 生成 tpr 文件"""
    print("\t\t--running GROMPP to make a protein tpr..")
    
    # 第一次 GROMPP
    cmd1 = [
        "gmx", "grompp", "-v", "-f", mdp_name, "-c", f"{conf_name}_starting_protein.pdb",
        "-p", f"{top_file}_protein.top", "-o", f"{root_name}/{conf_name}.tpr", "-maxwarn", "1"
    ]
    run_command(cmd1)
    
    # 第二次 GROMPP
    cmd2 = [
        "gmx", "grompp", "-v", "-f", mdp_name, "-c", f"{conf_name}_starting_protein.pdb",
        "-p", f"{top_file}_protein.top", "-o", f"{root_name}/{conf_name}_newGRO.tpr", "-maxwarn", "1"
    ]
    run_command(cmd2)
    print("\t\t--GROMPP completed successfully!")

def count_his_residues(pdb_file):
    """计算 HIS 残基数量并生成字符串"""
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
        subprocess.run(make_ndx_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Something wrong during running MAKE_NDX!")
        return 
    
    fileNameGRO = starting_gro_file
    pathPDB = os.path.dirname(repository_pdb_file)
    pdb_name_with_extension = os.path.basename(repository_pdb_file) #xxxx.pdb
    pdb_name_without_extension = os.path.splitext(pdb_name_with_extension)[0] #xxxx
    fileNamePDB = pdb_name_without_extension
    

    GRO_to_PDB(fileNameGRO, pathPDB, fileNamePDB, vmd_function_folder, temp_files_folder, pathGRO = ".")
    remove_pbc(trj_file, tpr_file, startingFrameGMXPBSA, root_name, conf_name)
    make_index(conf_name, root_name, receptor_frag, ab_chains)

    print("Files created successfully!")

# 示例运行
# make_files_gmxt_pbsa("system_Compl_MD", "traj", "input.tpr", "topology.top", "energy_comp.mdp")
