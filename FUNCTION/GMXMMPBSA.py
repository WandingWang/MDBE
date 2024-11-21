import os
import subprocess
import sys

def check_file(file_path):
    """检查指定的文件是否存在"""
    if not os.path.isfile(file_path):
        print(f"ERROR: File {file_path} not found!")
        return False
    return True

def GRO_to_PDB(pathPDB, fileNameGRO, fileNamePDB, vmd_function_folder, vmd_executable, temp_files_folder, pathGRO = "."):
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
        [vmd_executable, "-dispdev", "none", "-e", vmd_script],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )

    # Log results
    with open(log_file, 'w') as log:
        log.write(result.stdout + result.stderr)

    if result.returncode != 0:
        raise RuntimeError(f"VMD execution failed. Check log: {log_file}")

    print(f"GRO to PDB conversion completed successfully: {fileNamePDB}")

#GRO_to_PDB(pathGRO, fileNameGRO, pathPDB, fileNamePDB, vmd_function_folder, vmd_executable, temp_files_folder)

def files_gmxmmpbsa(starting_gro_file, trj_file, tpr_file, top_file, mdp_name, root_name="config1"):
    """根据输入文件生成所需的 GROMACS 相关文件"""
    
    # 检查输入文件是否存在
    if not check_file(f"{starting_gro_file}.gro") or not check_file(f"{trj_file}.xtc") or not check_file(f"{tpr_file}.tpr") or not check_file(f"{top_file}.top"):
        return   # 如果任何一个输入文件不存在，返回错误
    
    # 创建工作目录并清空已有文件
    if not os.path.exists(root_name):
        os.makedirs(root_name)
    for filename in os.listdir(root_name):
        file_path = os.path.join(root_name, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

    # 使用 make_ndx 生成索引文件
    print("\t\t--running MAKE_NDX to create index.ndx..")
    make_ndx_command = f"echo -e 'keep 1\n\nq\n' | gmx make_ndx -f {starting_gro_file}.gro -o index.ndx"
    try:
        subprocess.run(make_ndx_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error running MAKE_NDX!")
        return 
    
    # 使用 TRJCONV 处理轨迹文件，生成去除 PBC 的轨迹文件
    print("\t\t--running TRJCONV to remove PBC from the trajectory..")
    trjconv_command = f"echo '0' | trjconv -n index.ndx -f {trj_file}.xtc -o {root_name}/{root_name}_noPBC.xtc -s {tpr_file}.tpr -pbc nojump -b 0"
    try:
        subprocess.run(trjconv_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error running TRJCONV!")
        return 
    
    # 使用 GROMPP 生成仅含蛋白质的 tpr 文件
    print("\t\t--using GROMPP to create a protein tpr file..")
    grompp_command = f"grompp -v -f {mdp_name} -c {root_name}_starting_protein.pdb -p {top_file}_protein.top -o {root_name}/{root_name}.tpr -maxwarn 1"
    try:
        subprocess.run(grompp_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error running GROMPP!")
        return 

    print("Files created successfully!")

# 示例运行
# make_files_gmxt_pbsa("system_Compl_MD", "traj", "input.tpr", "topology.top", "energy_comp.mdp")
