import os
import subprocess
import sys

def check_file(file_path):
    """检查指定的文件是否存在"""
    if not os.path.isfile(file_path):
        print(f"ERROR: File {file_path} not found!")
        return False
    return True

def GRO_to_PDB(pathGRO=".", FileNameGRO="", pathPDB=".", FileNamePDB="", FileNamePDB_OUT=""):
    # Ensure ABiSS_LIB, VMD, and fatal are defined as global variables if needed
    # GLOBAL VARIABLES NEEDED:
    # ABiSS_LIB, VMD, fatal (These should be defined somewhere in your environment or Python code)

    # LOCAL VARIABLES:
    if not FileNamePDB_OUT:
        FileNamePDB_OUT = FileNameGRO

    # Writing TCL script for VMD execution
    vmd_tcl_script = f"""
    variable _pathGRO "{pathGRO}"
    variable _FileNameGRO "{FileNameGRO}"
    variable _pathPDB "{pathPDB}"
    variable _FileNamePDB "{FileNamePDB}"
    variable _FileNamePDB_OUT "{FileNamePDB_OUT}"
    """
    
    # Path to TEMP_FILES_FOLDER and ABiSS_LIB should be set somewhere in your environment or code
    TEMP_FILES_FOLDER = "/path/to/temp_files"  # Set to your actual temporary folder
    ABiSS_LIB = "/path/to/ABiSS_LIB"  # Set to your actual ABiSS_LIB folder
    VMD = "/path/to/vmd"  # Set to your VMD executable path

    vmd_script_file = os.path.join(TEMP_FILES_FOLDER, "vmd_GRO_to_PDB.tcl")
    with open(vmd_script_file, 'w') as f:
        f.write(vmd_tcl_script)
        with open(os.path.join(ABiSS_LIB, "VMD_function_GRO_to_PDB.tcl"), 'r') as lib_file:
            f.write(lib_file.read())

    # Running VMD command
    try:
        result = subprocess.run([VMD, "-dispdev", "none", "-e", vmd_script_file], 
                                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(os.path.join(TEMP_FILES_FOLDER, "vmd_GRO_to_PDB.out"), 'wb') as log_file:
            log_file.write(result.stdout)
            log_file.write(result.stderr)

        if result.returncode != 0:
            raise Exception(f"GRO_to_PDB -> {VMD} failed.. check {TEMP_FILES_FOLDER}/vmd_GRO_to_PDB.out")
    
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(33)

# Example usage (you can set the parameters as needed):
# GRO_to_PDB(pathGRO=".", FileNameGRO="example", pathPDB=".", FileNamePDB="output", FileNamePDB_OUT="final_output")

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
