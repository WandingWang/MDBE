import os
import subprocess

# 检查工具是否可用的函数
def check_tool(tool_name):
    """检查是否能找到指定的工具（例如grompp, trjconv等）"""
    result = subprocess.run([tool_name, '-h'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"ERROR: Cannot find {tool_name}!")
        return False
    return True

# 检查文件是否存在的函数
def check_file(file_path):
    """检查指定的文件是否存在"""
    if not os.path.isfile(file_path):
        print(f"ERROR: File {file_path} not found!")
        return False
    return True

# 主函数
def make_files_gmxt_pbsa(starting_gro_file, trj_file, tpr_file, top_file, mdp_name, root_name="config1"):
    """根据输入文件生成所需的 GROMACS 相关文件"""
    
    # 检查必要的工具是否存在
    required_tools = ["grompp", "trjconv", "make_ndx", "check"]
    for tool in required_tools:
        if not check_tool(tool):
            return 1  # 如果任何一个工具找不到，返回错误
    
    # 检查输入文件是否存在
    if not check_file(f"{starting_gro_file}.gro") or not check_file(f"{trj_file}.xtc") or not check_file(f"{tpr_file}.tpr") or not check_file(f"{top_file}.top"):
        return 2  # 如果任何一个输入文件不存在，返回错误
    
    # 创建工作目录并清空已有文件
    if not os.path.exists(root_name):
        os.makedirs(root_name)
    for filename in os.listdir(root_name):
        file_path = os.path.join(root_name, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

    # 使用 make_ndx 生成索引文件
    print("\t\t--running MAKE_NDX to create index.ndx..")
    make_ndx_command = f"echo -e 'keep 1\n\nq\n' | make_ndx -f {starting_gro_file}.gro -o index.ndx"
    try:
        subprocess.run(make_ndx_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error running MAKE_NDX!")
        return 3
    
    # 使用 TRJCONV 处理轨迹文件，生成去除 PBC 的轨迹文件
    print("\t\t--running TRJCONV to remove PBC from the trajectory..")
    trjconv_command = f"echo '0' | trjconv -n index.ndx -f {trj_file}.xtc -o {root_name}/{root_name}_noPBC.xtc -s {tpr_file}.tpr -pbc nojump -b 0"
    try:
        subprocess.run(trjconv_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error running TRJCONV!")
        return 4
    
    # 使用 GROMPP 生成仅含蛋白质的 tpr 文件
    print("\t\t--using GROMPP to create a protein tpr file..")
    grompp_command = f"grompp -v -f {mdp_name} -c {root_name}_starting_protein.pdb -p {top_file}_protein.top -o {root_name}/{root_name}.tpr -maxwarn 1"
    try:
        subprocess.run(grompp_command, shell=True, check=True)
    except subprocess.CalledProcessError:
        print("Error running GROMPP!")
        return 5

    print("Files created successfully!")

# 示例运行
make_files_gmxt_pbsa("start", "traj", "input.tpr", "topology.top", "energy_comp.mdp")