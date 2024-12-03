import os
import shutil
import logging



# 清理文件 NO SUCH FILE NEED TO IMPROVE
def clean_up(removed_files_folder):
    logging.info("Cleaning files.")
    try:
        for pattern in ["*#*", "*~*", "*temp*"]:
            for filename in os.listdir('.'):
                if filename.startswith(pattern):
                    file_path = os.path.join('.', filename)
                    if os.path.isdir(file_path):
                        shutil.rmtree(file_path)
                    else:
                        os.remove(file_path)
        # 记录清理操作的日志
        with open(os.path.join(removed_files_folder, 'rm.out'), 'w') as f:
            f.write("Cleanup completed.")
    except Exception as e:
        logging.error(f"Error during cleanup: {e}")

# 移动和复制文件
def move_and_copy_files(cycle_number, repository_folder, cycle_number_MD_folder, root_name):
    logging.info("Moving files.")
    try:
        # 移动文件
        for filename in os.listdir('.'):
            if filename.startswith(f'RUN1_cycle{cycle_number}_prot') or filename.endswith('.out') or filename == root_name:
                shutil.move(filename, os.path.join(repository_folder, filename))
        
        # 复制文件
        shutil.copy('system_Compl_MDstart.gro', os.path.join(repository_folder, f'system_cycle{cycle_number}_MDstart.gro'))
        if os.path.isdir(cycle_number_MD_folder):
            shutil.copytree(cycle_number_MD_folder, os.path.join(repository_folder, f'cycle{cycle_number}_MD'))
        
    except Exception as e:
        logging.error(f"Error during file operations: {e}")

def clean_for_each_cycle(cycle_number, repository_folder, cycle_number_MD_folder, root_name, removed_files_folder, current_conf_path):
    clean_up(removed_files_folder)
    move_and_copy_files(cycle_number, repository_folder, cycle_number_MD_folder, root_name)

    # 假设这些变量已定义
    #removed_files_folder = folders["REMOVED_FILES_FOLDER"] # 被删除文件的目标文件夹
    #current_conf_path = configuration_path  # 当前配置路径
    #cycle_number = 1  # 示例初始周期号

    # 移动文件到 removed_files_folder
    for filename in os.listdir('.'):
        if os.path.isfile(filename):  # 检查是否是普通文件
            shutil.move(filename, os.path.join(removed_files_folder, filename))

    # 删除当前目录下的所有文件
    for filename in os.listdir('.'):
        file_path = os.path.join('.', filename)
        if os.path.isfile(file_path):
            os.remove(file_path)

    # 删除特定的目录
    cycle_number_md_folder = os.path.join(current_conf_path, f"cycle{cycle_number}_MD")
    if os.path.exists(cycle_number_md_folder):
        shutil.rmtree(cycle_number_md_folder)

