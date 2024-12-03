import csv
import os
import shutil
import logging
import numpy as np
import pandas as pd
import math

def extract_delta_energy_terms(input_filename, output_filename):
    try:
        with open(input_filename, 'r') as infile:
            lines = infile.readlines()  # 以行读取文件内容
    except FileNotFoundError:
        print(f"Error: {input_filename} file not found.")
        return

    # 查找 Delta Energy Terms 部分的开始和结束
    delta_start = False
    delta_lines = []

    for line in lines:
        # 查找 Delta Energy Terms 标题
        if "Delta Energy Terms" in line:
            delta_start = True
            continue  # 跳过该行
        # 当找到 Delta Energy Terms 部分后开始收集数据
        if delta_start:
            if line.strip() == "":  # 遇到空行结束
                break
            delta_lines.append(line.strip().split(','))  # 用逗号分割每行数据

    if not delta_lines:
        print("Error: 'Delta Energy Terms' data not found.")
        return

    # 将提取的数据写入新的 CSV 文件
    try:
        with open(output_filename, 'w', newline='') as outfile:
            writer = csv.writer(outfile)
            writer.writerows(delta_lines)  # 将提取的数据写入文件
        print(f"Delta Energy Terms extracted and saved to {output_filename}")
    except Exception as e:
        print(f"Error saving the file: {e}")

def append_results(cycle_number, results_folder):
    # 定义输出文件路径
    output_file = os.path.join(results_folder, "MoleculesResults.dat")
    # 定义输入文件路径
    input_file = os.path.join(results_folder, f"cycle{cycle_number}_results.dat")
    
    # 打开输出文件，追加写入模式
    with open(output_file, "a") as outfile:
        # 写入 Cycle_Number，左对齐，占10个字符宽度
        outfile.write(f"{cycle_number:<5}\t")
        
        # 从输入文件中提取包含 "#AVG" 的行并处理
        if os.path.exists(input_file):
            with open(input_file, "r") as infile:
                for line in infile:
                    if "#AVG" in line:
                        # 提取从第二列开始的内容（假设以制表符分隔）
                        columns = line.strip().split("\t")[1:]
                        outfile.write("\t".join(columns) + "\n")

def Data_Analysis_Pre(cycle_number_MD_FOLDER, REMOVED_FILES_FOLDER, NUMframe = "all"):
    input_filename = "gmx_MMPBSA_plot.csv"  # 输入的 CSV 文件
    output_filename = "delta_energy_terms.csv"  # 输出的 CSV 文件
    extract_delta_energy_terms(input_filename, output_filename)
    # read the csv file
    df = pd.read_csv('delta_energy_terms.csv')

    # create df
    output_data = []
    for index, row in df.iterrows():
        # for each line
    
        # Frame #
        frame = row['Frame #']  
    
        # DeltaG(kJ/mol)
        delta_g = row['TOTAL']  
    
        # Coul(kJ/mol) = EEL + EEL14
        EEL = row['EEL']  
        EEL14 = row['1-4 EEL']  
        coul = EEL + EEL14
    
        # vdW(kJ/mol) = VDW + VDW14
        VDW = row['VDWAALS']  
        VDW14 = row['1-4 VDW']  
        vdW = VDW + VDW14
    
        # PolSol(kJ/mol) = EPB + ENPOLAR (if have)
        # for this infile the NpoSol is the value of [ENPOLAR], so it is included in the PolSol

        if pd.notna(row['EPB']) and pd.notna(row['ENPOLAR']):
            EPB = row['EPB']
            ENPOLAR = row['ENPOLAR']
            pol_sol = EPB + ENPOLAR
        else:
            pol_sol = row['EGB']  # if have
    
        # NpoSol(kJ/mol) = EDISPER or ESURF (if have)
        non_pol_solv = 0
        if pd.notna(row['EDISPER']):
            non_pol_solv = row['EDISPER']
        elif pd.notna(row['ESURF']):
            non_pol_solv = row['ESURF']
    
        # add the values 
        output_data.append({
            '# frame': frame,
            'DeltaG(kJ/mol)': delta_g,
            'Coul(kJ/mol)': coul,
            'vdW(kJ/mol)': vdW,
            'PolSol(kJ/mol)': pol_sol,
            'NpoSol(kJ/mol)': non_pol_solv
        })

    # transfer
    output_df = pd.DataFrame(output_data)

    # write it in
    output_df.to_csv('energy_plot_temp.csv', sep='\t', index=False)

    # 配置参数
    #NUMframe = "all"  # 或者一个特定的数字

    #logging.info("Cleaning files.")
    # 处理 energy_plot_temp 文件
    with open('./energy_plot_temp.csv', 'r') as f:
        lines = f.readlines()

    # 获取表头
    #header = lines[0]  # 第一行是表头

    # 根据 NUMframe 判断截取行
    if NUMframe == "all":
        selected_lines = lines[1:]  # 如果是 "all"，保留所有内容，包括表头
    else:
        selected_lines = lines[int(NUMframe):]  # 保留表头，并从指定行开始截取

    # 保存处理后的文件
    with open('../energy_plot_temp.csv', 'w') as f:
        f.writelines(selected_lines)


    # 清理文件
    # 移动匹配的文件到指定目录
    for filename in os.listdir('.'):
        if filename.startswith('_GMXMMPBSA'):
            shutil.move(filename, os.path.join(REMOVED_FILES_FOLDER, filename))
        elif '_temp' in filename:
            shutil.move(filename, os.path.join(REMOVED_FILES_FOLDER, filename))

    # 删除备份和 .ff 文件
    for filename in os.listdir('.'):
        if filename.endswith(('#', '~', '.ff')):
            os.remove(filename)

    # 记录清理操作
    with open(os.path.join(REMOVED_FILES_FOLDER, 'rm.out'), 'w') as rm_out_file:
        rm_out_file.write("Removed backup and .ff files.")

    # 切换到指定目录
    try:
        os.chdir(cycle_number_MD_FOLDER)
    except FileNotFoundError:
        raise SystemExit(f"Cannot enter '{cycle_number_MD_FOLDER}' folder")

def Data_Analysis_Cal(cycle_number, results_folder):
    
    logging.info("Data Analysis.")
    # 读取数据
    input_file = 'energy_plot_temp.csv'  
    # NEED TO CHANGE cycle${Cycle_Number}_results.dat
    output_file = f'cycle{cycle_number}_results.dat'

    # 假设文件数据为 tab-separated
    df = pd.read_csv(input_file, delim_whitespace=True, header=None)

    # 为了清晰，假设文件列对应以下字段
    df.columns = ['frame', 'DeltaG', 'Coul', 'VdW', 'PolSol', 'NpoSol']

    # 初始化变量
    Population = len(df)
    DeltaG = df['DeltaG'].sum()
    Coul = df['Coul'].sum()
    VdW = df['VdW'].sum()
    PolSol = df['PolSol'].sum()
    NpoSol = df['NpoSol'].sum()

    # 计算 SF1 和 SF2
    df['SF1'] = (df['Coul'] / 10) - (df['PolSol'] / 10) + (df['NpoSol'] * 10)
    df['SF2'] = (3 * df['Coul']) + df['PolSol']

    # 初始化变量
    Canonical_AVG = 0.0
    Canonical_AVG_w = 0.0

    for _, row in df.iterrows():
        # 尝试将 DeltaG 转换为浮动类型
        DeltaG_temp = row[1]  # 使用列名访问 DeltaG
        # 计算指数加权值
        weight = math.exp(-int(DeltaG_temp / 2.479))
        # 累加计算
        Canonical_AVG += DeltaG_temp * weight
        Canonical_AVG_w += weight

    # 归一化
    if Canonical_AVG_w != 0:
        Canonical_AVG /= Canonical_AVG_w

    # 平均值计算
    mean_DeltaG = DeltaG / Population
    mean_Coul = Coul / Population
    mean_VdW = VdW / Population
    mean_PolSol = PolSol / Population
    mean_NpoSol = NpoSol / Population

    mean_SF1 = mean_Coul / 10 - mean_PolSol / 10 + mean_NpoSol * 10
    mean_SF2 = (3 * mean_Coul) + mean_PolSol


    # 计算标准差
    std_DeltaG = np.std(df['DeltaG'], ddof=1)
    std_Coul = np.std(df['Coul'], ddof=1)
    std_VdW = np.std(df['VdW'], ddof=1)
    std_PolSol = np.std(df['PolSol'], ddof=1)
    std_NpoSol = np.std(df['NpoSol'], ddof=1)
    std_SF1 = std_Coul / 10 - std_PolSol / 10 + std_NpoSol * 10
    std_SF2 = (3 * std_Coul) + std_PolSol

    # 计算 DeltaG 在2sigma内的值
    DeltaG_2s = df[(np.abs(df['DeltaG'] - mean_DeltaG) < 2 * std_DeltaG)]['DeltaG'].mean()

    # 计算中位数
    median_DeltaG = df['DeltaG'].median()
    # 计算标准差
    std_DeltaG = np.std(df['DeltaG'], ddof=1)  # ddof=1 表示使用样本标准差公式（贝塞尔校正）


    # 创建警告信息列表
    warnings = []

    # 遍历每一行，检查是否超出 2 sigma
    for i, row in df.iterrows():
        var = np.abs(row['DeltaG'] - mean_DeltaG)
        if var >= 2 * std_DeltaG:
            warnings.append(f"# WARNING: frame {row['frame']} is out of 2 sigma!!")

    # 输出处理结果
    with open(output_file, 'w') as f:
        # 写入头部信息
        f.write("# SF1=Coulomb/10-PolarSolvation/10+Non-PolarSolvation*10\n")
        f.write("# SF2=3*Coulomb+PolarSolvation\n")
        f.write("# C_AVG=norm(SUM Gi*e^BGi)\n")
        f.write(f"#frame\tDeltaG(kJ/mol)\tCoul(kJ/mol)\tVdW(kJ/mol)\tPolSol(kJ/mol)\tNpoSol(kJ/mol)\tSF1\tSF2\n")

        # 写入数据
        for _, row in df.iterrows():
            f.write(f"{row['frame']:<10}{row['DeltaG']:>12.3f}{row['Coul']:>13.3f}{row['VdW']:>13.3f}{row['PolSol']:>13.3f}{row['NpoSol']:>13.3f}{row['SF1']:>13.3f}{row['SF2']:>13.3f}\n")
    
         # 写入警告信息
        for warning in warnings:
            f.write(f"{warning}\n")   
   
        # 写入汇总信息表头
        f.write("\n# FINAL RESULTS\n")
        f.write(f"#frame\t{'DeltaG(kJ/mol)':>15}\t{'Coul(kJ/mol)':>15}\t{'VdW(kJ/mol)':>15}\t{'PolSol(kJ/mol)':>15}\t{'NpoSol(kJ/mol)':>15}\t{'SF1':>15}\t{'SF2':>15}\t{'Canonical_AVG':>15}\t{'MedianDeltaG(kJ/mol)':>15}\t{'DeltaG_2s(kJ/mol)':>15}\n")
        # 写入汇总结果
        f.write(f"#AVG\t{mean_DeltaG:>15.1f}\t{mean_Coul:>15.1f}\t{mean_VdW:>15.1f}\t{mean_PolSol:>15.1f}\t{mean_NpoSol:>15.1f}\t{mean_SF1:>15.1f}\t{mean_SF2:>15.1f}\t{Canonical_AVG:>15.1f}\t{median_DeltaG:>15.1f}\t{DeltaG_2s:>15.1f}\n")
        f.write(f"#STD\t{std_DeltaG:>15.1f}\t{std_Coul:>15.1f}\t{std_VdW:>15.1f}\t{std_PolSol:>15.1f}\t{std_NpoSol:>15.1f}\t{std_SF1:>15.1f}\t{std_SF2:>15.1f}\t{"nan":>15}\t{"nan":>15}\t{std_DeltaG:>15.1f}\n")


    try:
        with open(f'cycle{cycle_number}_results.dat', 'r') as file:
            lines = file.readlines()
    except FileNotFoundError:
        logging.error(f"File not found: cycle{cycle_number}_results.dat")
        exit()

    # 初始化变量存储表头、平均值和标准差
    header = None
    avg_values = None
    std_values = None

    # 逐行解析文件内容
    for line in lines:
        if line.startswith("#frame"):
            header = line.strip().lstrip("#").split()
        elif line.startswith("#AVG"):
            avg_values = line.strip().lstrip("#").split()
        elif line.startswith("#STD"):
            std_values = line.strip().lstrip("#").split()

    # 验证数据完整性
    if not header or not avg_values or not std_values:
        logging.error("Missing required data (header, AVG, or STD).")
        exit()

    # 格式化日志输出
    log_message = (
        f"Results for cycle{cycle_number}:\n"
        f"Headers: {' | '.join(header)}\n"
        f"AVG Values: {' | '.join(avg_values)}\n"
        f"STD Values: {' | '.join(std_values)}")

    # 输出合并后的日志信息
    logging.info(log_message)
    # move outpufile (cycle1_results.dat) to results folder
    shutil.move(output_file, os.path.join(results_folder, output_file))
    # append cycle data on MoleculeResults.dat
    append_results(cycle_number, results_folder)