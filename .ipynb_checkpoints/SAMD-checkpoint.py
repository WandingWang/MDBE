
想看这个
看呗
给你买票
感觉像画毁了
我小时候画的比这好
给你买票
看呗
# Energy checks
try:
subprocess.run(f"echo 'Temperature' | gmx energy -f NVT.edr -o temp_NVT.xvg", check=True, shell=True)
subprocess.run(f"echo 'Pressure' | gmx energy -f NPT.edr -o press_NPT.xvg", check=True, shell=True)
except subprocess.CalledProcessError:
print("Something went wrong on the energy check...")
return
# Energy checks
try:
    # Correctly select Temperature and Pressure for NVT energy file
    subprocess.run(f"echo -e 'Temperature\n0' | gmx energy -f NVT.edr -o temp_NVT.xvg", check=True, shell=True)
    
    # Correctly select Pressure for NPT energy file
    subprocess.run(f"echo -e 'Pressure\n0' | gmx energy -f NPT.edr -o press_NPT.xvg", check=True, shell=True)
    
except subprocess.CalledProcessError:
    print("Something went wrong on the energy check...")
    return
import subprocess
import os
import time

def run_gromacs_command(command, error_message, pipe_file, output_file=None):
    """
    Runs a GROMACS command and logs the output to a file.
    
    Args:
        command (str): The GROMACS command to execute.
        error_message (str): Error message to display if the command fails.
        pipe_file (str): File to write "exit" if the command fails.
        output_file (str): File to log command output.
    """
    try:
        if output_file:
            with open(output_file, "w") as out:
                subprocess.run(command, check=True, shell=True, stdout=out, stderr=subprocess.STDOUT)
        else:
            subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError:
        print(f"Something went wrong: {error_message}")
        with open(pipe_file, "w") as f:
            f.write("exit")
        raise

def make_new_minim_config_samd(input_structure_file, samd_mdp, output_gro, top_name="topol", pipe_file="out.out"):
    """
    Performs SAMD equilibration using GROMACS commands, and logs outputs to files.
    
    Args:
        input_structure_file (str): Path to the input structure file.
        samd_mdp (str): Path to the SAMD MDP file.
        output_gro (str): Path to the output GRO file.
        top_name (str): Topology file name without extension. Default is "topol".
        pipe_file (str): Pipe file for signaling errors. Default is "out.out".
    """
    # SAMD equilibration
    print(f"{time.strftime('%H:%M:%S')} -- Running SAMD MD for adaptive molecular dynamics...")
    grompp_samd_out = "gromppSAMD.out"
    mdrun_samd_out = "mdrun_SAMD.out"

    # GROMPP step for SAMD
    samd_grompp_command = (
        f"gmx grompp -f {samd_mdp} -c {input_structure_file} -r {input_structure_file} "
        f"-p {top_name}.top -o system_SAMD_MD.tpr"
    )
    run_gromacs_command(samd_grompp_command, "Something wrong on SAMD GROMPP", pipe_file, output_file=grompp_samd_out)

    # MDRUN step for SAMD
    samd_mdrun_command = (
        f"gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -bonded gpu -update gpu "
        f"-s system_SAMD_MD.tpr -c {output_gro}.gro -cpo state_SAMD_MD.cpt -x traj_SAMD_MD.xtc -e SAMD.edr -v"
    )
    run_gromacs_command(samd_mdrun_command, "Something wrong on SAMD MD RUN", pipe_file, output_file=mdrun_samd_out)

    # Energy checks for SAMD
    try:
        subprocess.run(f"printf 'Potential\n0\n' | gmx energy -f SAMD.edr -o pot_SAMD.xvg", check=True, shell=True)
        subprocess.run(f"printf 'Temperature\n0\n' | gmx energy -f SAMD.edr -o temp_SAMD.xvg", check=True, shell=True)
        subprocess.run(f"printf 'Pressure\nDensity\n0\n' | gmx energy -f SAMD.edr -o press_SAMD.xvg", check=True, shell=True)
    except subprocess.CalledProcessError:
        print("Something went wrong during the energy check...")iii        return

    print("SAMD equilibration completed successfully!")

    # Copy output files to the specified folder for checking
    os.makedirs("DOUBLE_CHECK_FOLDER", exist_ok=True)
    subprocess.run(f"cp ./grompp*.out ./*edr ./*xvg DOUBLE_CHECK_FOLDER", check=True,
