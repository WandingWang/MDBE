import subprocess
import os
import time

def run_gromacs_command(command, error_message, pipe_file):
    try:
        subprocess.run(command, check=True, shell=True)
    except subprocess.CalledProcessError:
        print(f"Something went wrong: {error_message}")
        with open(pipe_file, "w") as f:
            f.write("exit")
        exit()

def make_new_minim_nvt_npt(input_structure_file, nvt_mdp, npt_mdp, top_name="topol", output_gro, pipe_file="out.out"):
    # Local variables
    #top_name = top_name.rstrip(".top")
    #output_gro = output_gro.rstrip(".gro")

    # Simulated Annealing-MD (this part is commented out in the original script)
    # Run simulated annealing - simulated annealing code can be added here
    # Example:
    # samd_command = f"$GROMPP -f {samd_mdp} -c {input_structure_file} -r {input_structure_file} -p {top_name}.top -o system_SAMD.tpr -maxwarn 1"
    # run_gromacs_command(samd_command, "Something wrong on SA-MD GROMPP", pipe_file)
    # run_gromacs_command(f"$MDRUN_md -s system_SAMD.tpr -c system_SAMD.gro -cpo state_SAMD.cpt -x traj_SAMD.xtc", "Something wrong on SA-MD RUN", pipe_file)

    # NVT
    print(f"{time.strftime('%H:%M:%S')} -- Running NVT MD for temperature equilibration...")
    nvt_command = f"gmx grompp -f {nvt_mdp} -c {input_structure_file} -r {input_structure_file} -p {top_name}.top -o system_NVT_MD.tpr"
    run_gromacs_command(nvt_command, "Something wrong on NVT GROMPP", pipe_file)

    run_gromacs_command(f"gmx mdrun -ntmpi 1 -ntomp 8 -nb -gpu -pme -bonded gpu -update gpu -s system_NVT_MD.tpr -c system_NVT_MD.gro -cpo state_NVT_MD.cpt -e NVT.edr -v", "Something wrong on NVT MD RUN", pipe_file)

    # NPT 
    print(f"{time.strftime('%H:%M:%S')} -- Running NPT MD for pressure equilibration...")
    npt_command = f"gmx grompp -f {npt_mdp} -c system_NVT_MD.gro -r system_NVT_MD.gro -p {top_name}.top -o system_NPT_MD.tpr -t state_NVT_MD.cpt -maxwarn 1"
    run_gromacs_command(npt_command, "Something wrong on NPT GROMPP", pipe_file)

    run_gromacs_command(f"gmx mdrun -ntmpi 1 -ntomp 8 -nb -gpu -pme -bonded gpu -update gpu -s system_NPT_MD.tpr -c {output_gro}.gro -cpo state_NPT_MD.cpt -x traj_NPT_MD.xtc -e NPT.edr -v", "Something wrong on NPT MD RUN", pipe_file)

    # Energy checks
    try:
        subprocess.run(f"printf 'Temperature\n0\n' | gmx energy -f NVT.edr -o temp_NVT.xvg", check=True, shell=True)
        subprocess.run(f"printf 'Pressure\nDensity\n0\n' | gmx energy -f NPT.edr -o press_NPT.xvg", check=True, shell=True)
    except subprocess.CalledProcessError:
        print("Something went wrong on the energy check...")
        exit()

    # Copy output files to the specified folder for checking
    os.makedirs("DOUBLE_CHECK_FOLDER", exist_ok=True)
    subprocess.run(f"cp ./grompp*_seq*out ./*edr ./*xvg DOUBLE_CHECK_FOLDER", check=True, shell=True)
