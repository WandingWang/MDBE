import subprocess
import os

out_file = "fill_water_ions.out"  

def fill_water_ions(starting_system, topology, mdp_file):
    """
    Function to add water and ions to a system.
    - starting_system: GRO file from last step
    - topology: Topology file for the system.
    - mdp_file: MDP file.
    """

    # Local variables
    mdp_file = mdp_file 

    # Check file readability
    if not os.path.exists(starting_system + ".gro"):
        print(f"ERROR: I cannot read SYSTEM_file({starting_system}.gro)")
        return 1
    if not os.path.exists(topology + ".top"):
        print(f"ERROR: I cannot read TOPOLOGY_file({topology}.top)")
        return 1
    if not os.path.exists(mdp_file):
        print(f"ERROR: I cannot read MDP_file({mdp_file})")
        return 1
    

    # Add water to the system using gmx solvate
    genbox_cmd = f"gmx solvate  -cp {starting_system}.gro -cs spc216.gro -o system_water.gro -p {topology}.top"
    try:
        with open(out_file, 'a') as log:
            subprocess.run(genbox_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong on gmx solvate!")
        return 2

    # preparing input file for gmx genion by gmx grommp
    grompp_cmd = f"gmx grompp -f {mdp_file} -c system_water.gro -p {topology}.top -o system_ions.tpr -maxwarn 1"
    try:
        with open(out_file, 'a') as log:
            subprocess.run(grompp_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong on gmx grompp!")
        return 3

    #add ions
    genion_cmd = f'echo "SOL" | gmx genion -s system_ions.tpr -o system_ions.gro -p {topology}.top -pname K -nname CL -neutral'

    try:
        with open(out_file, 'a') as log:
            subprocess.run(genion_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong on genion!")
        return 4
    print("Water and ions successfully added to the system!")

# Example usage:
# fill_water_ions("system.gro", "topol.top", "ions.mdp")