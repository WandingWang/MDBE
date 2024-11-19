import subprocess
import os

mdp_NAME = "ions.mdp"  
out_file = "fill_water_ions.out"  

def fill_water_ions(starting_system, topology, mdp_file):
    """
    Function to add water and ions to a system.
    Args:
    - starting_system: GRO file from last step
    - topology: Topology file for the system.
    - mdp_file: MDP file.
    """

    # Local variables
    mdp_file = mdp_file or mdp_NAME

    # Check file readability
    if not os.path.exists(starting_system):
        print(f"ERROR: I cannot read SYSTEM_file({starting_system})")
        return 1
    if not os.path.exists(topology):
        print(f"ERROR: I cannot read TOPOLOGY_file({topology})")
        return 1
    if not os.path.exists(mdp_file):
        print(f"ERROR: I cannot read MDP_file({mdp_file})")
        return 1
    
    topology = topology.rstrip('.top')

    # Add water to the system using gmx solvate
    genbox_cmd = f"gmx solvate  -cp {starting_system} -cs spc216.gro -o system_water.gro -p {topology}.top"
    try:
        with open(out_file, 'a') as log:
            subprocess.run(genbox_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong on gmx solvate!! Exiting...")
        return 2

    # preparing input file for gmx genion by gmx grommp
    grompp_cmd = f"gmx grompp -f {mdp_file} -c system_water.gro -p {topology}.top -o system_ions.tpr -maxwarn 1"
    try:
        with open(out_file, 'a') as log:
            subprocess.run(grompp_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong on gmx grompp!! Exiting...")
        return 3
    """

    # Step 3: Extract charge value 
    with open(out_file, 'r') as f:
        for line in f:
            if "System has non-zero total charge" in line:
                ChargeValue = float(line.split()[-1])
                ChargeValue = ChargeValue + 0.5 if ChargeValue > 0 else ChargeValue - 0.5
                break
        else:
            ChargeValue = 0  # Default to 0 if no charge value is found

    # Step 4: Calculate ions number
    with open(f"{topology}.top", 'r') as f:
        WaterMolecules = sum(1 for line in f if "SOL" in line)

    IonsNumber = int((WaterMolecules * 0.00271) // 1)
    KNumber = int(((IonsNumber - 1 * ChargeValue / 2) // 1))
    ClNumber = int(((IonsNumber + 1 * ChargeValue / 2 + 1 * ChargeValue % 2) // 1))

    # Step 5: Check net charge
    netCHARGE = int((KNumber - ClNumber + 1 * ChargeValue) // 1)
    NEWnetCHARGE = "not computed"
    print(f"ChargeValue={ChargeValue} WaterMolecules={WaterMolecules} IonsNumber={IonsNumber} "
          f"netCHARGE={netCHARGE} KNumber={KNumber} ClNumber={ClNumber} NEWnetCHARGE={NEWnetCHARGE}", file=open("SystemCharge.out", 'w'))

    if netCHARGE != 0:
        print("\n>>> WARNING!! Check the charge of the system!! (SystemCharge.out) <<<", file=open("SystemCharge.out", 'a'))
        ClNumber = int((ClNumber + netCHARGE) // 1)
        NEWnetCHARGE = int((KNumber - ClNumber + 1 * ChargeValue) // 1)

    print(f"After netCHARGE CHECK: netCHARGE={netCHARGE} KNumber={KNumber} ClNumber={ClNumber} NEWnetCHARGE={NEWnetCHARGE}",
          file=open("SystemCharge.out", 'a'))

    # Step 6: Generate index file with MAKE_NDX
    make_ndx_cmd = f"echo -e 'keep 0\nr SOL\nkeep 1\n\nq\n' | gmx make_ndx -f system_water.gro -o index_SOL.ndx"
    try:
        with open(out_file, 'a') as log:
            subprocess.run(make_ndx_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong with MAKE_NDX!! Exiting...")
        return 4

    # Step 7: Add ions using GENION
    genion_cmd = f"echo -e '0\n' | gmx genion -s system_ions.tpr -n index_SOL.ndx -o system_ions.gro -p {topology}.top -nn {ClNumber} -nname CL -np {KNumber} -pname K"
    try:
        with open(out_file, 'a') as log:
            subprocess.run(genion_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong on GENION!! Exiting...")
        return 5

    # Clean up temporary index file
    os.remove("index_SOL.ndx")
    """
    #add ions
    genion_cmd = f'echo "SOL" | gmx genion -s system_ions.tpr -o system_ions.gro -p {topology}.top -pname K -nname CL -neutral'

    try:
        with open(out_file, 'a') as log:
            subprocess.run(genion_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong on genion!! Exiting...")
        return 4
    print("Water and ions successfully added to the system!")

# Example usage:
# fill_water_ions("system.gro", "topol.top", "ions.mdp")