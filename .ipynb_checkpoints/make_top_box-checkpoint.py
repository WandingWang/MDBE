import os
import subprocess
import glob

def make_top_protein(protein_infile, forcefield, watermodel, protein_outfile, topfile, merge=""):
    """
    generate GROMACS topology (.top and .itp) from .pdb
    buid box to get .gro

    parameters:
        protein_infile (str): input file PDB name(no extension);
        forcefield (str): forcrfild used;
        watermodel (str): water model used;
        protein_outfile (str): output filr GRO name (no extension);
        topfile (str): TOP name (no extension);
        merge (str): wheather to interactively merge histidine residues" 'interactive' or ''
    """
   
    # output logging 
    out_file = "MakeTOP_protein.out"

    merge_string = "y\nn\n" if merge == "interactive" else ""
    pdb2gmx_string = "\n1\n"
    num_his = 0

    # open .pdb
    with open(f"./{protein_infile}.pdb", 'r') as infile:
        pdb_data = infile.readlines()

    # get the number of histidine
    his_residues = ["HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP"]
    for line in pdb_data:
        if any(his in line for his in his_residues):
            num_his += 1
    
    # set pdb2gmx_string
    if merge == "interactive":
        pdb2gmx_string += merge_string
    for _ in range(num_his):
        pdb2gmx_string += "1\n"

    with open(out_file, 'a') as log:
        log.write(f"_protein_OutFileNAME: {protein_outfile}\n")
        log.write(f"_ForceField: {forcefield}\n")
        log.write(f"_protein_InFileNAME: {protein_infile}\n")
        log.write(f"_topFileNAME: {topfile}\n")
        log.write(f"_merge: {merge}\n")
        log.write(f"mergeSTRING: {merge_string}\n")
        log.write(f"num hist: {num_his}\npdb2gmx_string: {pdb2gmx_string}\n\n")

    # use pdb2gmx to generate topology
    pdb2gmx_cmd = f"gmx pdb2gmx -f {protein_infile}.pdb -o system.pdb -p {topfile}.top -ignh -ff {forcefield} -water {watermodel}"
    try:
        subprocess.run(pdb2gmx_cmd, shell=True, check=True, stdout=open(out_file, 'a'), stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong with pdb2gmx!")
        return


    print("Topology generation completed successfully.")

    # input and ouput for buiding box
    input_file = "system.pdb"
    output_file = f"{protein_outfile}.gro"

    # box type
    editconf_option = "-bt triclinic -d 1.5"  

    try:
        # run editconf to buid box
        editconf_cmd = f"gmx editconf -f {input_file} {editconf_option} -o {output_file}"
        # log
        with open(out_file, 'a') as log:
            subprocess.run(editconf_cmd, shell=True, check=True, stdout=log, stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError:
        print("Something went wrong with editconf!")

    print("Simulation box definition completed successfully.")