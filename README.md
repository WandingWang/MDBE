# MDFreeEnergy  
## **Overview**
This pipeline is designed for **Binding Free Energy Calculation Using Multiple Short MD Simulations**. 
### **The calculation method**  
The method involves running 10 short MD simulations. The binding free energy is then calculated for each simulation, and the results are averaged to obtain a final estimate of the binding energy.  

It automates the process from input preparation to result generation.

---

## **Installation**
### **Dependencies**
- VMD: use for generating files for gmx_MMPBSA
- GROMACS (GPU version)
- gmx_MMPBSA: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/installation/  
- Other dependencies: Pandas, Numpy, Yaml, if you don't have them, please use 
   ```bash
   conda install xxxx (like pandas)  

### **Installation Steps**
1. Clone the repository:
   ```bash
   git clone https://github.com/WandingWang/MDFreeEnergy.git
   cd MDFreeEnergy
2. Make sure you have all dependencies: check conda, python, vmd, gromacs, gmx_mmpbsa;
3. Open the input file **[input.yaml]**, change parameters: **input setting, IMPORTANT Basic setting, gmx_mmpbsa, run**, You need to set the input file path, conda activation path, etc., as well as the number of chains, etc., based on your own requirements.
4. When you make sure that all dependencies and setting are ok, run main script:
   ```bash
   python MDFreeEnergy.py  
5. When it finished, you can get the energy information in the **[OUTPUT.out]**

---

### **Steps to Prepare Your Structure for the pipeline**  
1. **Obtain a Stable Conformation**: You will have a stable conformation of your system (e.g., the last frame of your NPT simulation). This is typically stored in a .gro file.  
2. **Convert the GRO File to PDB**: You can convert the .gro file (which does not contain chain information) to a .pdb file using 'gmx editconf' or other tools like VMD. This command 'gmx editconf' will generate a PDB file from the GRO file, but without chain information.  
3. **Remove Water and Ions**: Please remove water and ions. We will apply a uniform standard to avoid errors.   
4. **Ensure Chain ID Information**: Ensure that chain IDs are properly assigned in the PDB file. The PDB format uses the CHAIN field to denote which molecules belong to which chain. If the CHAIN information is missing or improperly assigned, you can manually edit the PDB file to add chain IDs or use VMD or PyMOL to check and add chain IDs.  
5. **Make Sure Receptor firstly and then Antibody in your PDB file**  
