# MDFreeEnergy  
## **Overview**
This pipeline is designed for **Binding Free Energy Calculation Using Multiple Short MD Simulations**. 
### **The calculation method**  
The method involves running 10 short MD simulations. The binding free energy is then calculated for each simulation, and the results are averaged to obtain a final estimate of the binding energy.  

It automates the process from input preparation to result generation.

---

## **Installation**
### **Dependencies**
- **VMD**: use for generating files for gmx_MMPBSA
- **GROMACS** (GPU version)
- **gmx_MMPBSA**: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/installation/ (**please do not use gmx_mpi when you run it in clusters**, cuase when running gmx_MMPBSA with MPI, GROMACS's gmx_mpi can't be used. This is probably because of gmx_mpi conflicts with mpirun.)    
- **Other dependencies**: Pandas, Numpy, pyyaml, if you don't have them, please use 
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
1. **Obtain a Stable Conformation**: You will have a stable conformation of your system (e.g., the last frame or the equilibrated frame of your NPT simulation). Please prepare a PDB file with **chain information**.   
2. **Remove Water and Ions**: Please remove water and ions. We will apply a uniform standard to avoid errors.   
3. **Ensure Chain ID Information**: Ensure that chain IDs are properly assigned in the PDB file. The PDB format uses the CHAIN field to denote which molecules belong to which chain. If the CHAIN information is missing or improperly assigned, you can manually edit the PDB file to add chain IDs or use VMD or PyMOL to check and add chain IDs.  
4. **Make Sure Receptor firstly and then Antibody in your PDB file**  
