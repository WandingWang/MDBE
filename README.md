# antibody_test
# **Pipeline Name**

## **Overview**
This pipeline is designed for **[multiple short MD simulations for binding energy calculation]**. 
### **The whole process includes**  
**1. The initial MD：Make topolofy, Buid box, Add water and ions, Energy minimization, NVT and NTP;**   
**2. The production MD;**  
**3. Binding free energy calculation: gmx_MMPBSA;**  
**4. Mutation.**  
It automates the process from input preparation to result generation.

---

## **Installation**
### **Dependencies**
- CONDA
- Python
- VMD
- GROMACS (GPU version)
- gmx_MMPBSA
- Modeller
- [Other dependencies]

### **Installation Steps**
1. Clone the repository:
   ```bash
   git clone https://github.com/antibody_test.git
   cd antibody_test
2. Make sure you have all dependencies
3. Open the input file **[input.yaml]**, change parameters: **input setting, IMPORTANT Basic setting, gmx_mmpbsa, Modeller and run**, You need to set the input file path, conda activation path, etc., as well as the number of chains, mutation locations, etc., based on your own requirements.
4. When you make sure that all dependencies and setting are ok, run main script:
   ```bash
   main.py
   
