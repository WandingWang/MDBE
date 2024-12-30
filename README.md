# antibody_test

## **Overview**
This pipeline is designed for **Antibody design**. 
### **The whole process includes**  
**1. The initial MD：Make topolofy, Buid box, Add water and ions, Energy minimization, NVT and NTP;**   
**2. The production MD;**  
**3. Binding free energy calculation: gmx_MMPBSA;**  
**4. Mutation.**  
It automates the process from input preparation to result generation.

---

## **Installation**
### **Dependencies**
- VMD
- GROMACS (GPU version)
- gmx_MMPBSA: https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/installation/  
- Modeller: https://salilab.org/modeller/download_installation.html  
- Other dependencies: requirements.txt

### **Installation Steps**
1. Clone the repository:
   ```bash
   git clone https://github.com/antibody_test.git
   cd antibody_test
2. Make sure you have all dependencies: check conda, python, vmd, gromacs, gmx_mmpbsa and modeller;
3. Make sure you have all Python packages:
   ```bash
   pip install -r requirements.txt 
4. Open the input file **[input.yaml]**, change parameters: **input setting, IMPORTANT Basic setting, gmx_mmpbsa, Modeller and run**, You need to set the input file path, conda activation path, etc., as well as the number of chains, mutation locations, etc., based on your own requirements.
5. When you make sure that all dependencies and setting are ok, run main script:
   ```bash
   python AutoBindingMD.py
   
