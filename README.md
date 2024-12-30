# antibody_test
# **Pipeline Name**

## **Overview**
This pipeline is designed for **[multiple short MD simulations for binding energy calculation]**. 
### **The whole process includes**
**1. The initial MD：Make topolofy, Buid box, Add water and ions, Energy minimization, NVT and NTP; The prodcution MD;** 
**2. The production MD;**
**3. Binding free energy calculation: gmx_MMPBSA;**
**4. Mutation.**
It automates the process from input preparation to result generation.

---

## **Installation**
### **Dependencies**
- Python >= 3.x
- GROMACS >= 2024.4
- gmx_MMPBSA >= 1.x
- [Other dependencies]

### **Installation Steps**
1. Clone the repository:
   ```bash
   git clone https://github.com/your-repo-name.git
   cd your-repo-name
