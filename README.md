## There are 2 main examples: the Van der Pol system and Duffing system are shown in the paper, and an addtional example of the normal oscillator system.
### The file VDP.m, Duffing.m and Normal_Oscillator are corresponding to the different system with Min-Max K-RMPC. The function of the dynamics should be changed if running KMPC.m.
### The folder with 'Quasi' contains the Van der Pol system and Duffing system applied with the Quasi-Min-Max K-RMPC.

---

# Corrigendum  
**Correction to Fig. 1(b) in the Paper**  
*"Robust Nonlinear Model Predictive Control Based on the LPV System with Polytopic Uncertainty Identified Using Koopman Operator"* (CDC 2024)  

---

## Corrected Figure  
![Corrected Fig 1(b)](https://github.com/MichaelMillerCSU/RNMPC-LPV-KO/blob/main/Corrigendum_CDC2024_Paper/Corrigendum_of_Fig_1_b.png)  

---

## Details of the Correction  
This is a correction to **Fig. 1(b)** in the above-mentioned paper.  
[The reason for this correction is that, due to a typo in the source code, the newly added vertices were not reversed. As a result, vertices 3 and 4 were identical to vertices 1 and 2. The corrected figure reflects the validation results after flipping vertices 3 and 4 by adding a negative sign, leading to improved performance.It is worth noting that whether the vertices are flipped or not does not affect the solution of the controller, as the negative sign is directly eliminated in the LMI formulation.]


---

Thank you for your understanding.
