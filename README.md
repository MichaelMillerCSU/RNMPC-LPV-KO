## There are 2 main examples: the Van der Pol system and Duffing system are shown in the paper, and an addtional example of the normal oscillator system.
### The file VDP.m, Duffing.m and Normal_Oscillator are corresponding to the different system with Min-Max K-RMPC. The function of the dynamics should be changed if running KMPC.m.
### The folder with 'Quasi' contains the Van der Pol system and Duffing system applied with the Quasi-Min-Max K-RMPC.


#Erratum: VDP2.m is the simulation file of CDC 2024 version, and the figure 1(b) is obtained by mistakes such that the negative sign was missing out. However, the missing negative sign DO NOT influence the performance of the control signal, since the negative sign would be canceled in this equality : $P - (A + BF)^{T}P(A+BF) = P - (-A - BF)^{T}P(- A - BF) \ge Q + F^{T}RF$.
