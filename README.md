# TDA_script
Contains main script to TDA application as utilized by Argiris et al. (2023) Topological Data Analysis Predicts Longitudinal Behavioral Change of Fluid Reasoning in the RANN Cohort


This script depends on a TDA tutorial developed by Centeno et al. (2022)
Centeno, E.G.Z., Moreni, G., Vriend, C., Douw, L., & Santos, F.A.N. (2022). A hands-on tutorial on network and topological neuroscience. Brain Structure and Function, 227(3), 741-762. 

TDA_ex.py is the main script and requires functions graph_density.py and Betti_k.py
  (A modified version of the graph_density.py function can be found after line 45 of the script)
data.mat contains the .mat MATLAB file of the Fisher's Z values of the connectivity matrices reshaped as vectors (edges x subjects x tasks)
