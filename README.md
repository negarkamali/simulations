# Abaqus User Elements/Material for Enriched Finite Element Method

This Fortran code is for an enriched finite element method, Abaqus user subroutine. Using both Lagrangian and harmonic functions to discretize the solution field, more degrees of freedom are added, which improves the accuracy and efficiency for solving wave propogation problems. This method can then be used to perform material damage detection via higher order harmonics. 

The enriched method, as outlined in this publication - https://www.sciencedirect.com/science/article/abs/pii/S0045794912000028 - is specifically intended for linear, elastic wave propogation applications. To utilize the method for hyperelastic applications, code for the UELMAT user defined subroutine was written to define user elements and materials for the nonlinear engine native to Abaqus. 
