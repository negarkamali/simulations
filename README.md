# Abaqus User Element/Material for Enriched Finite Element Method

This Fortran code was developed to perform material damage detection via higher order harmonics by interfacing with the commercially available simulation tool, Abaqus. 

The enriched method, as outlined in this publication - https://www.sciencedirect.com/science/article/abs/pii/S0045794912000028 - is specifically intended for linear, elastic applications. To utilize the method for hyperelastic applications, code was written to define user elements and materials for the nonlinear elastic engine native to Abaqus. 
