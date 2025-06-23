# MatLab Implementation of EKV 2.6 by Nicola Scolari

Inthis folder you find the MatLab implementation of EKV 2.6 developed by Dr. Nicola Scolari.

The ekv.m file is the model code and takes 6 input parameters:
  1. The model parameters.
  2. The transistor gometry with d=struct('l',1e-6,'w',1e-6,'ns',1,'m',1);
  3. The gate voltage(s).
  4. The source voltage(s).
  5. The drain voltage(s).
  6. The temperature in degree Celsius.
 
The function returns 5 variables:
  1. Drain current Ids.
  2. Gate transconductance gm.
  3. Source transconductance gms.
  4. Drain transconductance gmd.
  5. A structure with all the other parameters (n, IS, mode, qI, qB, beta, ...).

The ekvint.m and ekvint.fig files are graphical interface helping for parameter extraction. You can launch it in MatLab by typing ekvint.
