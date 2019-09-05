# 2D_Ising
2D Ising model on a square lattice.

## Set constants in 2Dising.cpp.
L: lattice length, LxL is lattice size.  
B: magnetic field.  
T_sweep: # of warm up sweeps.  
M_sweep: # of measurement sweeps.  
T_min: minimum value of T.  
T_max: maximum value of T.  
T_step: increment of T.  

### Default values:
const int _L = 10;  
const int _B = 0;  
const int T_sweep = 1e5;  
const int M_sweep = 3e5;  
const double T_min = 0.05;  
const double T_max = 4.0;  
const double T_step = 0.05;  

## How to run:

[1] make  
    Generate binary file: ising.  
[2] ./ising

One can plot Binder ratio v.s. T, and find that the curves from different sizes of lattice cross at T = 2.26 ~ 2.27, where the analytic critical T is 2.269.
![Binder ratios v.s. T](/2D_Ising_binder.png)
