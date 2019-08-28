# 2D_Ising
2D Ising model on a square lattice.

## How to run:

[1] Set parameters at the header of 2DIsing.cpp.

[2] Compile with:
g++ 2DIsing.cpp -I ./ -o a.out

[3] Run with:
./a.out

[4] Generate output file:
output_L_{}.txt

One can plot Binder ratio v.s. T, and find that the curves from different sizes of lattice cross at T = 2.26 ~ 2.27.
![Binder ratios v.s. T](/2D_Ising_binder.png)
