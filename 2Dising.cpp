// Author: Wei-ting Chiu, Aug. 2019
// This Program calculate 2D Ising model, can be used to calculate
// E, M^2, M^4, and Binder ratio of function of T with different L (latice size)
// and B field.
// Wolff Cluster algorithm is applied. 
#include <cstdlib>
#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <array>
#include "Ising_func.h"

using namespace std;
// Define constant parameters.
// L: lattice length, LxL is lattice size.
// B: magnetic field.
// T_sweep: # of warm up sweeps.
// M_sweep: # of measurement sweeps.
// T_min: minimum value of T.
// T_max: maximum value of T.
// T_step: increment of T.
const int _L = 10;
const int _B = 0;
const int T_sweep = 1e5;
const int M_sweep = 3e5;
const double T_min = 0.05;
const double T_max = 4.0;
const double T_step = 0.05;

void print_result(Ising2D model, ofstream& myfile, double T) {
    
    double N = (double) (_L*_L); // lattice size
    double averageE, averageEsq, averageChi, averageChisq, C, Chi;
    double averagemsq, averagemquad;
    double binder;
    double m;
    
    averageE = model.Etot / double(M_sweep) * 10.;
    averageEsq = model.Esquare / double(M_sweep) * 10.;
    averagemsq = model.msq / double(M_sweep) * 10.;
    averagemquad = model.mquad / double(M_sweep) * 10.;
    binder = 0.5 * (3. - averagemquad / ( pow(averagemsq, 2) ) ) ;
    myfile << T << " " << averageE << " " << averagemsq << " " << averagemquad << " " << binder << "\n";
}

int main() {
    Ising2D model;
    double T;
    int i, j, k, m, n, o, p, count;
    int temp;
    double r, Nd; //Nd is for double type N.
    std::cout.precision(7);
    
    ofstream myfile;
    myfile.open ("output_L_" + std::to_string(_L) + ".txt");
    std::cout << "Start L = " << _L << endl;
    
    srand (time(NULL));
  
    for (T=T_min; T<=T_max; T=T+T_step) {
        model.T = T;
        model.L = _L;
        model.B = _B;
        model.init();
        std::cout << "Start Thermalization, T = "<< T << endl;
        model.Thermalization(T_sweep);
        std::cout << "Start Measurements, T = " << T << endl;
        model.Measurement(M_sweep, false, myfile);
        print_result(model, myfile,  T);
    }
    
    myfile.close();
    std::cout << "Finished ..." << endl;
    
    return 0;
}
