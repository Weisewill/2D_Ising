#include <cstdlib>
#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <array>

using namespace std;

const double T_c = 2.269;

// 2D Ising model with J = 1, k_b = 1;
class Ising2D {
    public:
        Ising2D();
        double T;
        int L;
        double B;
        int **bond, **S;
        double E;
        double Etot;
        double Esquare;
        double Chitot;
        double sqChi;
        double msq;
        double mquad;
        double calculate_energy(int, int);
        double calculate_total_energy();
        void init();
        void metropolis(int, int);
        void oneMonteCarloStep();
        void Thermalization(int);
        void Measurement(int, bool, ofstream&);
        void print_lattice(ofstream&);
        // Functions for Wolff algorithm.
        void growCluster(int, int, int);
        void judgeCluster(int, int, int);
};


class neighbors {
    public:
        int iloc, jloc;
        int right, up, left, down;
        void find_neighbors(int);
};

