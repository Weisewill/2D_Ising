#include <cstdlib>
#include <iostream>
#include <cmath>
#include <time.h>
#include <fstream>
#include <array>

using namespace std;

#define T_c 2.269

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


void Ising2D::growCluster(int i, int j, int clusterspin) {
    bond[i][j] = 1;
    S[i][j] = - S[i][j];
    
    neighbors N;
    N.iloc = i; N.jloc = j;
    N.find_neighbors(L);
    
    if(bond[i][N.up]==0)
    judgeCluster(i,N.up,clusterspin);
    if(bond[N.right][j]==0)
    judgeCluster(N.right,j,clusterspin);
    if(bond[i][N.down]==0)
    judgeCluster(i,N.down,clusterspin);
    if(bond[N.left][j]==0)
    judgeCluster(N.left,j,clusterspin);
}

void Ising2D::judgeCluster(int i, int j, int clusterspin){
    double r;
    double addProb = 1. - exp(-2./T);
    if(S[i][j] == clusterspin) {
        r = rand();
        r = r / RAND_MAX;
        if(addProb > r)
        growCluster(i,j,clusterspin);
    }
}

double Ising2D::calculate_energy(int i, int j){
    
    neighbors N;
    N.iloc = i; N.jloc = j;
    N.find_neighbors(L);
    
    //Calculate energy of 4 adjecent bonds of S[i][j]
    E = - S[i][j]*( S[N.right][j] + S[i][N.up] + S[i][N.down] + S[N.left][j]) - B*S[i][j];
    return E;
}

void neighbors::find_neighbors(int L){
    // This function find the neighbors of a site at (i, j)
    // w/ periodic boundary condition.
    right = iloc + 1; up = jloc + 1;
    left = iloc - 1; down = jloc - 1;
    if ( right == L && up == L ){
        right = 0; up = 0;
    }
    else if ( right == L && up != L ){
        right = 0;
    }
    else if ( right != L && up == L ){
        up = 0;
    }
    if ( iloc == 0 && jloc == 0 ){
        left = L - 1; down = L - 1;
    }
    else if ( iloc == 0 && jloc != 0 ){
        left = L - 1;
    }
    else if ( iloc != 0 && jloc == 0 ){
        down = L - 1;
    }
}

Ising2D::Ising2D() {
    int i, j;
    L = 10;
    B = 0;
    T = 2.0;
    E = 0;
    Etot = 0;
    Esquare = 0;
    Chitot = 0;
    sqChi = 0;
    msq = 0;
    mquad = 0;
    bond = new int*[L];
    for(int i = 0; i < L; ++i)
        bond[i] = new int[L];
    S = new int*[L];
    for(int i = 0; i < L; ++i)
        S[i] = new int[L];
    for (i=0;i<L;i++) {
        for (j=0;j<L;j++) {
            S[i][j] = rand() % 2;
            bond[i][j] = 0;
            if (S[i][j] == 0) {
                S[i][j] = - 1;
            }
        }
    }
}

void Ising2D::init() {
    int i, j;
    E = 0;
    Etot = 0;
    Esquare = 0;
    Chitot = 0;
    sqChi = 0;
    msq = 0;
    mquad = 0;
    bond = new int*[L];
    for(int i = 0; i < L; ++i)
    bond[i] = new int[L];
    S = new int*[L];
    for(int i = 0; i < L; ++i)
    S[i] = new int[L];
    for (i=0;i<L;i++) {
        for (j=0;j<L;j++) {
            S[i][j] = rand() % 2;
            bond[i][j] = 0;
            if (S[i][j] == 0) {
                S[i][j] = - 1;
            }
        }
    }
}

void Ising2D::metropolis(int i, int j){
    double Eold, Enew, deltaE, r, Pjudge;
    Eold = calculate_energy(i,j);
    S[i][j] = - S[i][j];
    Enew = calculate_energy(i,j);
    deltaE = Enew - Eold;
    r = rand();
    r = r / RAND_MAX;
    Pjudge = exp(-deltaE/T) / ( 1. + exp(-deltaE/T) ) ;
    if (Pjudge >= r) {
    }
    else {
        S[i][j] = - S[i][j];
    }
}

void Ising2D::oneMonteCarloStep() {
    int i, j, k;
    //srand (time(NULL));
    for (i=0;i<L;i++) {
        for (j=0;j<L;j++) {
            metropolis(i, j);
        }
    }
    // Wolff algorithm
    // Reset
    for (i=0;i<L;i++)
    for (j=0;j<L;j++)
    bond[i][j] = 0;
    // Choose a random site.
    i = rand() % L;
    j = rand() % L;
    growCluster(i,j,S[i][j]);
}

void Ising2D::Thermalization(int T_sweep){
    // Thermalization
    int k;
    for (k=0;k<T_sweep;k++) {
        oneMonteCarloStep();
    }
}

void Ising2D::Measurement(int M_sweep, bool print_config, ofstream& myfile) {
    int count, i, j, k;
    double Nd = (double) (L*L);
    neighbors N;
    double Chitemp = 0, Etottemp = 0, msqtemp = 0, mquadtemp = 0, mtemp;
    for (k=0;k<M_sweep;k++) {
        oneMonteCarloStep();
        if ( k % 10 == 0) {
            for (i=0;i<L;i++) {
                for (j=0;j<L;j++) {
                    N.iloc = i; N.jloc = j;
                    N.find_neighbors(L);
                    Etottemp += - S[i][j]*S[N.right][j] - S[i][j]*S[i][N.up] - B*S[i][j];
                    Chitemp += S[i][j];
                }
            }
            Etot += Etottemp / Nd;
            Esquare += Etottemp * Etottemp / Nd;
            Chitemp /= Nd;
            msq += Chitemp * Chitemp ;
            mquad += pow(Chitemp, 4) ;
            //Chitot = Chitot + abs(1./Nd*Chitemp);
            //sqChi = sqChi + 1. /Nd/Nd*(Chitemp*Chitemp);
            Chitemp = 0;
            Etottemp = 0;
            msqtemp = 0;
            mquadtemp = 0;
            if(print_config == true)
            print_lattice(myfile);
        }
    }
}

void Ising2D::print_lattice(ofstream& myfile) {
    int i, j;
    int label;
    label = ( T < T_c) ? 1 : 0; // 1 for ordered phase; 0 for non-ordered phase.
    myfile << label << '\t';
    for(i=0; i<L; i++)
    for(j=0; j<L; j++)
    myfile << S[i][j] << '\t';
    myfile << endl;
}
