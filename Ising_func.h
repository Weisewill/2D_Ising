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
        double T_;
        int L_;
        double B_;
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
        int iloc_, jloc_;
        int right_, up_, left_, down_;
        void find_neighbors(int);
};


void Ising2D::growCluster(int i, int j, int clusterspin) {
    bond[i][j] = 1;
    S[i][j] = - S[i][j];
    
    neighbors N;
    N.iloc_ = i; N.jloc_ = j;
    N.find_neighbors(L_);
    
    if(bond[i][N.up_]==0)
    judgeCluster(i,N.up_,clusterspin);
    if(bond[N.right_][j]==0)
    judgeCluster(N.right_,j,clusterspin);
    if(bond[i][N.down_]==0)
    judgeCluster(i,N.down_,clusterspin);
    if(bond[N.left_][j]==0)
    judgeCluster(N.left_,j,clusterspin);
}

void Ising2D::judgeCluster(int i, int j, int clusterspin){
    double r;
    double addProb = 1. - exp(-2./T_);
    if(S[i][j] == clusterspin) {
        r = rand();
        r = r / RAND_MAX;
        if(addProb > r)
        growCluster(i,j,clusterspin);
    }
}

double Ising2D::calculate_energy(int i, int j){
    
    neighbors N;
    N.iloc_ = i; N.jloc_ = j;
    N.find_neighbors(L_);
    
    //Calculate energy of 4 adjecent bonds of S[i][j]
    E = - S[i][j]*( S[N.right_][j] + S[i][N.up_] + S[i][N.down_] + S[N.left_][j]) - B_*S[i][j];
    return E;
}

void neighbors::find_neighbors(int L_){
    // This function find the neighbors of a site at (i, j)
    // w/ periodic boundary condition.
    right_ = iloc_ + 1; up_ = jloc_ + 1;
    left_ = iloc_ - 1; down_ = jloc_ - 1;
    if ( right_ == L_ && up_ == L_ ){
        right_ = 0; up_ = 0;
    }
    else if ( right_ == L_ && up_ != L_ ){
        right_ = 0;
    }
    else if ( right_ != L_ && up_ == L_ ){
        up_ = 0;
    }
    if ( iloc_ == 0 && jloc_ == 0 ){
        left_ = L_ - 1; down_ = L_ - 1;
    }
    else if ( iloc_ == 0 && jloc_ != 0 ){
        left_ = L_ - 1;
    }
    else if ( iloc_ != 0 && jloc_ == 0 ){
        down_ = L_ - 1;
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
    bond = new int*[L_];
    for(int i = 0; i < L_; ++i)
        bond[i] = new int[L_];
    S = new int*[L_];
    for(int i = 0; i < L_; ++i)
        S[i] = new int[L_];
    for (i=0;i<L_;i++) {
        for (j=0;j<L_;j++) {
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
    Pjudge = exp(-deltaE/T_) / ( 1. + exp(-deltaE/T_) ) ;
    if (Pjudge >= r) {
    }
    else {
        S[i][j] = - S[i][j];
    }
}

void Ising2D::oneMonteCarloStep() {
    int i, j, k;
    //srand (time(NULL));
    for (i=0;i<L_;i++) {
        for (j=0;j<L_;j++) {
            metropolis(i, j);
        }
    }
    // Wolff algorithm
    // Reset
    for (i=0;i<L_;i++)
    for (j=0;j<L_;j++)
    bond[i][j] = 0;
    // Choose a random site.
    i = rand() % L_;
    j = rand() % L_;
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
    double Nd = (double) (L_*L_);
    neighbors N;
    double Chitemp = 0, Etottemp = 0, msqtemp = 0, mquadtemp = 0, mtemp;
    for (k=0;k<M_sweep;k++) {
        oneMonteCarloStep();
        if ( k % 10 == 0) {
            for (i=0;i<L_;i++) {
                for (j=0;j<L_;j++) {
                    N.iloc_ = i; N.jloc_ = j;
                    N.find_neighbors(L_);
                    Etottemp += - S[i][j]*S[N.right_][j] - S[i][j]*S[i][N.up_] - B_*S[i][j];
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
    label = ( T_ < T_c) ? 1 : 0; // 1 for ordered phase; 0 for non-ordered phase.
    myfile << label << '\t';
    for(i=0; i<L_; i++)
    for(j=0; j<L_; j++)
    myfile << S[i][j] << '\t';
    myfile << endl;
}
