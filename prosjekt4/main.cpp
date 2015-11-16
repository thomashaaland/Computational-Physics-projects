
#include <iostream>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <random>
#include <typeinfo>

using namespace arma;
using namespace std;

double Energyfind(double L, double J, mat Ising) {
    int k0, l0;
    double E = 0;
    for (int k = 0; k < L; k++) {
        for (int l = 0; l < L; l++) {
            k0 = k+1;
            l0 = l+1;
            if (k == L-1) {
                k0 = 0;
            }
            if (l == L-1) {
                l0 = 0;
            }

            E += -J*(Ising(k,l)*(Ising(k,l0)+Ising(k0,l)));

        }
    }
    return E;
}

double Magnetisation(mat Ising) {
    return accu(Ising);
}


// Had to create my own dot product since couldnt get Arma's internal to work
double dot(vec a, vec b) {
    int n = a.size();
    double value = 0;
    for (int i = 0; i < n; i++) {
        value += a(i)*b(i);
    }
    return value;
}

int main()
{
    // Declaration of variables
    int L = 2;
    imat ij;
    int MCiterations = 200000;
    double beta = 1;
    double J = 1;
    double E_initial, E_final, Delta_E;
    double w;
    vec Ebar, Esqbar, Mbar, Msqbar, Cv, Xi;
    vec r;
    vec E = zeros(MCiterations);
    vec M = zeros(MCiterations);
    vec Mabs = zeros(MCiterations);
    // The Ising grid filled with spins with values +/-1
    mat Ising = ones<mat>(L,L);

    cout << Ising << endl;

    int ltot = 7;
    int l = 0;
    Ebar = zeros(ltot);
    Esqbar = zeros(ltot);
    Mbar = zeros(ltot);
    Msqbar = zeros(ltot);
    Cv = zeros(ltot);
    Xi = zeros(ltot);

    for (int MC = 1e4; MC < 1e10; MC *= 10) {
        E = zeros(MC);
        M = zeros(MC);
        Mabs = zeros(MC);
        for (int i = 0; i < MC; i++) {
            // Finding the initial energy value
            E_initial = Energyfind(L,J,Ising);
            // Flipping some random spin position
            ij = randi<imat>(1,2, distr_param(0,L-1));
            Ising(ij(0),ij(1)) *= -1;

            // Finding the change in energy
            E_final = Energyfind(L,J,Ising);
            Delta_E = E_final-E_initial;

            // If the energychange is favourable, then flip, if not flip only with some probability
            // dependent on temperature, that is only flip if r > w
            if (Delta_E > 0) {
                w = exp(-beta*Delta_E);
                r = randu<vec>(1);
                if (r(0) > w) {
                    Ising(ij(0),ij(1)) *= -1;
                }
            }
            E(i) = Energyfind(L,J,Ising);
            M(i) = Magnetisation(Ising);
            Mabs(i) = abs(Magnetisation(Ising));
        }


        // Finding the energy and magnetication
        Ebar(l) = accu(E)/MCiterations;
        Esqbar(l) = dot(E,E)/MCiterations;
        Mbar(l) = accu(Mabs)/MCiterations;
        Msqbar(l) = dot(M,M)/MCiterations;
        Cv(l) = Esqbar(l) - Ebar(l)*Ebar(l);
        Xi(l) = Msqbar(l) - Mbar(l)*Mbar(l);
        l += 1;
    }

    // Write to a file
    ofstream outfile;
    outfile.open("Table_of_values.txt");
    outfile << "<E>= " << Ebar << "  <M>= " << Mbar << endl;
    outfile << "<E^2>= " << Esqbar << "  <M^2>= " << Msqbar << endl;
    outfile << "<Cv>=  " << Cv << "  <Xi>=  " << Xi << endl;
    outfile.close();

    return 0;
}

