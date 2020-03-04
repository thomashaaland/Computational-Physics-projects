#include "lennardjones.h"
#include<iostream>
#include <omp.h>
#include <math.h>
#include <time.h>

LennardJones::LennardJones(double sigma, double epsilon) :
    m_sigma(sigma),
    m_epsilon(epsilon)
{

}

void LennardJones::calculateForces(System *system)
{

    m_potentialEnergy = 0; // Remember to compute this in the loop
    double s, sigma12, sigma6, rIJSquare, uIJ, uIJFirstArgument, uIJSecondArgument, prefix;
    vector<Atom *> atom = system->atoms();
    vec3 rI, rJ, rIJ, fIJ;
    vec3 zerovec;
    zerovec[0] = 0; zerovec[1] = 0; zerovec[2] = 0;
    s = atom.size();
    sigma12 = pow(m_sigma,12);
    sigma6 = pow(m_sigma,6);
    prefix = 4 * m_epsilon;
    vector<vec3> totalFi(atom.size());

//    #pragma omp parallel for

    for(int i = 0; i < s; i++) {

        rI = atom[i]->position;

        for(int j = i+1; j < s; j++) {

            rJ = atom[j]->position;
            rIJ = rI-rJ;
            system->applyPeriodicBoundaryConditions(rIJ, zerovec);
            rIJSquare = rIJ.lengthSquared();

            // The potential
            uIJFirstArgument = sigma12/pow(rIJSquare,6);
            uIJSecondArgument = sigma6/pow(rIJSquare,3);
            uIJ = prefix * (uIJFirstArgument - uIJSecondArgument);

            m_potentialEnergy += uIJ;

            // The forces
            fIJ = prefix * 6 * (2 * uIJFirstArgument - uIJSecondArgument)/rIJSquare*rIJ;

            totalFi[i] += fIJ;
            totalFi[j] -= fIJ;
        }
        atom[i]->force = totalFi[i];
    }
}

// A slower method
/*
    for(Atom *atomi : system->atoms()){
        for(Atom *atomj : system->atoms()){

            if(atomi != atomj){
                rIJ = atomi->position - atomj->position;
                system->applyPeriodicBoundaryConditions(rIJ);
                rIJSquare = rIJ.lengthSquared();
                uIJFirstArgument = sigma12/pow(rIJSquare,6);
                uIJSecondArgument = sigma6/pow(rIJSquare,3);
                uIJ = 0.5 * prefix * (uIJFirstArgument - uIJSecondArgument);
                m_potentialEnergy += uIJ;

                atomi->force += (6 * prefix * (2 * uIJFirstArgument - uIJSecondArgument)/rIJSquare)*rIJ;
            }
        }

    }*/
