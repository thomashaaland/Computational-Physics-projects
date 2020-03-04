#include "eulercromer.h"
#include "../system.h"
#include <iostream>

//using namespace std;

void EulerCromer::integrate(System *system, double dt)
{
    system->calculateForces();
    for(Atom *atom : system->atoms()) {
        atom->velocity += atom->force*dt / atom->mass();
        atom->position += atom->velocity*dt;
        system->applyPeriodicBoundaryConditions(atom->position);
    }

    //system->applyPeriodicBoundaryConditions();
}


