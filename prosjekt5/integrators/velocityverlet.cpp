#include "velocityverlet.h"
#include "../system.h"
#include "../atom.h"
#include<iostream>

void VelocityVerlet::integrate(System *system, double dt)
{
    if (system->steps() == 0) {
        system->calculateForces();
        std::cout << "Calculating initial forces!" << std::endl;
    }

    for(Atom *atom : system->atoms()) {
        atom->velocity += 0.5 * atom->force*dt / atom->mass();
        atom->position += atom->velocity*dt;
        system->applyPeriodicBoundaryConditions(atom->position, atom->initialPosition);

    }
    system->calculateForces();
    for(Atom *atom : system->atoms()) {
        atom->velocity += 0.5 * atom->force*dt / atom->mass();
    }
}
