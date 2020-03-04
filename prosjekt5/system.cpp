#include "system.h"
#include "integrators/integrator.h"
#include "potentials/potential.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "math/random.h"
#include <iostream>

using namespace std;


System::System()
{

}

System::~System()
{
    delete m_potential;
    delete m_integrator;
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();
}

void System::applyPeriodicBoundaryConditions(vec3 &position) {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    for (int i = 0; i < 3; i++) {
        if (position(i) > 0.5 * m_systemSize[i]) {position(i) -= m_systemSize[i];}
        if (position(i) < -0.5 * m_systemSize[i]) {position(i) += m_systemSize[i];}
    }
}

void System::applyPeriodicBoundaryConditions(vec3 &position, vec3 &initialPosition) {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    for (int i = 0; i < 3; i++) {
        if (position(i) > 0.5 * m_systemSize[i]) {initialPosition(i) -= m_systemSize[i];}
        if (position(i) < -0.5 * m_systemSize[i]) {initialPosition(i) += m_systemSize[i];}
        if (position(i) > 0.5 * m_systemSize[i]) {position(i) -= m_systemSize[i];}
        if (position(i) < -0.5 * m_systemSize[i]) {position(i) += m_systemSize[i];}
    }
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    vec3 totalMomentum;
    vec3 averageMomentum;
    for(Atom *atom : m_atoms) {
        totalMomentum += atom->mass() * atom->velocity;
    }
    averageMomentum = totalMomentum / m_atoms.size();
    for(Atom *atom : m_atoms) {
        atom->velocity -= averageMomentum / atom->mass();
    }
}

void System::resetForcesOnAllAtoms() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
}

void System::createFCCLattice(int numberOfUnitCellsEachDimension, double latticeConstant, double temperature) {
    // You should implement this function properly. Right now, 100 atoms are created uniformly placed in the system of size (10, 10, 10).

    // Create FCC lattice
    int N = numberOfUnitCellsEachDimension;
    double b = latticeConstant;
    double x[4] = {0, 0.5 * b, 0,       0.5 * b};
    double y[4] = {0, 0.5 * b, 0.5 * b, 0};
    double z[4] = {0, 0,       0.5 * b, 0.5 * b};

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                for (int l = 0; l < 4; l++) {
                    Atom *atom = new Atom(UnitConverter::massFromSI(6.63352088e-26));
                    atom->position.set(b * i + x[l],b * j + y[l],b * k + z[l]);
                    atom->initialPosition = atom->position;
                    atom->resetVelocityMaxwellian(temperature);
                    m_atoms.push_back(atom);
                }
            }
        }
        setSystemSize(vec3(b * N,b * N,b * N));

    }
}

void System::calculateForces() {
    resetForcesOnAllAtoms();
    m_potential->calculateForces(this);
}

void System::step(double dt) {
    m_integrator->integrate(this, dt);
    m_steps++;
    m_time += dt;
}
