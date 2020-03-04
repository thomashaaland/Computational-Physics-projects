#include "unittests.h"

void UnitTests::test_removeTotalMomentum(System *system) {
    double epsilon = 1e-10;
    vec3 totalMomentum;
    for(Atom *atom : system->atoms()) {
        totalMomentum += atom->mass() * atom->velocity;
    }
    if (totalMomentum.length() < epsilon) {cout << "The total momentum has been successfully calbrated to zero" << endl;}
    else {cout << "ERROR: The total momentum has not been successfully recalibated to zero";}
}
