#ifndef UNITTESTS_H
#define UNITTESTS_H
#include "math/random.h"
#include "potentials/lennardjones.h"
#include "integrators/eulercromer.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include <iostream>

using namespace std;

class UnitTests
{
public:
    void test_removeTotalMomentum(System *system);

//signals:

//public slots:
};

#endif // UNITTESTS_H
