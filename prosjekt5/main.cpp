#include "math/random.h"
#include "potentials/lennardjones.h"
#include "integrators/eulercromer.h"
#include "integrators/velocityverlet.h"
#include "system.h"
#include "statisticssampler.h"
#include "atom.h"
#include "io.h"
#include "unitconverter.h"
#include "unittests.h"
#include <omp.h>
#include <iostream>

using namespace std;

/*
int main(int numberOfArguments, char **argumentList)
{

    int numberOfUnitCells = 5;
    double initialTemperature = UnitConverter::temperatureFromSI(8*pow(10,4)); // measured in Kelvin
    double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

    // If a first argument is provided, it is the number of unit cells
    if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
    // If a second argument is provided, it is the initial temperature (measured in kelvin)
    if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
    // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
    if(numberOfArguments > 3) initialTemperature = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

    double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds

    cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
    cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
    cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
    cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
    cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

    System system;
    system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
    system.setPotential(new LennardJones(3.405, 119.8)); // You must insert correct parameters here
    system.setIntegrator(new VelocityVerlet());
    system.removeTotalMomentum();

    UnitTests test;
    test.test_removeTotalMomentum(&system);

    StatisticsSampler statisticsSampler;
    IO movie; // To write the state to file
    movie.open("movie.xyz");

    ofstream samples;
    samples.open("StatisticsFile.txt");


    cout << "Timestep Time Temperature KineticEnergy PotentialEnergy TotalEnergy" << endl;
    for(int timestep=0; timestep<1000; timestep++) {
        system.step(dt);
        statisticsSampler.sample(system);
        if( !(timestep % 100) ) {
            // Print the timestep every 100 timesteps
            cout << system.steps() << "      " << system.time() << "      " << UnitConverter::temperatureToSI(statisticsSampler.temperature()) << "      " << statisticsSampler.kineticEnergy() << "      "
                 << statisticsSampler.potentialEnergy() << "      " << statisticsSampler.totalEnergy() << "     " << UnitConverter::diffusionToSI(statisticsSampler.diffusionConstant()) << endl;
        }
        movie.saveState(&system);
        statisticsSampler.saveToFile(system, samples);
    }



    //movie.close();
    samples.close();

    return 0;
}*/


int main(int numberOfArguments, char **argumentList)
{
    double temp = 40000;
    double dTemp = 1000;
    for (int q = 0; q < 101; q++) {
        int numberOfUnitCells = 5;
        double initialTemperature = UnitConverter::temperatureFromSI(temp); // measured in Kelvin
        double latticeConstant = UnitConverter::lengthFromAngstroms(5.26); // measured in angstroms

        // If a first argument is provided, it is the number of unit cells
        if(numberOfArguments > 1) numberOfUnitCells = atoi(argumentList[1]);
        // If a second argument is provided, it is the initial temperature (measured in kelvin)
        if(numberOfArguments > 2) initialTemperature = UnitConverter::temperatureFromSI(atof(argumentList[2]));
        // If a third argument is provided, it is the lattice constant determining the density (measured in angstroms)
        if(numberOfArguments > 3) initialTemperature = UnitConverter::lengthFromAngstroms(atof(argumentList[3]));

        double dt = UnitConverter::timeFromSI(1e-15); // Measured in seconds

        cout << "One unit of length is " << UnitConverter::lengthToSI(1.0) << " meters" << endl;
        cout << "One unit of velocity is " << UnitConverter::velocityToSI(1.0) << " meters/second" << endl;
        cout << "One unit of time is " << UnitConverter::timeToSI(1.0) << " seconds" << endl;
        cout << "One unit of mass is " << UnitConverter::massToSI(1.0) << " kg" << endl;
        cout << "One unit of temperature is " << UnitConverter::temperatureToSI(1.0) << " K" << endl;

        System system;
        system.createFCCLattice(numberOfUnitCells, latticeConstant, initialTemperature);
        system.setPotential(new LennardJones(3.405, 119.8)); // You must insert correct parameters here
        system.setIntegrator(new VelocityVerlet());
        system.removeTotalMomentum();

        UnitTests test;
        test.test_removeTotalMomentum(&system);

        StatisticsSampler statisticsSampler;
        IO movie; // To write the state to file
        movie.open("movie.xyz");

        ofstream samples;
        samples.open("StatisticsFile" + to_string(int(temp)) + ".txt");


        cout << "Timestep Time Temperature KineticEnergy PotentialEnergy TotalEnergy" << endl;
        for(int timestep=0; timestep<1000; timestep++) {
            system.step(dt);
            statisticsSampler.sample(system);
            if( !(timestep % 100) ) {
                // Print the timestep every 100 timesteps
                cout << system.steps() << "      " << system.time() << "      " << UnitConverter::temperatureToSI(statisticsSampler.temperature()) << "      " << statisticsSampler.kineticEnergy() << "      "
                     << statisticsSampler.potentialEnergy() << "      " << statisticsSampler.totalEnergy() << "     " << statisticsSampler.diffusionConstant() << endl;
            }
            //movie.saveState(&system);
            statisticsSampler.saveToFile(system, samples);
        }
        temp += dTemp;
        samples.close();
    }




    //movie.close();


    return 0;
}
