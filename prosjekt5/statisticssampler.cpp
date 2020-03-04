#include "system.h"
#include "statisticssampler.h"
#include "potentials/potential.h"
#include <io.h>
using namespace std;

StatisticsSampler::StatisticsSampler()
{

}

void StatisticsSampler::saveToFile(System &system, ofstream &samples)
{
    if (system.steps() == 1){samples << "IterationStep Time Temperature KineticEnergy PotentialEnergy TotalEnergy" << " " << "D" << endl;}
    samples << system.steps() << " " << system.time() << " " << UnitConverter::temperatureToSI(m_temperature) << " " << UnitConverter::energyToEv(m_kineticEnergy) << " " << UnitConverter::energyToEv(m_potentialEnergy) << " " << UnitConverter::energyToEv(m_kineticEnergy+m_potentialEnergy) << " " << m_diffusionConstant << endl;
}

void StatisticsSampler::sample(System &system)
{
    // Here you should measure different kinds of statistical properties and save it to a file.
    sampleKineticEnergy(system);
    samplePotentialEnergy(system);
    sampleTemperature(system);
    sampleDensity(system);
    sampleDiffusionConstant(system);
}

void StatisticsSampler::sampleKineticEnergy(System &system)
{
    m_kineticEnergy = 0; // Remember to reset the value from the previous timestep
    for(Atom *atom : system.atoms()) {
        m_kineticEnergy += 0.5 * atom->mass() * atom->velocity.lengthSquared();
    }
}

void StatisticsSampler::samplePotentialEnergy(System &system)
{
    m_potentialEnergy = system.potential()->potentialEnergy();
}

void StatisticsSampler::sampleTemperature(System &system)
{
    m_temperature = 2. / 3 * m_kineticEnergy / (system.atoms().size());
}

void StatisticsSampler::sampleDensity(System &system)
{

}

void StatisticsSampler::sampleDiffusionConstant(System &system)
{
    m_diffusionConstant = 0;
    for (Atom *atom : system.atoms()) {
        m_diffusionConstant += (atom->position-atom->initialPosition).lengthSquared();
    }
    m_diffusionConstant /= (6 * system.time() * system.atoms().size());
}
