#include <iostream>

#include "Registry.h"
#include "IO/Logging.h"
#include "IO/FileInput.h"
#include "IO/CheckpointIO.h"
#include "potentials/LJ12_6.h"
#include "potentials/FENE.h"
#include "util/PhasespaceGenerator.h"
#include "potentials/Limit.h"
#include "sensors/LJ12_6_Sensor.h"
#include "sensors/FENE_Sensor.h"
#include "sensors/RDFSensor.h"
#include "thermostats/VelocityScaling.h"

static void init();
static void finalize();
static void quit(int code);

int main(int argc, char** argv) {
    init();

    if (argc <= 1) {
        std::cout << "Run program with: KoMD <input file>" << std::endl;
        quit(0);
    }

    if(!FileInput::readFile(argv[1])) quit(0);
    Registry::instance->moleculeContainer_ptr() = std::make_shared<MoleculeContainer>();
    Registry::instance->boundary_ptr() = std::make_shared<Boundary>();
    Registry::instance->simulation_ptr() = std::make_shared<Simulation>();
    Registry::instance->forceFunctors().push_back(std::make_unique<LJ12_6>());
    Registry::instance->forceFunctors().push_back(std::make_unique<FENE>());
    Registry::instance->forceFunctors().push_back(std::make_unique<Limit>());
    Registry::instance->integrators().push_back(std::make_unique<Integrator>());

    if (Registry::instance->configuration()->loadCheckpoint) CheckpointIO::loadCheckpoint();
    else PhasespaceGenerator::generate();

    Registry::instance->temperature_sensor_ptr() = std::make_shared<TemperatureSensor>();
    Registry::instance->sensors().push_back(std::dynamic_pointer_cast<Sensor>(Registry::instance->temperature_sensor()));
    if (Registry::instance->configuration()->enable_sensor_lj) Registry::instance->sensors().push_back(std::make_shared<LJ12_6_Sensor>());
    if (Registry::instance->configuration()->enable_sensor_fene) Registry::instance->sensors().push_back(std::make_shared<FENE_Sensor>());
    if (Registry::instance->configuration()->enable_sensor_rdf) Registry::instance->sensors().push_back(std::make_shared<RDFSensor>());

    Registry::instance->thermostats().push_back(std::make_unique<VelocityScaling>());

    Registry::instance->boundary()->setup();
    Registry::instance->simulation()->run();

    finalize();
    return 0;
}

static void init() {
    Registry::instance = std::make_unique<Registry>();
    Registry::instance->configuration_ptr() = std::make_shared<Configuration>();
    Registry::instance->vtkWriter_ptr() = std::make_shared<VTKWriter>();
    Log::init();
}

static void finalize() {
    Log::finalize();
}

static void quit(int code) {
    finalize();
    exit(code);
}