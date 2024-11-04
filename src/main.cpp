#include <iostream>

#include "Registry.h"
#include "IO/Logging.h"
#include "IO/FileInput.h"
#include "IO/CheckpointIO.h"

#include "potentials/LJ12_6.h"
#include "potentials/FENE.h"
#include "potentials/Limit.h"
#include "potentials/ATM.h"
//#include "potentials/Grav.h"

#include "plugins/ADR/AdResS.h"

#include "util/PhasespaceGenerator.h"

#include "sensors/LJ12_6_Sensor.h"
#include "sensors/FENE_Sensor.h"
#include "sensors/RDFSensor.h"

#include "thermostats/VelocityScaling.h"
#include "container/LinkedCells.h"

#include "Kokkos_Core.hpp"

static void init(int argc, char** argv);
static void finalize();
static void quit(int code);

int main(int argc, char** argv) {
    init(argc, argv);

    if (argc <= 1) {
        std::cout << "Run program with: KoMD <input file>" << std::endl;
        quit(0);
    }

    if(!FileInput::readFile(argv[1])) quit(0);
    if (Registry::instance->configuration()->enable_one_cell) throw std::runtime_error("not supported anymore");
    else Registry::instance->moleculeContainer_ptr() = std::make_shared<LinkedCells>();
    Registry::instance->simulation_ptr() = std::make_shared<Simulation>();
    if (Registry::instance->configuration()->enable_one_cell) throw std::runtime_error("not supported anymore");
    else Registry::instance->forceFunctors().push_back(std::make_unique<LJ12_6>());
    Registry::instance->forceFunctors().push_back(std::make_unique<FENE>());
    Registry::instance->forceFunctors().push_back(std::make_unique<Limit>());
    //Registry::instance->forceFunctors().push_back(std::make_unique<Grav>());
    if (Registry::instance->configuration()->enable_3b) Registry::instance->forceFunctors3b().push_back(std::make_unique<ATM>());
    Registry::instance->integrators().push_back(std::make_unique<Integrator>());

    CheckpointIO::loadCheckpoint();
    PhasespaceGenerator::generate();

    Registry::instance->temperature_sensor_ptr() = std::make_shared<TemperatureSensor>();
    Registry::instance->sensors().push_back(std::dynamic_pointer_cast<Sensor>(Registry::instance->temperature_sensor()));
    if (Registry::instance->configuration()->enable_sensor_lj) Registry::instance->potential_sensor_ptr() = std::make_shared<LJ12_6_Sensor>();
    if (Registry::instance->configuration()->enable_sensor_lj) Registry::instance->sensors().push_back(std::dynamic_pointer_cast<Sensor>(Registry::instance->potential_sensor()));
    if (Registry::instance->configuration()->enable_sensor_fene) Registry::instance->sensors().push_back(std::make_shared<FENE_Sensor>());
    if (Registry::instance->configuration()->enable_sensor_rdf) Registry::instance->sensors().push_back(std::make_shared<RDFSensor>());
    if (Registry::instance->configuration()->enable_sensor_disp) Registry::instance->displacement_sensor_ptr() = std::make_shared<DisplacementSensor>();
    if (Registry::instance->configuration()->enable_sensor_disp) Registry::instance->sensors().push_back(std::dynamic_pointer_cast<Sensor>(Registry::instance->displacement_sensor()));

    // Plugins
    if (Registry::instance->configuration()->ADR_enable) Registry::instance->plugins().push_back(std::make_unique<AdResS>());

    Registry::instance->thermostats().push_back(std::make_unique<VelocityScaling>());

    Registry::instance->moleculeContainer()->init();
    Registry::instance->simulation()->run();

    finalize();
    return 0;
}

static void init(int argc, char** argv) {
    Kokkos::initialize(argc, argv);
    Registry::instance = std::make_unique<Registry>();
    Registry::instance->configuration_ptr() = std::make_shared<Configuration>();
    Registry::instance->vtkWriter_ptr() = std::make_shared<VTKWriter>();
    Log::init();
}

static void finalize() {
    Registry::instance = nullptr;
    Log::finalize();
    Kokkos::finalize();
}

static void quit(int code) {
    finalize();
    exit(code);
}