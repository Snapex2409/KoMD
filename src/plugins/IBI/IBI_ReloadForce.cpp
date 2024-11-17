//
// Created by alex on 11/11/24.
//

#include "IBI_ReloadForce.h"

#include "ExcluisionLJ.h"
#include "ExclusionATM.h"
#include "IBI_ForceFunctor.h"
#include "Registry.h"

IBI_ReloadForce::IBI_ReloadForce() : Plugin("IBI-Reload") {
    auto config = Registry::instance->configuration();

    if (config->IBI_enable) throw std::runtime_error("Cannot have IBI and IBI_Reload at the same time!");

    FunctionPL force{config->IBI_bins, 1e+6, 0};
    FunctionPL pot{config->IBI_bins, 1e+6, 0};
    force.read(config->IBI_reload_fpath);
    pot.read(config->IBI_reload_ppath);

    const auto volume = (config->IBI_exclusion_high - config->IBI_exclusion_low).product();
    auto& force_functors = Registry::instance->forceFunctors();
    force_functors.clear();
    auto& force_functors3b = Registry::instance->forceFunctors3b();
    force_functors3b.clear();

    if (volume == 0) {
        auto ptr = std::make_shared<IBI_ForceFunctor<IBI_Default>>();
        ptr->getForceFunction().setXValues(force.getXValues());
        ptr->getForceFunction().setYValues(force.getYValues());
        ptr->getForceFunction().setLowerDefault(force.getLowerDefault());
        ptr->getForceFunction().setUpperDefault(force.getUpperDefault());

        ptr->getPotentialFunction().setXValues(pot.getXValues());
        ptr->getPotentialFunction().setYValues(pot.getYValues());
        ptr->getPotentialFunction().setLowerDefault(pot.getLowerDefault());
        ptr->getPotentialFunction().setUpperDefault(pot.getUpperDefault());

        force_functors.push_back(ptr);
        Registry::instance->moleculeContainer()->disable3B();
    }
    else {
        if (!config->enable_3b) throw std::runtime_error("Three body potential must be configured, when not using as full replacement!");
        auto ptr = std::make_shared<IBI_ForceFunctor<IBI_Reload>>();
        ptr->getForceFunction().setXValues(force.getXValues());
        ptr->getForceFunction().setYValues(force.getYValues());
        ptr->getForceFunction().setLowerDefault(force.getLowerDefault());
        ptr->getForceFunction().setUpperDefault(force.getUpperDefault());

        ptr->getPotentialFunction().setXValues(pot.getXValues());
        ptr->getPotentialFunction().setYValues(pot.getYValues());
        ptr->getPotentialFunction().setLowerDefault(pot.getLowerDefault());
        ptr->getPotentialFunction().setUpperDefault(pot.getUpperDefault());

        force_functors.push_back(ptr);
        force_functors.push_back(std::make_shared<ExcluisionLJ>());
        force_functors3b.push_back(std::make_shared<ExclusionATM>());
        Registry::instance->moleculeContainer()->enable3B();
    }

}
