//
// Created by alex on 7/31/24.
//

#include "ForceFunctor.h"
#include "Registry.h"

ForceFunctor::ForceFunctor() : p_soa(Registry::instance->moleculeContainer()->getSOA()) {}

void ForceFunctor::operator()() {
    auto container = Registry::instance->moleculeContainer();
    handlePairList(container->getPairList());
    Kokkos::fence("Force Functor - Pair List");
}
