//
// Created by alex on 11/4/24.
//

#include "ForceFunctor3B.h"
#include "Registry.h"

ForceFunctor3B::ForceFunctor3B() : p_soa(Registry::instance->moleculeContainer()->getSOA()) { }

void ForceFunctor3B::operator()() {
    auto container = Registry::instance->moleculeContainer();
    handleTripleList(container->getTripleList());
    Kokkos::fence("Force Functor - Triple List");
}
