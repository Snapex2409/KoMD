//
// Created by alex on 7/31/24.
//

#ifndef KOMD_FORCEFUNCTOR_H
#define KOMD_FORCEFUNCTOR_H

class SOA;
#include "math/Array.h"
#include "container/PairList.h"

/**
 * General force function interface
 * */
class ForceFunctor {
public:
    ForceFunctor();
    virtual ~ForceFunctor() = default;

    /**
     * This should compute some sort of forces.
     * */
    void operator()();

protected:
    /**
     * Handles all possible interactions
     * */
    virtual void handlePairList(PairList& pairList) = 0;

    /// Reference to MoleculeContainer::soa
    SOA& p_soa;
};


#endif //KOMD_FORCEFUNCTOR_H
