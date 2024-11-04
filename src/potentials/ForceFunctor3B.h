//
// Created by alex on 11/4/24.
//

#ifndef FORCEFUNCTOR3B_H
#define FORCEFUNCTOR3B_H

class SOA;
#include "container/TripleList.h"

/**
 * 3 body force function interface
 */
class ForceFunctor3B {
public:
    ForceFunctor3B();
    virtual ~ForceFunctor3B() = default;

    /**
     * This should compute some sort of forces.
     * */
    void operator()();

protected:

    /**
     * Handles all possible interactions
     * @param tripleList active triple list from particle container
     */
    virtual void handleTripleList(TripleList& tripleList) = 0;

    /// Reference to MoleculeContainer::soa
    SOA& p_soa;
};



#endif //FORCEFUNCTOR3B_H
