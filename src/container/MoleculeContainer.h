//
// Created by alex on 7/31/24.
//

#ifndef KOMD_MOLECULECONTAINER_H
#define KOMD_MOLECULECONTAINER_H

#include <vector>
#include <memory>

#include "Vec3D.h"
#include "Cell.h"
#include "molecule/Molecule.h"

class MoleculeContainer {
public:
    MoleculeContainer();
    virtual ~MoleculeContainer() = default;

    virtual void init() = 0;
    virtual void addMolecule(const Molecule& molecule) = 0;
    virtual void updateContainer() = 0;

    uint64_t getNumMolecules() { return p_molecule_count; }

    virtual void getCenterOfMassPositions(SOA::vec_t<math::d3>& buffer) = 0;
    virtual void writeSOA2AOS() = 0;

    enum IteratorType {
        MOLECULE,
        SITE
    };

    enum IteratorRegion {
        DOMAIN,
        DOMAIN_HALO
    };

    class Iterator {
    public:
        virtual ~Iterator() = default;
        virtual void operator++() = 0;
        virtual bool isValid() const = 0;

        virtual math::d3& f() = 0;
        virtual math::d3& r() = 0;
        virtual math::d3& v() = 0;
        virtual double epsilon() = 0;
        virtual double sigma() = 0;
        virtual double mass() = 0;
        virtual uint64_t ID() = 0;
        virtual Molecule& molecule() = 0;
    };

    class CellIterator {
    public:
        virtual ~CellIterator() = default;
        virtual void operator++() = 0;
        virtual bool isValid() const = 0;

        virtual Cell& cell() = 0;
    };

    class CellPairIterator {
    public:
        virtual ~CellPairIterator() = default;
        virtual void operator++() = 0;
        virtual bool isValid() const = 0;
        virtual bool colorSwitched() = 0;

        virtual Cell& cell0() = 0;
        virtual Cell& cell1() = 0;
    };

    /**
     * Only to be used for IO functionality on CPU side
     * */
    virtual std::unique_ptr<Iterator> iterator(const IteratorType type, const IteratorRegion region) = 0;
    /**
     * Only to be used on CPU side
     * */
    virtual std::unique_ptr<CellIterator> iteratorCell(const IteratorRegion region) = 0;
    /**
     * Only to be used on CPU side
     * */
    virtual std::unique_ptr<CellPairIterator> iteratorC08() = 0;
private:
    virtual void constructSOAs() = 0;
    virtual void clearForces() = 0;
protected:
    /// total number of molecules
    uint64_t p_molecule_count;
};


#endif //KOMD_MOLECULECONTAINER_H
