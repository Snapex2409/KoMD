//
// Created by alex on 7/31/24.
//

#ifndef KOMD_MOLECULECONTAINER_H
#define KOMD_MOLECULECONTAINER_H

#include <vector>
#include <memory>

#include "Vec3D.h"
#include "Cell.h"
#include "SOA.h"
#include "molecule/Molecule.h"
#include "util/Kokkos_Wrapper.h"
#include "PairList.h"

class MoleculeContainer {
public:
    MoleculeContainer();
    virtual ~MoleculeContainer() = default;

    virtual void init() = 0;
    void addMolecule(const Molecule& molecule);
    virtual void updateContainer() = 0;

    uint64_t getNumMolecules() { return p_molecule_count; }

    void getCenterOfMassPositions(KW::vec_t<math::d3>& buffer);

    void writeSOA2AOS();
    void constructSOAs();

    KW::vec_t<math::d3> getCOM() { return p_com; }
    SOA& getSOA() {return p_soa; }
    PairList& getPairList() { return p_pair_list; }

    enum IteratorType {
        MOLECULE,
        SITE
    };

    class Iterator {
    public:
        Iterator(bool only_mol, std::vector<Molecule>& molecules, SOA& soa);
        ~Iterator() = default;
        void operator++();
        bool isValid() const;

        math::d3& f();
        math::d3& r();
        math::d3& v();
        double epsilon();
        double sigma();
        double mass();
        uint64_t ID();
        Molecule& molecule();
    private:
        const bool m_only_molecule;
        uint64_t m_site_idx;
        uint64_t m_molecule_idx;
        uint64_t m_visited_sites;
        std::vector<Molecule>& m_molecules;
        SOA& m_soa;
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
        virtual math::d3 getCell0Shift() { return {0, 0, 0}; }
        virtual math::d3 getCell1Shift() { return {0, 0, 0}; }
    };

    /**
     * Only to be used for IO functionality on CPU side
     * */
    Iterator iterator(IteratorType type);
    /**
     * Only to be used on CPU side
     * */
    virtual std::unique_ptr<CellIterator> iteratorCell() = 0;
    /**
     * Only to be used on CPU side
     * */
    virtual std::unique_ptr<CellPairIterator> iteratorC08() = 0;
protected:
    /**
     * Updates p_COM buffer to current values by only accessing SOA data
     * */
    void updateCOM();

    /// total number of indices
    uint64_t p_molecule_count;
    /// container for all AoS indices
    std::vector<Molecule> p_molecules;
    /// soa for all indices
    SOA p_soa;
    /// center of mass positions for all indices (has size of all sites)
    KW::vec_t<math::d3> p_com;
    /// Pair List for Kokkos force computation
    PairList p_pair_list;
};


#endif //KOMD_MOLECULECONTAINER_H
