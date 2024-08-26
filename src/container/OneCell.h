//
// Created by alex on 8/20/24.
//

#ifndef KOMD_ONECELL_H
#define KOMD_ONECELL_H

#include "MoleculeContainer.h"
#include "util/Kokkos_Wrapper.h"
#include "SOA.h"

class OneCell : public MoleculeContainer{
public:
    OneCell();
    ~OneCell() override = default;

    void addMolecule(const Molecule &molecule) override;
    void updateContainer() override;
    void getCenterOfMassPositions(KW::vec_t<math::d3> &buffer) override;
    void writeSOA2AOS() override;
    void init() override {};

    class OCIterator : public Iterator {
    public:
        OCIterator(bool only_mol, Cell& cell);
        ~OCIterator() override = default;
        void operator++() override;
        bool isValid() const override;

        math::d3& f() override;
        math::d3& r() override;
        math::d3& v() override;
        double epsilon() override;
        double sigma() override;
        double mass() override;
        uint64_t ID() override;

        Molecule& molecule() override;

    private:
        const bool m_only_molecule;
        uint64_t m_site_idx;
        uint64_t m_molecule_idx;
        uint64_t m_visited_sites_cell;
        Cell& m_cell;
    };

    class OCCellIterator : public CellIterator {
    public:
        OCCellIterator(Cell& cell) : m_cell(cell), m_valid(true) {};
        ~OCCellIterator() override = default;
        void operator++() override { m_valid = false; };
        bool isValid() const override { return m_valid; };
        Cell& cell() override { return m_cell; };
    private:
        Cell& m_cell;
        bool m_valid;
    };

    class OCC08Iterator : public CellPairIterator {
    public:
        OCC08Iterator(Cell& cell) : m_cell(cell) {};
        ~OCC08Iterator() override = default;
        void operator++() override { };
        bool isValid() const override {return false; };
        Cell& cell0() override {return m_cell; };
        Cell& cell1() override {return m_cell; };
        bool colorSwitched() override { return false; };
    private:
        Cell& m_cell;
    };

    std::unique_ptr<Iterator> iterator(const IteratorType type, const IteratorRegion region) override;
    std::unique_ptr<CellIterator> iteratorCell(const IteratorRegion region) override;
    std::unique_ptr<CellPairIterator> iteratorC08() override;

    struct Periodic_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> com;
        KW::vec_t<math::d3> r;
        const math::d3 low;
        const math::d3 high;
        const math::d3 domain_size;
    };

    struct FReset_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> f;
    };

    KW::vec_t<math::d3> getCoM() { return m_com; }
private:
    Cell m_data;
    bool m_first_update;
    /// center of mass positions for all molecules (has size of all sites)
    KW::vec_t<math::d3> m_com;
    void constructSOAs() override;
    void clearForces() override;
    void updateCOM();
};


#endif //KOMD_ONECELL_H
