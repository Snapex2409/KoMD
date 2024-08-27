//
// Created by alex on 8/20/24.
//

#ifndef KOMD_ONECELL_H
#define KOMD_ONECELL_H

#include "MoleculeContainer.h"
#include "util/Kokkos_Wrapper.h"
#include "SOA.h"

class OneCell : public MoleculeContainer {
public:
    OneCell();
    ~OneCell() override = default;

    void updateContainer() override;
    void init() override;

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

    std::unique_ptr<CellIterator> iteratorCell() override;
    std::unique_ptr<CellPairIterator> iteratorC08() override;

    struct FReset_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> f;
    };

    struct Periodic_Kernel {
        KOKKOS_FUNCTION void operator()(int idx) const;
        KW::vec_t<math::d3> com;
        KW::vec_t<math::d3> r;
        const math::d3 low;
        const math::d3 high;
        const math::d3 domain_size;
    };
private:
    Cell m_data;
};


#endif //KOMD_ONECELL_H
