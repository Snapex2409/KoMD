/*
 * VTKMoleculeWriterImplementation.cpp from ls-Mardyn
 *
 * @Date: 25.08.2010
 * @Author: eckhardw
 */


#include "VTKWriter.h"
#include "Logging.h"
#include "Registry.h"
#include "vtk/vtk-punstructured.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void VTKWriter::write(const std::string &filename, int iteration) {
    auto container = Registry::instance->moleculeContainer();
    auto& cells = container->getCells();
    const auto cell_dims = cells.dims();

    uint64_t sites = 0;
    for (uint64_t z = 0; z < cell_dims.z(); z++) {
        for (uint64_t y = 0; y < cell_dims.y(); y++) {
            for (uint64_t x = 0; x < cell_dims.x(); x++) {
                Cell& cell = cells[x, y, z];
                for (Molecule& molecule : cell.molecules()) sites += molecule.getSites().size();
            }
        }
    }

    initializeOutput(sites);

    for (uint64_t z = 0; z < cell_dims.z(); z++) {
        for (uint64_t y = 0; y < cell_dims.y(); y++) {
            for (uint64_t x = 0; x < cell_dims.x(); x++) {
                Cell& cell = cells[x, y, z];
                for (Molecule& molecule : cell.molecules()) plotMolecule(molecule);
            }
        }
    }

    writeFile(filename, iteration);
}

void VTKWriter::initializeOutput(int numParticles) {

    vtkFile = new VTKFile_t("UnstructuredGrid");

    // per point, we add type, position, velocity and force
    PointData pointData;
    DataArray_t mass(type::Float32, "mass", 1);
    DataArray_t velocity(type::Float32, "velocity", 3);
    DataArray_t forces(type::Float32, "force", 3);
    DataArray_t id(type::Int32, "id", 1);
    pointData.DataArray().push_back(mass);
    pointData.DataArray().push_back(velocity);
    pointData.DataArray().push_back(forces);
    pointData.DataArray().push_back(id);

    CellData cellData; // we don't have cell data => leave it empty

    // 3 coordinates
    Points points;
    DataArray_t pointCoordinates(type::Float32, "points", 3);
    points.DataArray().push_back(pointCoordinates);

    Cells cells; // we don't have cells, => leave it empty
    // for some reasons, we have to add a dummy entry for paraview
    DataArray_t cells_data(type::Float32, "types", 0);
    cells.DataArray().push_back(cells_data);

    PieceUnstructuredGrid_t piece(pointData, cellData, points, cells, numParticles, 0);
    UnstructuredGrid_t unstructuredGrid(piece);
    vtkFile->UnstructuredGrid(unstructuredGrid);
}

void VTKWriter::writeFile(const std::string &filename, int iteration) {
    std::stringstream strstr;
    strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".vtu";

    std::ofstream file(strstr.str().c_str());
    VTKFile(file, *vtkFile);
    delete vtkFile;
}

void VTKWriter::plotMolecule(Molecule &m) {
    if (!vtkFile->UnstructuredGrid().present()) {
        Log::io->error() << "No UnstructuredGrid present" << std::endl;
    }

    for (Site& site : m.getSites()) {
        PointData::DataArray_sequence &pointDataSequence = vtkFile->UnstructuredGrid()->Piece().PointData().DataArray();
        PointData::DataArray_iterator dataIterator = pointDataSequence.begin();

        dataIterator->push_back(site.getMass());
        //loggers::general->debug("Appended mass data in: {}", dataIterator->Name());

        dataIterator++;
        dataIterator->push_back(site.v_arr().x());
        dataIterator->push_back(site.v_arr().y());
        dataIterator->push_back(site.v_arr().z());
        //loggers::general->debug("Appended velocity data in: {}", dataIterator->Name());

        dataIterator++;
        dataIterator->push_back(site.f_arr().x());
        dataIterator->push_back(site.f_arr().y());
        dataIterator->push_back(site.f_arr().z());
        //loggers::general->debug("Appended force data in: {}", dataIterator->Name());

        dataIterator++;
        dataIterator->push_back(m.ID());

        Points::DataArray_sequence &pointsSequence =
                vtkFile->UnstructuredGrid()->Piece().Points().DataArray();
        Points::DataArray_iterator pointsIterator = pointsSequence.begin();
        pointsIterator->push_back(site.r_arr().x());
        pointsIterator->push_back(site.r_arr().y());
        pointsIterator->push_back(site.r_arr().z());
    }
}

