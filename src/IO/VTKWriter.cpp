/*
 * VTKMoleculeWriterImplementation.cpp from ls-Mardyn
 *
 * @Date: 25.08.2010
 * @Author: eckhardw
 */


#include "VTKWriter.h"
#if defined(USE_VTK)
#include "Logging.h"
#include "Registry.h"
#include "vtk/vtk-punstructured.h"
#endif

#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void VTKWriter::write(const std::string &filename, int iteration) {
#if defined(USE_VTK)
    auto container = Registry::instance->moleculeContainer();

    uint64_t sites = 0;
    for (auto it = container->iterator(MoleculeContainer::MOLECULE); it.isValid(); ++it) {
        Molecule& molecule = it.molecule();
        sites += molecule.getSites().size();
    }

    initializeOutput(sites);

    for (auto it = container->iterator(MoleculeContainer::MOLECULE); it.isValid(); ++it) {
        Molecule& molecule = it.molecule();
        plotMolecule(molecule);
    }

    writeFile(filename, iteration);
#endif
}

void VTKWriter::initializeOutput(int numParticles) {
#if defined(USE_VTK)
    vtkFile = new VTKFile_t("UnstructuredGrid");

    // per point, we add type, position, velocity and force
    PointData pointData;
    DataArray_t mass(type::Float32, "mass", 1);
    DataArray_t velocity(type::Float32, "velocity", 3);
    DataArray_t forces(type::Float32, "force", 3);
    DataArray_t id(type::Int32, "id", 1);
    DataArray_t cid(type::Int32, "cid", 1);
    pointData.DataArray().push_back(mass);
    pointData.DataArray().push_back(velocity);
    pointData.DataArray().push_back(forces);
    pointData.DataArray().push_back(id);
    pointData.DataArray().push_back(cid);

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
#endif
}

void VTKWriter::writeFile(const std::string &filename, int iteration) {
#if defined(USE_VTK)
    std::stringstream strstr;
    strstr << filename << "_" << std::setfill('0') << std::setw(4) << iteration << ".vtu";

    std::ofstream file(strstr.str().c_str());
    VTKFile(file, *vtkFile);
    delete vtkFile;
#endif
}

void VTKWriter::plotMolecule(Molecule &m) {
#if defined(USE_VTK)
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

        dataIterator++;
        dataIterator->push_back(m.CID());

        Points::DataArray_sequence &pointsSequence =
                vtkFile->UnstructuredGrid()->Piece().Points().DataArray();
        Points::DataArray_iterator pointsIterator = pointsSequence.begin();
        pointsIterator->push_back(site.r_arr().x());
        pointsIterator->push_back(site.r_arr().y());
        pointsIterator->push_back(site.r_arr().z());
    }
#endif
}

