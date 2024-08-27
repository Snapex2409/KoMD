//
// Created by alex on 8/6/24.
//

#ifndef KOMD_VTKWRITER_H
#define KOMD_VTKWRITER_H

#include <string>
#include "molecule/Molecule.h"

class VTKFile_t;

/**
 * This class implements the functionality to generate vtk output from
 * particles.
 */
class VTKWriter {
public:
    void write(const std::string &filename, int iteration);

private:
    /**
     * set up internal data structures and prepare to plot a particle.
     */
    void initializeOutput(int numParticles);

    /**
     * plot type, mass, position, velocity and force of a particle.
     *
     * @note: initializeOutput() must have been called before.
     */
    void plotMolecule(Molecule &m);

    /**
     * writes the final output file.
     *
     * @param filename the base name of the file to be written.
     * @param iteration the number of the current iteration,
     *        which is used to generate an unique filename
     */
    void writeFile(const std::string &filename, int iteration);

    VTKFile_t *vtkFile;
};


#endif //KOMD_VTKWRITER_H
