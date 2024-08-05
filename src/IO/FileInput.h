//
// Created by alex on 8/4/24.
//

#ifndef KOMD_FILEINPUT_H
#define KOMD_FILEINPUT_H


#include <string>

class FileInput {
public:
/**
     * @brief Reading input parameters from given input data file
     *
     * @param[in] input data file
     */
    static bool readFile(const std::string &filename);
};


#endif //KOMD_FILEINPUT_H
