//
// Created by alex on 8/8/24.
//

#ifndef KOMD_CHECKPOINTIO_H
#define KOMD_CHECKPOINTIO_H


#include <cstdint>

class CheckpointIO {
public:
    static void writeCheckpoint(uint64_t simstep);
    static void loadCheckpoint();
};


#endif //KOMD_CHECKPOINTIO_H
