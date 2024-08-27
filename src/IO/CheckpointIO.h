//
// Created by alex on 8/8/24.
//

#ifndef KOMD_CHECKPOINTIO_H
#define KOMD_CHECKPOINTIO_H


#include <cstdint>

class CheckpointIO {
public:
    /// writes a single checkpoint file
    static void writeCheckpoint(uint64_t simstep);
    /// loads the checkpoint defined by configuration and stores all indices in the molecule container
    static void loadCheckpoint();
};


#endif //KOMD_CHECKPOINTIO_H
