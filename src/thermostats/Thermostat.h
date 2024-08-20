//
// Created by alex on 8/10/24.
//

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

/**
 * General thermostat interface
 * */
class Thermostat {
public:
    virtual ~Thermostat() = default;
    /// Gets called from simulation at some point
    virtual void apply() = 0;
};



#endif //THERMOSTAT_H
