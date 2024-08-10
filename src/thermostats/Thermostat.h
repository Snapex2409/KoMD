//
// Created by alex on 8/10/24.
//

#ifndef THERMOSTAT_H
#define THERMOSTAT_H

class Thermostat {
public:
    virtual ~Thermostat() = default;
    virtual void apply() = 0;
};



#endif //THERMOSTAT_H
