//
// Created by alex on 7/31/24.
//

#ifndef KOMD_INTEGRATOR_H
#define KOMD_INTEGRATOR_H

/**
 * 2 Step Leapfrog integration
 * */
class Integrator {
public:
    /**
     * Creates new Integrator
     * @param delta_t time step size
     * */
    explicit Integrator(double delta_t);

    /**
     * First half step of integration
     * */
    void integrate0();

    /**
     * Second half step of integration
     * */
    void integrate1();
private:
    double m_delta_t;
    bool m_use_soa;
};


#endif //KOMD_INTEGRATOR_H
