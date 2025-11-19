#pragma once
#include <vector>
#include "Particle.h"

class SPHsolver
{
public:
    SPHsolver();
    void InitSPH();
    void step();

    std::vector<Particle>& GetParticles()  { return particles; }
    std::vector<Particle> particles;

    float h; // kernel radius
    float rest_dens; // rest density
    float gas_const; // gas constant
    float visc; // viscosity constant
    float mass; // particle mass
    float dt; // time step
    float hsq; // radius squared for optimization

    // smoothing kernels (这些变量需要在类外初始化)
    float poly6;
    float spiky_grad;
    float visc_lap;

    // simulation parameters 
    float eps;
    float bound_damp;
    static int window_width;
    static int window_height;
    static double view_width;
    static double view_height;

    // interaction parameters 
    int max_particles;
    int dam_particles;
    int block_particles;

    void ComputeDensityPressure();
    void ComputeForces();
    void Integrate();
};
