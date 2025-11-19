#include "SPHsolver.h"
#include <cmath>
#include <iostream>


int SPHsolver::window_width = 800;
int SPHsolver::window_height = 600;
double SPHsolver::view_width = 1.5 * SPHsolver::window_width;
double SPHsolver::view_height = 1.5 * SPHsolver::window_height;

SPHsolver::SPHsolver()
{
    h = 16.f;
    rest_dens = 300.f;
    gas_const = 2000.f;
    visc = 200.f;
    mass = 2.5f;
    dt = 0.0007f;
    hsq = h * h;
    bound_damp = -0.5f;
    max_particles = 2500;
    dam_particles = 500;
    block_particles = 250;
    poly6 = 4.f / (M_PI * pow(h,8.f));
    spiky_grad = -10.f / (M_PI * pow(h,5.f));
    visc_lap = 40.f / (M_PI * pow(h,5.f));
    eps = h;
}

void SPHsolver::InitSPH()
{
    std::cout << "Initializing Dam Break with " << dam_particles << " particles." << std::endl;
    for (float y = eps; y < view_height - eps * 2.f; y += h)
    {
        for (float x = view_width / 4; x <= view_width / 2; x += h)
        {
            if (particles.size() < dam_particles)
            {
                float jitter = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
                particles.push_back(Particle(x + jitter, y));
            }
            else
            {
                return;
            }
        }
    }
}

void SPHsolver::step()
{
    ComputeDensityPressure();
    ComputeForces();
    Integrate();
}

void SPHsolver::ComputeDensityPressure()
{
    for (auto& pi : particles)
    {
        pi.rho = 0.f;
        for (auto& pj : particles)
        {
            Eigen::Vector2d rij = pj.x - pi.x;
            float r2 = rij.squaredNorm();

            if (r2 < hsq)
            {
                pi.rho+= mass * poly6 * pow(hsq - r2, 3.f);
            }
        }
        pi.p = gas_const * (pi.rho - rest_dens);
    }
}

void SPHsolver::ComputeForces()
{
    for (auto& pi : particles)
    {
        Eigen::Vector2d fpress(0.f, 0.f);
        Eigen::Vector2d fvisc(0.f, 0.f);
        for (auto& pj : particles)
        {
            if (&pi == &pj) continue;
            Eigen::Vector2d rij = pj.x - pi.x;
            float r = rij.norm();

            if (r < h)
            {
                fpress += -rij.normalized() * mass * (pi.p + pj.p) / (2.f * pj.rho) * spiky_grad * pow(h - r, 3.f);
                fvisc += visc * mass * (pj.v - pi.v) / pj.rho * visc_lap * (h - r);
            }
        }
        Eigen::Vector2d G(0.f, -10.f);
        Eigen::Vector2d fgrav = mass * G/pi.rho;
        pi.f = fpress + fvisc + fgrav;
    }
}

void SPHsolver::Integrate()
{
    for (auto& p : particles)
    {
        p.v += dt * p.f / p.rho;
        p.x += dt * p.v;

        // boundary conditions
        if (p.x(0) - eps < 0.f) // left wall
        {
            p.v(0) *= bound_damp;
            p.x(0) = eps;
        }
        if (p.x(0) + eps > view_width) // right wall
        {
            p.v(0) *= bound_damp;
            p.x(0) = view_width - eps;
        }
        if (p.x(1) - eps < 0.f) // bottom wall
        {
            p.v(1) *= bound_damp;
            p.x(1) = eps;
        }
        if (p.x(1) + eps > view_height) // top wall
        {
            p.v(1) *= bound_damp;
            p.x(1) = view_height - eps;
        }
    }
}
