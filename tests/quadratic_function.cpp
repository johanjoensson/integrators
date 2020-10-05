#include <gtest/gtest.h>
#include <verlet.h>
#include <fstream>
#include <vector>
#include <cmath>

TEST(VerletIntegrator, QuadraticFunction5x2px)
{
    auto F =  [](const double& x){ return 10;};
    Verlet_integrator<> integrator(0, 5e-8 + 2e-4, 1e-4);
    std::vector<double> res;
    for(auto i = 0; i < 10000; i++){
        ASSERT_NEAR(integrator.advance(F), 5*std::pow((i + 2)*1e-4, 2) + 2e-4*(i + 2), 5e-8);
    }
}

TEST(VelocityVerletIntegrator, QuadraticFunction5x2p2x)
{
    auto F =  [](const double& x){ return 10;};
    Velocity_Verlet_integrator<> integrator(0, 2, 1e-4, F(0));
    std::vector<double> res;
    for(auto i = 0; i < 100000; i++){
        ASSERT_NEAR(integrator.advance(F), 5*std::pow((i + 1)*1e-4, 2) + 2e-4*(i + 1), 5e-8);
    }
}

TEST(LeapfrogIntegrator, QuadraticFunction5x2px)
{
    auto F =  [](const double& x){ return 10;};
    Leapfrog_integrator<> integrator(0, 5e-8 + 2e-4, 1e-4);
    std::vector<double> res;
    for(auto i = 0; i < 10000; i++){
        ASSERT_NEAR(integrator.advance(F), 5*std::pow((i + 2)*1e-4, 2) + 2e-4*(i + 2), 5e-8);
    }
}
