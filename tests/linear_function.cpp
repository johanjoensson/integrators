#include <gtest/gtest.h>
#include <verlet.h>
#include <fstream>
#include <vector>

TEST(VerletIntegrator, LinearFunction5x)
{
    auto F =  [](const double&){ return 0;};
    Verlet_integrator<> integrator(0, 5e-3, 1e-3);
    std::vector<double> res;
    for(int i = 0; i < 10000; i++){
        ASSERT_NEAR(integrator.advance(F), 5e-3*(i + 2), 1e-10);
    }
}

TEST(VerletIntegrator, LinearFunction5xm2)
{
    auto F =  [](const double&){ return 0;};
    Verlet_integrator<> integrator(-2, -2 + 5e-3, 1e-3);
    std::vector<double> res;
    for(auto i = 0; i < 10000; i++){
        ASSERT_NEAR(integrator.advance(F), 5e-3*(i + 2) - 2, 1e-10);
    }
}

TEST(VelocityVerletIntegrator, LinearFunction5x)
{
    auto F =  [](const double&){ return 0;};
    Velocity_Verlet_integrator<> integrator(0, 5, 1e-5);
    std::vector<double> res;
    for(int i = 0; i < 10000; i++){
        ASSERT_NEAR(integrator.advance(F), 5e-5*(i + 1), 1e-10);
    }
}

TEST(VelocityVerletIntegrator, LinearFunction5xm2)
{
    auto F =  [](const double&){ return 0;};
    Velocity_Verlet_integrator<> integrator(-2, 5, 1e-5);
    std::vector<double> res;
    for(int i = 0; i < 10000; i++){
        ASSERT_NEAR(integrator.advance(F), 5e-5*(i + 1) - 2, 1e-10);
    }
}

TEST(LeapfrogIntegrator, LinearFunction5x)
{
    auto F =  [](const double&){ return 0;};
    Leapfrog_integrator<> integrator(0, 5e-3, 1e-3);
    std::vector<double> res;
    for(int i = 0; i < 10000; i++){
        ASSERT_NEAR(integrator.advance(F), 5e-3*(i + 2), 1e-10);
    }
}

TEST(LeapfrogIntegrator, LinearFunction5xm2)
{
    auto F =  [](const double&){ return 0;};
    Leapfrog_integrator<> integrator(-2, -2 + 5e-3, 1e-3);
    std::vector<double> res;
    for(auto i = 0; i < 10000; i++){
        ASSERT_NEAR(integrator.advance(F), 5e-3*(i + 2) - 2, 1e-10);
    }
}
