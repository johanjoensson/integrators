#include <gtest/gtest.h>
#include <verlet.h>
#include <fstream>
#include <vector>
#include <cmath>

TEST(VerletIntegrator, TrigSin)
{
    auto F =  [](const double& x){ return -x;};
    Verlet_integrator<> integrator(std::sin(0), std::sin(1e-4), 1e-4);
    std::vector<double> res;
    for(auto i = 0; i < 100000; i++){
        ASSERT_NEAR(integrator.advance(F), std::sin(1e-4*(i + 2)), 1e-8);
    }
}

TEST(VerletIntegrator, TrigCos)
{
    auto F =  [](const double& x){ return -x;};
    Verlet_integrator<> integrator(std::cos(0), std::cos(1e-4), 1e-4);
    std::vector<double> res;
    for(auto i = 0; i < 100000; i++){
        ASSERT_NEAR(integrator.advance(F), std::cos(1e-4*(i + 2)), 1e-8);
    }
}

TEST(VelocityVerletIntegrator, TrigSin)
{
    auto F =  [](const double& x){ return -x;};
    Velocity_Verlet_integrator<> integrator(std::sin(0), std::cos(0), 1e-4, F(std::sin(0)));
    std::ofstream os;
    os.open("sin_integrated.dat", std::ios::out);
    for(auto i = 0; i < 100000; i++){
        ASSERT_NEAR(integrator.advance(F), std::sin(1e-4*(i + 1)), 1e-8);
        os << integrator.value() << " " << std::sin(1e-4*(i + 2)) << "\n";
    }
    os << "\n\n";
    os.close();
}

TEST(VelocityVerletIntegrator, TrigCos)
{
    auto F =  [](const double& x){ return -x;};
    Velocity_Verlet_integrator<> integrator(std::cos(0), -std::sin(0), 1e-3, F(std::cos(0)));
    std::ofstream os;
    os.open("cos_integrated.dat", std::ios::out);

    for(auto i = 0; i < 100000; i++){
        ASSERT_NEAR(integrator.advance(F), std::cos(1e-3*(i + 1)), 5e-6) << " i = " << i;
        os << integrator.value() << " " << std::cos(1e-4*(i + 2)) << "\n";
    }
    os << "\n\n";
    os.close();
}

TEST(LeapfrogIntegrator, TrigSin)
{
    auto F =  [](const double& x){ return -x;};
    Verlet_integrator<> integrator(std::sin(0), std::sin(1e-4), 1e-4);
    for(auto i = 0; i < 100000; i++){
        ASSERT_NEAR(integrator.advance(F), std::sin(1e-4*(i + 2)), 1e-8);
    }
}

TEST(LeapfrogIntegrator, TrigCos)
{
    auto F =  [](const double& x){ return -x;};
    Verlet_integrator<> integrator(std::cos(0), std::cos(1e-4), 1e-4);
    for(auto i = 0; i < 100000; i++){
        ASSERT_NEAR(integrator.advance(F), std::cos(1e-4*(i + 2)), 1e-8);
    }
}
