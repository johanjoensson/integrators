#include <gtest/gtest.h>
#include <verlet.h>
#include <fstream>
#include <vector>

TEST(VerletIntegrator, ConstantFunction0)
{
    auto F =  [](const double& F){ return 0;};
    Verlet_integrator<> integrator(0, 0, 1e-3);
    std::vector<double> res;
    for(auto i = 0; i < 1000; i++){
        res.push_back(integrator.value());
    }
    for(auto val : res){
        ASSERT_DOUBLE_EQ(val, 0);
    }
}

TEST(VerletIntegrator, ConstantFunction5)
{
    auto F =  [](const double& F){ return 0;};
    Verlet_integrator<> integrator(5, 5, 1e-3);
    std::vector<double> res;
    for(auto i = 0; i < 1000; i++){
        res.push_back(integrator.advance(F));
    }
    for(auto val : res){
        ASSERT_DOUBLE_EQ(val, 5);
    }
}
TEST(VelocityVerletIntegrator, ConstantFunction0)
{
    auto F =  [](const double& F){ return 0;};
    Velocity_Verlet_integrator<> integrator(0, 0, 1e-3);
    std::vector<double> res;
    for(auto i = 0; i < 1000; i++){
        res.push_back(integrator.advance(F));
    }
    for(auto val : res){
        ASSERT_DOUBLE_EQ(val, 0);
    }
}

TEST(VelocityVerletIntegrator, ConstantFunction5)
{
    auto F =  [](const double& F){ return 0;};
    Velocity_Verlet_integrator<> integrator(5, 0, 1e-3);
    std::vector<double> res;
    for(auto i = 0; i < 1000; i++){
        res.push_back(integrator.advance(F));
    }
    for(auto val : res){
        ASSERT_DOUBLE_EQ(val, 5);
    }
}

TEST(LeapfrogIntegrator, ConstantFunction0)
{
    auto F =  [](const double& F){ return 0;};
    Leapfrog_integrator<> integrator(0, 0, 1e-3);
    std::vector<double> res;
    for(auto i = 0; i < 1000; i++){
        res.push_back(integrator.value());
    }
    for(auto val : res){
        ASSERT_DOUBLE_EQ(val, 0);
    }
}

TEST(LeapfrogIntegrator, ConstantFunction5)
{
    auto F =  [](const double& F){ return 0;};
    Leapfrog_integrator<> integrator(5, 5, 1e-3);
    std::vector<double> res;
    for(auto i = 0; i < 1000; i++){
        res.push_back(integrator.advance(F));
    }
    for(auto val : res){
        ASSERT_DOUBLE_EQ(val, 5);
    }
}
