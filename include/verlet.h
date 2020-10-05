#ifndef VERLET_INTEGRATOR_H
#define VERLET_INTEGRATOR_H

template<class Val, class Timestep>
class Integrator
{
public:
  virtual ~Integrator() {};

  virtual Val position() const = 0;
  virtual Val velocity() const = 0;
  virtual Val advance(const std::function<Val(const Val&)>& F) = 0;
};



template<class Val = double, class Timestep = double>
class Verlet_integrator : public Integrator<Val, Timestep>
{
protected:
    Timestep h_m;
    Val ancient_position_m;
    Val old_position_m;
    Val position_m;

public:
    Verlet_integrator() = default;
    Verlet_integrator(const Verlet_integrator&) = default;
    Verlet_integrator(Verlet_integrator&&) = default;

    Verlet_integrator(const Val& x0, const Val& x1, const Timestep& h)
     : h_m(h), ancient_position_m(), old_position_m(x0), position_m(x1)
    {}
    ~Verlet_integrator() = default;

    Verlet_integrator& operator=(const Verlet_integrator&) = default;
    Verlet_integrator& operator=(Verlet_integrator&&) = default;

    virtual Val advance(const std::function<Val(const Val&)>& F) override
    {
        ancient_position_m = old_position_m;
        old_position_m = position_m;
        position_m = 2*position_m - ancient_position_m + h_m*h_m*F(old_position_m);
        return position_m;
    }

    Val position() const override
    {
        return position_m;
    }

    virtual Val velocity() const override
    {
        return (position_m - ancient_position_m)/(2*h_m);
    }
};

template<class Val = double, class Timestep = double>
class Velocity_Verlet_integrator : public Integrator<Val, Timestep>
{
private:
    Val position_m;
    Val velocity_m;
    Val force_m;
    Timestep h_m;
public:
    Velocity_Verlet_integrator() = default;
    Velocity_Verlet_integrator(const Velocity_Verlet_integrator&) = default;
    Velocity_Verlet_integrator(Velocity_Verlet_integrator&&) = default;
    ~Velocity_Verlet_integrator() = default;

    Velocity_Verlet_integrator(const Val& x0, const Val& v0, const Timestep& h, const Val& f0 = {})
     : position_m(x0), velocity_m(v0), force_m(f0), h_m(h)
    {}

    Velocity_Verlet_integrator& operator=(const Velocity_Verlet_integrator&) = default;
    Velocity_Verlet_integrator& operator=(Velocity_Verlet_integrator&&) = default;

    virtual Val position() const override
    {
        return position_m;
    }

    virtual Val velocity() const override
    {
        return velocity_m;
    }

    virtual Val advance(const std::function<Val(const Val&)>& F) override
    {
        this->position_m += this->velocity_m*this->h_m + this->force_m*this->h_m*this->h_m/2;
        Val new_force = F(this->position_m);

        this->velocity_m += this->h_m*(this->force_m + new_force)/2;
        this->force_m = new_force;
        return this->position_m;
    }
};

template<class Val = double, class Timestep = double>
class Leapfrog_integrator : public Integrator<Val, Timestep>
{
private:
    Val position_m;
    Val velocity_m;
    Timestep h_m;
public:
    Leapfrog_integrator() = default;
    Leapfrog_integrator(const Leapfrog_integrator&) = default;
    Leapfrog_integrator(Leapfrog_integrator&&) = default;
    Leapfrog_integrator(const Val& x0, const Val& x1, const Timestep& h)
     :  position_m(x1), velocity_m((x1 - x0)/h), h_m(h)
    {}
    ~Leapfrog_integrator() = default;

    Leapfrog_integrator& operator=(const Leapfrog_integrator&) = default;
    Leapfrog_integrator& operator=(Leapfrog_integrator&&) = default;

    virtual Val position() const override
    {
      return position_m;
    }

    virtual Val velocity() const override {return velocity_m;}

    virtual Val advance(const std::function<Val(const Val&)>& F) override
    {
        this->velocity_m = this->velocity_m + this->h_m*F(this->position_m);
        this->position_m += this->h_m*this->velocity_m;
        return this->position_m;
    }
};

#endif // VERLET_INTEGRATOR_H
