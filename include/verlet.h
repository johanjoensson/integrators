#ifndef VERLET_INTEGRATOR_H
#define VERLET_INTEGRATOR_H

template<class Val = double, class Timestep = double>
class Verlet_integrator
{
protected:
    Timestep h_m;
    Val ancient_val_m;
    Val old_val_m;
    Val newest_val_m;

public:
    Verlet_integrator() = default;
    Verlet_integrator(const Verlet_integrator&) = default;
    Verlet_integrator(Verlet_integrator&&) = default;

    Verlet_integrator(const Val& x0, const Val& x1, const Timestep& h)
     : h_m(h), ancient_val_m(), old_val_m(x0), newest_val_m(x1)
    {}
    virtual ~Verlet_integrator() = default;

    Verlet_integrator& operator=(const Verlet_integrator&) = default;
    Verlet_integrator& operator=(Verlet_integrator&&) = default;

    virtual Val advance(const std::function<Val(const Val&)>& F)
    {
        ancient_val_m = old_val_m;
        old_val_m = newest_val_m;
        newest_val_m = 2*old_val_m - ancient_val_m + h_m*h_m*F(old_val_m);
        return newest_val_m;
    }

    Val value() const
    {
        return newest_val_m;
    }

    virtual Val velocity() const
    {
        return (newest_val_m - ancient_val_m)/(2*h_m);
    }
};

template<class Val = double, class Timestep = double>
class Velocity_Verlet_integrator : public Verlet_integrator<Val, Timestep>
{
private:
    Val velocity_m;
    Val force_m;
public:
    Velocity_Verlet_integrator() = default;
    Velocity_Verlet_integrator(const Velocity_Verlet_integrator&) = default;
    Velocity_Verlet_integrator(Velocity_Verlet_integrator&&) = default;

    Velocity_Verlet_integrator(const Val& x0, const Val& v0, const Timestep& h, const Val& f0 = {})
     : Verlet_integrator<Val, Timestep>({}, x0, h), velocity_m(v0), force_m(f0)
    {}

    Velocity_Verlet_integrator& operator=(const Velocity_Verlet_integrator&) = default;
    Velocity_Verlet_integrator& operator=(Velocity_Verlet_integrator&&) = default;

    Val velocity() const
    {
        return velocity_m;
    }

    Val advance(std::function<Val(const Val&)> F)
    {
        this->ancient_val_m = this->old_val_m;
        this->old_val_m = this->newest_val_m;

        this->newest_val_m += this->velocity_m*this->h_m + this->force_m*this->h_m*this->h_m/2;
        Val new_force = F(this->newest_val_m);

        this->velocity_m += this->h_m*(this->force_m + new_force)/2;
        this->force_m = new_force;
        return this->newest_val_m;
    }
};

template<class Val = double, class Timestep = double>
class Leapfrog_integrator : public Verlet_integrator<Val, Timestep>
{
private:
    Val velocity_m;
public:
    Leapfrog_integrator() = default;
    Leapfrog_integrator(const Leapfrog_integrator&) = default;
    Leapfrog_integrator(Leapfrog_integrator&&) = default;

    Leapfrog_integrator(const Val& x0, const Val& x1, const Timestep& h)
     : Verlet_integrator<Val, Timestep>(x0, x1, h), velocity_m((x1 - x0)/h)
    {}

    Leapfrog_integrator& operator=(const Leapfrog_integrator&) = default;
    Leapfrog_integrator& operator=(Leapfrog_integrator&&) = default;

    Val advance(std::function<Val(const Val&)> F)
    {
        this->ancient_val_m = this->old_val_m;
        this->old_val_m = this->newest_val_m;
        this->velocity_m = this->velocity_m + this->h_m*F(this->newest_val_m);
        this->newest_val_m = this->old_val_m + this->h_m*this->velocity_m;
        return this->newest_val_m;
    }
};

#endif // VERLET_INTEGRATOR_H
