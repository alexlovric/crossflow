#pragma once

#include "fem.hpp"
#include <stdexcept>
#include <sstream>

/**
 * @class Time
 * @brief Handles time integration and adaptive time stepping for FEM
 * simulations
 *
 * Provides functionality for:
 * - Time step management
 * - Adaptive time stepping
 * - Time integration schemes
 * - Simulation time tracking
 */
class Time : public virtual Fem
{
   public:
    Time(void);

    virtual ~Time() {}

    /// Current simulation time
    double time;
    /// Current time step size
    double dt;
    /// Minimum allowed time step
    double dtMin;
    /// Maximum allowed time step
    double dtMax;
    /// Final simulation time
    double tFin;
    /// Spectral radius for generalized-alpha method
    double rhoInf;
    /// Time for next output
    double tpr;
    /// Theta parameter for time integration
    double theta;

    /// Current iteration number
    int it;
    /// Maximum number of iterations
    int itMax;

    /**
     * @brief Sets a time parameter
     * @param name Parameter name ("itMax", "dt", etc.)
     * @param val Parameter value
     */
    virtual void setTime(const char *name, double val);

    /**
     * @brief Updates solution variables after time step
     */
    virtual void updateVars(void);

    /**
     * @brief Performs adaptive time stepping
     * @param nrIter Number of iterations in current step
     * @param nrOpt Optimal number of iterations
     */
    virtual void adaptiveTimeStepping(int nrIter, int nrOpt);

    /**
     * @brief Recovers from solver crash by reducing time step
     */
    virtual void crashRecall(void);
};

inline Time::Time(void) : Fem()
{
    time = 0;
    tpr = 0;
    it = 1;
    incCo = 1;
    duCo = 1;
    rhoInf = 0.5;
    theta = 0.75;
}

inline void Time::setTime(const char *name, double val)
{
    if (strcmp(name, "itMax") == 0)
    {
        itMax = (int)val;
        return;
    }
    else if (strcmp(name, "dtMax") == 0)
    {
        dtMax = val;
        return;
    }
    else if (strcmp(name, "dtMin") == 0)
    {
        dtMin = val;
        return;
    }
    else if (strcmp(name, "tFin") == 0)
    {
        tFin = val;
        return;
    }
    else if (strcmp(name, "dt") == 0)
    {
        dt = val;
        return;
    }
    else if (strcmp(name, "steady") == 0)
    {
        steady = (bool)val;
        return;
    }
    else if (strcmp(name, "adaptive") == 0)
    {
        adaptive = (bool)val;
        return;
    }
    else if (strcmp(name, "rhoInf") == 0)
    {
        rhoInf = val;
        return;
    }
    else if (strcmp(name, "theta") == 0)
    {
        theta = val;
        return;
    }
    else
    {
        std::ostringstream oss;
        oss << "Time::setTime: Unknown parameter '" << name << '\'';
        throw std::runtime_error(oss.str());
    }
}

inline void Time::updateVars(void)
{
    int i;

    double coef, coef2;

    if (steady)
    {
        incCo = 1.;
        duCo = 0.;
    }
    else
    // {
    // backward Euler
    //     incCo = 1.;
    //     duCo = 1. / dt;

    //     for (i=0; i<ndof; i++)
    //         du[i] = (u[i] - un[i]) * duCo;

    //      ut = u;
    //     dut = du;
    // }
    {
        // gen alpha

        double af, am, gam, maf, mam;

        af = 1. / (1. + rhoInf);

        am = .5 * (3. - rhoInf) * af;

        gam = .5 + am - af;

        incCo = 1. / af;

        coef = 1. / (gam * dt);

        coef2 = -(1. - gam) / gam;

        duCo = am * coef / af;

        maf = 1. - af;
        mam = 1. - am;

        for (i = 0; i < ndof; i++)
        {
            du[i] = (u[i] - un[i]) * coef + dun[i] * coef2;
            ut[i] = af * u[i] + maf * un[i];
            dut[i] = am * du[i] + mam * dun[i];
        }
    }
}

inline void Time::adaptiveTimeStepping(int nrIter, int nrOpt)
{
    // check if adaptive
    if (!adaptive) return;

    double fact1 = nrIter - nrOpt, fact2 = time - tpr;

    // reduce dt
    if (fact1 > 0)
    {
        fact2 *= pow(theta, fact1);
        dt = max(dtMin, fact2);
    }
    // increase dt
    else if (fact1 < 0)
    {
        fact2 *= pow(1. / theta, -fact1);
        dt = min(dtMax, fact2);
    }
    // remain dt
    else
        return;

    cout << "    #############################################\n"
         << "    Time step size dt = " << dt << "\n\n";
}

inline void Time::crashRecall(void)
{
    // double cutFact = 0.5, lastcoeff = normR.lastCoeff();

    // // reset u
    // u = un;
    // du = dun;
    // mpapTime.curr = mpapTime.prev;

    // int i;
    // for (i=0; i<uDep.size(); i++) uDep[i].curr = uDep[i].prev;

    // // reduce time step size
    // mpapTime.dt = mpapTime.dt * cutFact;

    // mpapTime.hist.n -= 1;
    // // mpapTime.hist.popBack(); //FIXME: check

    // // conditions for exit
    // if (mpapTime.dt < mpapTime.dtMin)
    // 	mos.error("FemAlex::crashRecall: dt < dt_min ... exiting ...");

    // if (std::isnan(lastcoeff))
    // 	mos.error("FemAlex::crashRecall: residual is nan ... exiting ...");

    // // output message
    // cout << "    #############################################\n"
    //      << "    Program crash recovery ...\n"
    //      << "      Recalling previous u and t, and cutting dt ...\n"
    //      << "        New dt = " << mpapTime.dt << "\n\n";
}