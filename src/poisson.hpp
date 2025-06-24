#pragma once

#include "element.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "time.hpp"
#include "vtk.hpp"

/**
 * @class Poisson problem
 */
class Poisson : public Mesh, public Time, public Element, public Solver, public Vtk
{
   public:
    Poisson() : Mesh(), Time(), Element(), Solver(), Vtk() {
        // Initialize physics properties
        f = 1.0;
        mu = 1.0;
    }

    virtual ~Poisson() = default;

    // Physics methods
    void setProp(const char* name, double val);
    void assembleLocalMatrix(int e) override;

    // Physics properties
    double f;
    double mu;

};

inline void Poisson::setProp(const char *name, double val)
{
    if (strcmp(name, "f") == 0)
    {
        f = val;
        return;
    }
    else if (strcmp(name, "mu") == 0)
    {
        mu = val;
        return;
    }
    else
    {
        printf(
            "ERROR : Poisson::setProp\n"
            "\t-> Property does not exist.\n\n");
        exit(1);
    }
}

inline void Poisson::assembleLocalMatrix(int e)
{
    fill(Rl.begin(), Rl.end(), 0);  // TODO: maybe pointers are faster
    fill(Sl.begin(), Sl.end(), 0);

    getLocalCoords(e);
    getLocalVars(e);

    int i, j, g, ii;
    double coef;

    std::vector<double> gphi(ndm), xgphi(ndm);

    // gauss point loop
    for (g = 0; g < ngp; g++)
    {
        // get shape functions
        shape(g);

        ///////////////////////////////////////////////////////////////////////
        // elemental assembly start
        ///////////////////////////////////////////////////////////////////////

        // gradient of phi
        for (ii = 0; ii < ndm; ii++)
        {
            gphi[ii] = 0.;

            for (i = 0; i < nen; i++)
                gphi[ii] += ul[i] * shp[i * nshp + 1 + ii];
        }

        // laplace term
        coef = mu * dv;

        for (i = 0; i < nen; i++)
        {
            for (ii = 0; ii < ndm; ii++)
            {
                Rl[i] += coef * shp[i * nshp + 1 + ii] * gphi[ii];

                for (j = 0; j < nen; j++)
                {
                    Sl[j * nst + i] += coef * shp[i * nshp + 1 + ii] * shp[j * nshp + 1 + ii];
                }
            }
        }

        // source term
        coef = f * dv;

        for (i = 0; i < nen; i++)
        {
            Rl[i] += coef * shp[i * nshp];
        }

        ///////////////////////////////////////////////////////////////////////
        // elemental assembly end
        ///////////////////////////////////////////////////////////////////////
    }
}