#pragma once

#include "element.hpp"
#include "mesh.hpp"
#include "solver.hpp"
#include "time.hpp"
#include "vtk.hpp"

/**
 * @class Navier-Stokes Incompressible Flow
 */
class NavierStokesIncomp : public Mesh, public Time, public Element, public Solver, public Vtk
{
   public:
    NavierStokesIncomp() : Mesh(), Time(), Element(), Solver(), Vtk() {
        // Initialize physics properties
        rho = 1.0;
        mu = 0.01;
        linear = false;
    }

    virtual ~NavierStokesIncomp() = default;

    // Physics methods
    void setProp(const char* name, double val);
    void assembleLocalMatrix(int e) override;

    // Physics properties
    double rho;
    double mu;
    bool linear;

};

inline void NavierStokesIncomp::setProp(const char *name, double val)
{
    if (strcmp(name, "rho") == 0)
    {
        rho = val;
        return;
    }
    else if (strcmp(name, "mu") == 0)
    {
        mu = val;
        return;
    }
    else if (strcmp(name, "linear") == 0)
    {
        linear = (bool)val;
        return;
    }
    else
    {
        printf(
            "ERROR : NavierStokesIncomp::setProp\n"
            "\t-> Property does not exist.\n\n");
        exit(1);
    }
}

inline void NavierStokesIncomp::assembleLocalMatrix(int e)
{
    fill(Rl.begin(), Rl.end(), 0);  // TODO: move into local assembly
    fill(Sl.begin(), Sl.end(), 0);

    getLocalCoords(e);
    getLocalVars(e);

    // get shape functions 1gp
    gauss(1);
    shape(0);

    int i, j, k, ii, jj, g, nst2 = nst + nst, nst3 = nst2 + nst, nstx = nen * ndm;

    double h = getElemSize(e), pr, divU, coef, coef2, coef3, coef4, coef5,
           tauPSPG = 0.25 * h * h / mu, tauSUPG, nrmu, tauFact = 1.;
    std::vector<double> U(ndm), Un(ndm), dU(ndm), gU(ndm * ndm), gUn(ndm * ndm), gUs(ndm * ndm),
        gpr(ndm), UgU(ndm), r(ndm), dtauPSPG(nstx), dtauSUPG(nstx);
    if (steady) linear = false;

    for (i = 0; i < nstx; i++)
    {
        dtauPSPG[i] = 0;
    }

    if (!linear)
    {
        for (j = 0; j < ndm; j++)
        {
            U[j] = 0.;
            for (i = 0; i < nen; i++)
                U[j] += ul[i * ndf + j] * shp[i * nshp];
        }
    }
    else
    {
        for (j = 0; j < ndm; j++)
        {
            U[j] = 0.;
            for (i = 0; i < nen; i++)
            {
                U[j] += ul[i * ndf + j + nst2] * shp[i * nshp];
                //   mos << "ul " << ul[i*ndf+j+nst2] << "\n";
            }
        }
    }

    coef = 0.;
    for (i = 0; i < ndm; i++)
        coef += U[i] * U[i];
    nrmu = sqrt(coef);

    coef = 2. * rho * nrmu / h;

    if (!steady)

        tauSUPG =
            1. / sqrt(1. / (tauPSPG * tauPSPG) + coef * coef);  // + 1./(mpapTime.dt*mpapTime.dt));

    else

        tauSUPG = 1. / sqrt(1. / (tauPSPG * tauPSPG) + coef * coef);

    if (!linear)
    {
        coef = -tauSUPG * tauSUPG * 4. * rho * rho / (h * h);

        for (i = 0; i < nen; i++)
            for (j = 0; j < ndm; j++)
                dtauSUPG[i * ndm + j] = coef * U[j] * shp[i * nshp];
    }
    else
    {
        for (i = 0; i < nstx; i++)
        {
            dtauSUPG[i] = 0;
        }
    }
    tauPSPG *= tauFact;
    tauSUPG *= tauFact;

    gauss(ngp);

    // gauss point loop
    for (g = 0; g < ngp; g++)
    {
        // get shape functions
        shape(g);

        ///////////////////////////////////////////////////////////////////////
        // elemental assembly start
        ///////////////////////////////////////////////////////////////////////

        // #############################
        //  interpolations
        // #############################
        if (!steady)
        {
            for (j = 0; j < ndm; j++)
            {
                dU[j] = 0.;
                for (i = 0; i < nen; i++)
                    dU[j] += ul[j + i * ndf + nst3] * shp[i * nshp];
            }
        }

        // velocity

        for (j = 0; j < ndm; j++)
        {
            U[j] = 0.;
            for (i = 0; i < nen; i++)
                U[j] += ul[j + i * ndf] * shp[i * nshp];
        }

        // pressure

        pr = 0.;

        if (!steady)
            for (i = 0; i < nen; i++)
                pr += ul[ndm + i * ndf + nst] * shp[i * nshp];

        else
            for (i = 0; i < nen; i++)
                pr += ul[ndm + i * ndf] * shp[i * nshp];

        // velocity gradient

        for (j = 0; j < ndm; j++)
            for (k = 0; k < ndm; k++)
            {
                gU[j + k * ndm] = 0.;
                for (i = 0; i < nen; i++)
                    gU[j + k * ndm] += ul[j + i * ndf] * shp[i * nshp + 1 + k];
            }

        // pressure gradient

        if (!steady)

            for (k = 0; k < ndm; k++)
            {
                gpr[k] = 0.;
                for (i = 0; i < nen; i++)
                    gpr[k] += ul[ndm + i * ndf + nst] * shp[i * nshp + 1 + k];
            }

        else

            for (k = 0; k < ndm; k++)
            {
                gpr[k] = 0.;
                for (i = 0; i < nen; i++)
                    gpr[k] += ul[ndm + i * ndf] * shp[i * nshp + 1 + k];
            }

        // div(U)

        divU = gU[0];

        for (i = 1; i < ndm; i++)
            divU += gU[i * (ndm + 1)];

        // #############################
        //  Galerkin terms
        // #############################

        // convection

        if (!linear)
        {
            for (i = 0; i < ndm; i++)
            {
                UgU[i] = 0.;
                for (j = 0; j < ndm; j++)
                    UgU[i] += gU[j * ndm + i] * (U[j]);
            }

            coef = rho * dv;
            for (i = 0; i < nen; i++)
            {
                coef2 = coef * shp[i * nshp];

                for (ii = 0; ii < ndm; ii++)
                {
                    Rl[i * ndf + ii] += coef2 * UgU[ii];

                    for (j = 0; j < nen; j++)
                        for (jj = 0; jj < ndm; jj++)
                        {
                            Sl[(j * ndf + ii) * nst + i * ndf + ii] +=
                                coef2 * (U[jj]) * shp[j * nshp + 1 + jj];
                            Sl[(j * ndf + jj) * nst + i * ndf + ii] +=
                                coef2 * shp[j * nshp] * gU[jj * ndm + ii];
                        }
                }
            }
        }
        else
        {
            for (j = 0; j < ndm; j++)
            {
                Un[j] = 0.;
                for (i = 0; i < nen; i++)
                    Un[j] += ul[j + i * ndf + nst2] * shp[i * nshp];

                for (k = 0; k < ndm; k++)
                {
                    gUn[j + k * ndm] = 0.;
                    for (i = 0; i < nen; i++)
                        gUn[j + k * ndm] += ul[j + i * ndf + nst2] * shp[i * nshp + 1 + k];
                }
            }

            for (i = 0; i < ndm; i++)
            {
                UgU[i] = 0.;
                for (j = 0; j < ndm; j++)
                    UgU[i] += gU[j * ndm + i] * (Un[j]) + gUn[j * ndm + i] * (U[j] - Un[j]);
            }

            coef = rho * dv;
            for (i = 0; i < nen; i++)
            {
                coef2 = coef * shp[i * nshp];

                for (ii = 0; ii < ndm; ii++)
                {
                    Rl[i * ndf + ii] += coef2 * UgU[ii];

                    for (j = 0; j < nen; j++)
                        for (jj = 0; jj < ndm; jj++)
                        {
                            Sl[(j * ndf + ii) * nst + i * ndf + ii] +=
                                coef2 * (Un[jj]) * shp[j * nshp + 1 + jj];
                            Sl[(j * ndf + jj) * nst + i * ndf + ii] +=
                                coef2 * shp[j * nshp] * gUn[jj * ndm + ii];
                        }
                }
            }
        }

        // diffusion

        for (i = 0; i < ndm; i++)
            for (j = 0; j < ndm; j++)
                gUs[i * ndm + j] = gU[i * ndm + j] + gU[j * ndm + i];

        coef = mu * dv;
        for (i = 0; i < nen; i++)
            for (ii = 0; ii < ndm; ii++)
                for (jj = 0; jj < ndm; jj++)
                {
                    coef2 = coef * shp[i * nshp + 1 + jj];

                    Rl[i * ndf + ii] += coef2 * gUs[jj * ndm + ii];

                    for (j = 0; j < nen; j++)
                    {
                        Sl[(j * ndf + ii) * nst + i * ndf + ii] += coef2 * shp[j * nshp + 1 + jj];
                        Sl[(j * ndf + jj) * nst + i * ndf + ii] += coef2 * shp[j * nshp + 1 + ii];
                    }
                }

        // pressure and continuity

        for (i = 0; i < nen; i++)
        {
            Rl[i * ndf + ndm] += dv * shp[i * nshp] * divU;
            for (ii = 0; ii < ndm; ii++)
            {
                Rl[i * ndf + ii] -= dv * shp[i * nshp + 1 + ii] * pr;

                for (j = 0; j < nen; j++)
                {
                    Sl[(j * ndf + ii) * nst + i * ndf + ndm] +=
                        dv * shp[i * nshp] * shp[j * nshp + 1 + ii];
                    Sl[(j * ndf + ndm) * nst + i * ndf + ii] -=
                        dv * shp[j * nshp] * shp[i * nshp + 1 + ii];
                }
            }
        }

        // stabilisation terms

        for (i = 0; i < ndm; i++)
            r[i] = rho * UgU[i] + gpr[i];

        if (!steady)
            for (i = 0; i < ndm; i++)
                r[i] += rho * dU[i];

        if (!linear)
        {
            coef = dv * tauPSPG;

            for (i = 0; i < nen; i++)
            {
                for (ii = 0; ii < ndm; ii++)
                {
                    coef2 = coef * shp[i * nshp + 1 + ii];
                    coef3 = coef2 * rho;

                    Rl[i * ndf + ndm] += coef2 * r[ii];

                    for (j = 0; j < nen; j++)
                    {
                        Sl[(j * ndf + ndm) * nst + i * ndf + ndm] += coef2 * shp[j * nshp + 1 + ii];

                        for (jj = 0; jj < ndm; jj++)
                        {
                            Sl[(j * ndf + ii) * nst + i * ndf + ndm] +=
                                coef3 * (U[jj]) * shp[j * nshp + 1 + jj];
                            Sl[(j * ndf + jj) * nst + i * ndf + ndm] +=
                                coef3 * shp[j * nshp] * gU[jj * ndm + ii] +
                                coef2 * dtauPSPG[j * ndm + jj] * r[ii];
                        }
                    }
                }
            }
            coef = dv * tauSUPG;
            for (i = 0; i < nen; i++)
            {
                coef2 = 0.;
                for (jj = 0; jj < ndm; jj++)
                    coef2 += shp[i * nshp + 1 + jj] * (U[jj]);
                coef2 *= (coef * rho);
                coef3 = coef2 * rho;
                for (ii = 0; ii < ndm; ii++)
                {
                    coef4 = coef * rho * r[ii];

                    Rl[i * ndf + ii] += coef2 * r[ii];

                    for (j = 0; j < nen; j++)
                    {
                        Sl[(j * ndf + ndm) * nst + i * ndf + ii] += coef2 * shp[j * nshp + 1 + ii];

                        for (jj = 0; jj < ndm; jj++)
                        {
                            Sl[(j * ndf + ii) * nst + i * ndf + ii] +=
                                coef3 * (U[jj]) * shp[j * nshp + 1 + jj];
                            Sl[(j * ndf + jj) * nst + i * ndf + ii] +=
                                coef3 * shp[j * nshp] * gU[jj * ndm + ii] +
                                coef4 * shp[j * nshp] * shp[i * nshp + 1 + jj] +
                                coef2 * dtauSUPG[j * ndm + jj] * r[ii];
                        }
                    }
                }
            }
        }
        else
        {
            coef = dv * tauPSPG;
            for (i = 0; i < nen; i++)
            {
                for (ii = 0; ii < ndm; ii++)
                {
                    coef2 = coef * shp[i * nshp + 1 + ii];
                    coef3 = coef2 * rho;

                    Rl[i * ndf + ndm] += coef2 * r[ii];

                    for (j = 0; j < nen; j++)
                    {
                        Sl[(j * ndf + ndm) * nst + i * ndf + ndm] += coef2 * shp[j * nshp + 1 + ii];

                        for (jj = 0; jj < ndm; jj++)
                        {
                            Sl[(j * ndf + ii) * nst + i * ndf + ndm] +=
                                coef3 * (Un[jj]) * shp[j * nshp + 1 + jj];
                            Sl[(j * ndf + jj) * nst + i * ndf + ndm] +=
                                coef3 * shp[j * nshp] * gUn[jj * ndm + ii];
                        }
                    }
                }
            }

            coef = dv * tauSUPG;
            for (i = 0; i < nen; i++)
            {
                coef2 = 0.;
                for (jj = 0; jj < ndm; jj++)
                    coef2 += shp[i * nshp + 1 + jj] * (Un[jj]);
                coef2 *= (coef * rho);
                coef3 = coef2 * rho;

                for (ii = 0; ii < ndm; ii++)
                {
                    coef4 = coef * rho * r[ii];

                    Rl[i * ndf + ii] += coef2 * r[ii];

                    for (j = 0; j < nen; j++)
                    {
                        Sl[(j * ndf + ndm) * nst + i * ndf + ii] += coef2 * shp[j * nshp + 1 + ii];

                        for (jj = 0; jj < ndm; jj++)
                        {
                            Sl[(j * ndf + ii) * nst + i * ndf + ii] +=
                                coef3 * (Un[jj]) * shp[j * nshp + 1 + jj];
                            Sl[(j * ndf + jj) * nst + i * ndf + ii] +=
                                coef3 * shp[j * nshp] * gUn[jj * ndm + ii];
                        }
                    }
                }
            }
        }

        if (!steady)
        {
            // stabilisation

            coef3 = dv * rho * duCo;
            coef = coef3 * rho * tauSUPG;
            coef5 = coef3 * tauPSPG;
            for (i = 0; i < nen; i++)
            {
                coef2 = 0.;
                if (!linear)
                    for (jj = 0; jj < ndm; jj++)
                        coef2 += (U[jj]) * shp[i * nshp + 1 + jj];
                else
                    for (jj = 0; jj < ndm; jj++)
                        coef2 += (Un[jj]) * shp[i * nshp + 1 + jj];
                coef2 *= coef;
                for (j = 0; j < nen; j++)
                {
                    coef4 = coef2 * shp[j * nshp];

                    for (ii = 0; ii < ndm; ii++)
                    {
                        Sl[(j * ndf + ii) * nst + i * ndf + ii] += coef4;
                        Sl[(j * ndf + ii) * nst + i * ndf + ndm] +=
                            coef5 * shp[j * nshp] * shp[i * nshp + 1 + ii];
                    }
                }
            }

            // inertia

            coef = dv * rho;
            coef2 = coef * duCo;

            for (i = 0; i < nen; i++)
            {
                for (ii = 0; ii < ndm; ii++)
                {
                    Rl[i * ndf + ii] += coef * shp[i * nshp] * dU[ii];
                    for (j = 0; j < nen; j++)
                        Sl[(j * ndf + ii) * nst + (i * ndf + ii)] +=
                            coef2 * shp[i * nshp] * shp[j * nshp];
                }
            }
        }

        ///////////////////////////////////////////////////////////////////////
        // elemental assembly end
        ///////////////////////////////////////////////////////////////////////
    }

    // account for pressure not subject to time integration scheme
    for (i = 0; i < nen; i++)
    {
        for (j = 0; j < nst; j++)
        {
            Sl[(i * ndf + ndm) * nst + j] *= incCo;
        }
    }
}
