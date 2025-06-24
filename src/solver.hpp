#pragma once

#include <Eigen/OrderingMethods>
#include <Eigen/Sparse>

#include "fem.hpp"

using namespace Eigen;

/**
 * @class Solver
 * @brief Handles linear system solving for FEM simulations
 *
 * Provides functionality for:
 * - Assembling global stiffness matrix
 * - Solving linear systems
 * - Managing solver configurations
 *
 * Uses Eigen for sparse matrix operations.
 */
class Solver : public virtual Fem
{
   public:
    Solver(void);

    virtual ~Solver() {}

    /**
     * @brief Prepares the solver system
     *
     * Initializes matrices and vectors needed for solving.
     */
    virtual void prepareSolver(void);

    /**
     * @brief Solves the linear system
     *
     * Computes solution using configured solver method.
     */
    virtual void solve(void);

    /**
     * @brief Assembles global stiffness matrix
     * @param first Flag to reset matrix before assembly
     */
    virtual void assembleGlobalMatrix(bool first);

    /// Global stiffness matrix (sparse)
    SparseMatrix<double> Sg;
    /// Global right-hand side vector
    VectorXd Rg;
    /// Solution vector
    VectorXd sol;

    // solver type
    // SimplicialLDLT<SparseMatrix<double>> eigSolver;
    SparseLU<SparseMatrix<double>, COLAMDOrdering<int>>
        eigSolver;  // ordering method COLAMDOrdering<int>

    // BiCGSTAB<SparseMatrix<double> > eigSolver;
    // ConjugateGradient<SparseMatrix<double>,Upper> eigSolver;

    typedef Triplet<double> T;
    vector<T> tripletList;

    vector<double> resNorm;

    int nrMax, nrOpt;
    double nrTol;
};

inline Solver::Solver(void) : Fem()
{
    nrMax = 10;
    nrTol = 1.e-6;
    nrOpt = 4;
}

inline void Solver::prepareSolver(void)
{
    int ii, n = nnp * ndf;

    // get number of equations to solve
    for (ii = 0; ii < n; ii++)
        if (id[ii] == 0) id[ii] = neq++;

    // resize global matrices
    Sg.resize(neq, neq);
    Rg.resize(neq);
    sol.resize(neq);
    tripletList.reserve(3 * nel * nst * nst);

    // resize local matrices
    Sl.resize(nst * nst);
    Rl.resize(nst);

    // set shape functions size
    nshp = ndm + 1;
    gauss(ngp);
    shp.resize(nen * nshp);

    // set local solution variable and element point vectors
    if (steady)
        ul.resize(nen * ndf);
    else
        ul.resize(nen * ndf * 4);

    pl.resize(nen * ndm);
}

inline void Solver::solve(void)
{
    clock_t solve_time, assem_time;
    double time_tot;
    bool checkRes = true;
    int i, j, k, nrIt = 0;

    printf("\nStarting iteration:\n");

    resNorm.clear();

    for (i = 0; i < ufix.size(); i++)
    {
        uinc[i] = -ucur[i];
        ucur[i] = ufix[i];
        uinc[i] += ucur[i];
    }

    // Newton Raphson
    while (checkRes)
    {
        updateVars();

        sol.setZero();

        assembleGlobalMatrix((nrIt == 0));

        resNorm.push_back(Rg.norm());

        printf("%2d. norm(residual) = %9.3e\n", nrIt, resNorm.back());

        if (resNorm.back() < nrTol || nrIt == nrMax) break;

        solve_time = clock();

        eigSolver.compute(Sg);
        if (eigSolver.info() != Success)
        {
            cout << "\nxSolver: Failed decomposition!\n";
            return;
        }

        sol = eigSolver.solve(Rg);

        if (eigSolver.info() != Success)
        {
            cout << "\nxSolver: Failed matrix solving!\n";
            return;
        }

        time_tot = (clock() - solve_time) / (double)CLOCKS_PER_SEC;
        printf("solver time: %.6f\n", time_tot);

        // update degrees of freedom
        for (i = 0; i < ndof; i++)
        {
            // free
            k = id[i];
            if (k > -1) u[i] -= incCo * sol[k];

            // prescribed
            k = -2 - id[i];
            if (k > -1) u[i] = ucur[k];
        }

        nrIt++;
    }

    adaptiveTimeStepping(nrIt, nrOpt);
}

inline void Solver::assembleGlobalMatrix(bool first)
{
    int i, j, k, l, m, e, *el, loci;
    std::vector<int> loc(nst);

    double coef = 1. / incCo;

    Sg.setZero();
    Rg.setZero();

    // elemental loop
    for (e = 0; e < nel; e++)
    {
        // local assembly
        assembleLocalMatrix(e);

        el = &t[e * nen];

        // global assembly
        for (i = 0; i < nen; i++)
            for (j = 0; j < ndf; j++)
            {
                k = el[i] * ndf + j;
                l = id[k];
                m = i * ndf + j;
                loc[m] = l;
            }

        for (i = 0; i < nst; i++)
        {
            loci = loc[i];

            if (loci > -1)
            {
                for (j = 0; j < nst; j++)
                {
                    if (loc[j] > -1)
                        tripletList.push_back(T(loci, loc[j], Sl[j * nst + i]));

                    else if (loc[j] < -1)
                    {
                        k = -2 - loc[j];

                        if (first)
                        {
                            Rg[loci] += Sl[j * nst + i] * uinc[k] * coef;
                        }
                    }
                }
                Rg[loci] += Rl[i];
            }
        }
    }

    Sg.setFromTriplets(tripletList.begin(), tripletList.end());
    Sg.makeCompressed();
    tripletList.clear();
}