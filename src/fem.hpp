#pragma once

#include <stdlib.h>

#include <cmath>
#include <cstdio>
#include <cstring>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;

/**
 * @class Fem
 * @brief Base class for finite element method implementations
 *
 * Provides core functionality for finite element analysis including:
 * - Degree of freedom management
 * - Boundary condition application
 * - Solution initialization
 * - Common FEM operations
 */
class Fem
{
   public:
    Fem(void);

    virtual ~Fem() {}

    /**
     * @brief Sets the spatial dimension (2D/3D)
     * @param ndm Number of dimensions (2 or 3)
     */
    virtual void setDimension(int ndm) { this->ndm = ndm; }

    /**
     * @brief Sets degrees of freedom per node
     * @param ndf Number of degrees of freedom
     */
    virtual void setDof(int ndf) { this->ndf = ndf; }

    /**
     * @brief Prescribes degrees of freedom with specified values
     * @param dofs Vector of degree of freedom indices to prescribe
     * @param val_type Type of value to prescribe (0=zero, 1=non-zero)
     * @param val Value to prescribe (if val_type is non-zero)
     */
    virtual void prescribeDofs(vector<int> &dofs, int val_type, double val);

    /**
     * @brief Initializes variables
     */
    virtual void initialiseVars(void);

    /**
     * @brief Gets local points for the current element
     * @param e Element index
     * @return int Local points
     */
    virtual int getLocalPoints(int e);

    /**
     * @brief Gets local coordinates for the current element
     * @param e Element index
     */
    virtual void getLocalCoords(int e);

    /**
     * @brief Gets local solution variables for the current element
     * @param e Element index
     */
    virtual void getLocalVars(int e);

    /**
     * @brief Performs Gauss quadrature
     * @param e Element index
     */
    virtual void gauss(int e)
    {
        cout << "ERROR: Fem::gauss : not defined.";
        exit(1);
    }

    /**
     * @brief Evaluates shape functions
     * @param e Element index
     */
    virtual void shape(int e)
    {
        cout << "ERROR: Fem::shape : not defined.";
        exit(1);
    }

    /**
     * @brief Gets element size
     * @param e Element index
     * @return double Element size
     */
    virtual double getElemSize(int e)
    {
        cout << "ERROR: Fem::getElemSize : not defined.";
        exit(1);
    }

    /**
     * @brief Assembles local matrix
     * @param e Element index
     */
    virtual void assembleLocalMatrix(int e)
    {
        cout << "ERROR: Fem::assembleLocalMatrix : not defined.";
        exit(1);
    }

    /**
     * @brief Updates solution variables
     */
    virtual void updateVars(void)
    {
        cout << "ERROR: Fem::updateVars : not defined.";
        exit(1);
    }

    /**
     * @brief Performs adaptive time stepping
     * @param nrIter Number of iterations
     * @param nrOpt Number of optimization steps
     */
    virtual void adaptiveTimeStepping(int nrIter, int nrOpt)
    {
        cout << "ERROR: Fem::adaptiveTimeStepping : not defined.";
        exit(1);
    }

    /**
     * @brief Prints matrix
     * @tparam Type Type of matrix elements
     * @param mat Matrix to print
     */
    template <typename Type>
    void printMat(vector<vector<Type>> &mat);

    /**
     * @brief Prints vector
     * @tparam Type Type of vector elements
     * @param vec Vector to print
     */
    template <typename Type>
    void printVec(vector<Type> &vec);

    // mesh & element
    int ndm, ndf, nen, nel, nnp, ndof, ngp, nst, nshp;
    double dv;

    vector<int> t, id, ifix;
    vector<double> p, shp, pl, ul, ufix, uinc, ucur;

    // solver
    int neq = 0;
    vector<double> u, un, du, dun, ut, dut, Rl, Sl;

    // element
    vector<double> gp, w;

    // time
    double incCo, duCo;
    bool steady, adaptive;
};

inline Fem::Fem(void)
{
    printf(
        "\n========================================\n"
        "  C    -    Crossflow\n"
        "  F    -    Finite\n"
        "  E    -    Element\n"
        "  M    -    Method\n"
        "========================================\n");
}

inline void Fem::prescribeDofs(vector<int> &select, int df, double val)
{
    int ii, jj = 0, kk;

    // sort select
    // sort(select.begin(), select.end());

    // set points
    if (val == 0.)  // fix the points
    {
        for (ii = 0; ii < select.size(); ii++)
        {
            id[select[ii] * ndf + df] = -1;
        }
    }
    else  // prescribe fixed value
    {
        for (ii = 0; ii < select.size(); ii++)
        {
            // jj++;
            kk = select[ii] * ndf + df;

            //  u[select[ii]*ndf+df] = val;
            // dfix.push_back(kk);
            ufix.push_back(val);
            uinc.push_back(0.);  // FIXME: this needs to be fixed
            ucur.push_back(0.);

            id[kk] = -1 - ufix.size();

            // cout << "id " << -1 -jj << "\n" << "u " << val << "\n"; getchar();
        }
    }
}

inline void Fem::initialiseVars(void)
{
    // set id and u
    ndof = nnp * ndf;
    nst = nen * ndf;

    id.resize(ndof);

    u.resize(ndof);
    fill(u.begin(), u.end(), 0);
    ut = u;
    un = u;
    du = u;
    dun = u;
    dut = u;
}

inline int Fem::getLocalPoints(int e)
{
    std::vector<int> el(nen);

    for (int ii = 0; ii < nen; ii++)
        el[ii] = t[e * nen + ii];

    return *el.begin();
}

inline void Fem::getLocalCoords(int e)
{
    int *el = &t[e * nen], ii, jj, kk = 0;

    for (ii = 0; ii < nen; ii++)
    {
        for (jj = 0; jj < ndm; jj++)
            pl[kk++] = p[*el * ndm + jj];

        el++;
    }
}

inline void Fem::getLocalVars(int e)
{
    int *el = &t[e * nen], ii, jj, kk = 0, ll;

    if (steady)
    {
        for (ii = 0; ii < nen; ii++)
        {
            for (jj = 0; jj < ndf; jj++)
            {
                ul[kk++] = u[*el * ndf + jj];
            }
            el++;
        }
    }
    else
    {
        for (ii = 0; ii < nen; ii++)
        {
            for (jj = 0; jj < ndf; jj++)
            {
                ll = *el * ndf + jj;
                ul[kk] = ut[ll];
                ul[kk + nst] = u[ll];
                ul[kk + 2 * nst] = un[ll];
                ul[kk + 3 * nst] = dut[ll];
                kk++;
            }
            el++;
        }
    }
}

template <typename Type>
void Fem::printMat(vector<vector<Type>> &mat)
{
    cout << "\n{ ";

    if (mat.size() > 0)
    {
        for (int i = 0; i < mat.size(); i++)
        {
            if (mat[i].size() > 0)
            {
                cout << i << "_{" << mat[i][0];

                for (int j = 1; j < mat[i].size(); j++)
                {
                    cout << "," << mat[i][j];
                }
            }
            else
                cout << i << "_{ ";

            cout << "}" << " ";
            // cout << "}_" << mat[i].size() << " ";
        }

        cout << "}_" << mat.size() << "\n";
    }
    else
        cout << "}_0\n";
}

template <typename Type>
void Fem::printVec(vector<Type> &vec)
{
    cout << "\n{ ";

    if (vec.size() > 0)
    {
        cout << vec[0];
        for (int i = 1; i < vec.size(); i++)
            cout << "," << vec[i];

        cout << " }_" << vec.size() << "\n";
    }
    else
        cout << "}_0\n";
}