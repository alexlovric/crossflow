#pragma once

#include <sstream>
#include <stdexcept>

#include "fem.hpp"

/**
 * @class Element
 * @brief Handles element operations for finite element analysis
 *
 * Provides functionality for:
 * - Gauss quadrature integration
 * - Shape function calculations
 * - Element size computations
 */
class Element : public virtual Fem
{
   public:
    Element(void) : Fem() {}

    virtual ~Element() {}

    /**
     * @brief Sets the number of Gauss points for integration
     * @param ngpx Number of Gauss points
     */
    virtual void setNgp(int ngpx);

    /**
     * @brief Sets up Gauss integration points and weights for the element
     * @param ngpx Number of Gauss points
     */
    virtual void gauss(int ngpx);

    /**
     * @brief Sets up 1D Gauss integration points and weights
     */
    virtual void gauss1D(void);

    /**
     * @brief Sets up 2D quadrilateral Gauss integration points and weights
     * @param ngpx Number of Gauss points
     */
    virtual void gaussQuad2D(int ngpx);

    /**
     * @brief Sets up 2D triangular Gauss integration points and weights
     * @param ngpx Number of Gauss points
     */
    virtual void gaussTri2D(int ngpx);

    /**
     * @brief Computes shape functions for the element
     * @param g Gauss point index
     */
    virtual void shape(int g);

    /**
     * @brief Computes quadrilateral shape functions
     * @param g Gauss point index
     */
    virtual void shapeQuad2D(int g);

    /**
     * @brief Computes triangular shape functions
     * @param g Gauss point index
     */
    virtual void shapeTri2D(int g);

    /**
     * @brief Gets element size
     * @param e Element index
     * @return double Element size
     */
    virtual double getElemSize(int e);

    /**
     * @brief Gets quadrilateral element size
     * @param e Element index
     * @return double Element size
     */
    virtual double getElemSizeQuad2D(int e);

    /**
     * @brief Gets triangular element size
     * @param e Element index
     * @return double Element size
     */
    virtual double getElemSizeTri2D(int e);

    /**
     * @brief Computes area of a triangle given three points
     * @param x0 Coordinates of first point
     * @param x1 Coordinates of second point
     * @param x2 Coordinates of third point
     * @return double Area of the triangle
     */
    virtual double getAreaTri2D(double *x0, double *x1, double *x2);
};

inline void Element::gauss1D(void)
{
    gp.resize(ngp);
    w.resize(ngp);

    switch (ngp)
    {
        case 1:
            gp[0] = 0.0;
            w[0] = 2.0;
            break;
        case 2:
            gp[0] = -0.577350269189626;
            gp[1] = 0.577350269189626;
            w[0] = 1.0;
            w[1] = 1.0;
            break;
        case 3:
            gp[0] = -0.774596669241483;
            gp[1] = 0.0;
            gp[2] = 0.774596669241483;
            w[0] = 0.555555555555556;
            w[1] = 0.888888888888889;
            w[2] = 0.555555555555556;
            break;
        case 4:
            gp[0] = -0.861136311594053;
            gp[1] = -0.339981043584856;
            gp[2] = 0.339981043584856;
            gp[3] = 0.861136311594053;
            w[0] = 0.347854845137454;
            w[1] = 0.652145154862546;
            w[2] = 0.652145154862546;
            w[3] = 0.347854845137454;
            break;
        default:
            std::ostringstream oss;
            oss << "Element: Invalid number of Gauss points: " << ngp;
            throw std::runtime_error(oss.str());
    }
}

inline void Element::setNgp(int ngpx)
{
    ngp = ngpx;
}

inline void Element::gauss(int ngpx)
{
    if (ndm == 2) switch (nen)
        {
            // triangle
            case 3:
                gaussTri2D(ngpx);

                break;

            // quad
            case 4:
                gaussQuad2D(ngpx);

                break;

            default:
                throw std::runtime_error(
                    "Element::gauss: Unsupported number of element nodes (nen = " +
                    std::to_string(nen) + "). Only 3 (triangle) and 4 (quad) are supported.");
        }

    else
    {
        throw std::runtime_error("Element::gauss: 3D elements not yet implemented.");
    }
}

inline void Element::gaussQuad2D(int ngpx)
{
    gp.resize(ngpx * 2);
    w.resize(ngpx);

    switch (ngpx)
    {
        // 1 integration point
        case 1:
            gp[0] = 0;
            gp[1] = 0;

            w[0] = 4;

            break;

        // 2 integration point (4 in quadrilateral)
        case 4:
            gp[0] = -sqrt(1. / 3.);
            gp[1] = -sqrt(1. / 3.);
            gp[2] = -sqrt(1. / 3.);
            gp[3] = sqrt(1. / 3.);
            gp[4] = sqrt(1. / 3.);
            gp[5] = -sqrt(1. / 3.);
            gp[6] = sqrt(1. / 3.);
            gp[7] = sqrt(1. / 3.);

            w[0] = 1;
            w[1] = 1;
            w[2] = 1;
            w[3] = 1;

            break;

        // 3 integration point (9 in quadrilateral)
        case 9:
            gp[0] = -sqrt(3. / 5.);
            gp[1] = -sqrt(3. / 5.);
            gp[2] = -sqrt(3. / 5.);
            gp[3] = 0.;
            gp[4] = -sqrt(3. / 5.);
            gp[5] = sqrt(3. / 5.);
            gp[6] = 0.;
            gp[7] = -sqrt(3. / 5.);
            gp[8] = 0.;
            gp[9] = 0.;
            gp[10] = 0.;
            gp[11] = sqrt(3. / 5.);
            gp[12] = sqrt(3. / 5.);
            gp[13] = -sqrt(3. / 5.);
            gp[14] = sqrt(3. / 5.);
            gp[15] = 0.;
            gp[16] = sqrt(3. / 5.);
            gp[17] = sqrt(3. / 5.);

            w[0] = 0.308641975308642;
            w[1] = 0.493827160493828;
            w[2] = 0.308641975308642;
            w[3] = 0.493827160493828;
            w[4] = 0.790123456790124;
            w[5] = 0.493827160493828;
            w[6] = 0.308641975308642;
            w[7] = 0.493827160493828;
            w[8] = 0.308641975308642;

            break;

        default:
            throw std::runtime_error("Element::gaussQuad2D: Unsupported number of Gauss points (" +
                                     std::to_string(ngpx) +
                                     "). Only 1, 4, and 9 points are supported.");
    }
}

inline void Element::gaussTri2D(int ngpx)
{
    gp.resize(ngpx * 2);
    w.resize(ngpx);

    switch (ngpx)
    {
        case 1:
            gp[0] = 1. / 3.;
            gp[1] = gp[0];

            w[0] = .5;

            break;

        case 3:
            gp[0] = 1. / 6.;
            gp[1] = gp[0];
            gp[2] = 2. / 3.;
            gp[3] = gp[0];
            gp[4] = gp[0];
            gp[5] = gp[2];

            w[0] = 1. / 3.;
            w[1] = w[0];
            w[2] = w[0];

            break;

        case 6:
            gp[0] = 0.816847572980459;
            gp[1] = 0.091576213509771;
            gp[2] = gp[1];
            gp[3] = gp[0];
            gp[4] = gp[1];
            gp[5] = gp[1];
            gp[6] = 0.108103018168070;
            gp[7] = 0.445948490915965;
            gp[8] = gp[7];
            gp[9] = gp[6];
            gp[10] = gp[7];
            gp[11] = gp[7];

            w[0] = 0.109951743655322;
            w[1] = w[0];
            w[2] = w[0];
            w[3] = 0.223381589678011;
            w[4] = w[3];
            w[5] = w[3];

            break;

        default:
            throw std::runtime_error("Element::gaussTri2D: Unsupported number of Gauss points (" +
                                     std::to_string(ngpx) +
                                     "). Only 1, 3, 4, 6, 7, and 12 points are supported.");
    }
}

inline void Element::shape(int g)
{
    if (ndm == 2) switch (nen)
        {
            // triangle
            case 3:
                shapeTri2D(g);

                break;

            // quad
            case 4:
                shapeQuad2D(g);

                break;

            default:
                throw std::runtime_error(
                    "Element::shape: Unsupported number of element nodes (nen = " +
                    std::to_string(nen) + "). Only 3 (triangle) and 4 (quad) are supported.");
        }

    else
    {
        throw std::runtime_error("Element::shape: 3D elements not yet implemented.");
    }
}

inline void Element::shapeQuad2D(int g)
{
    int ii, jj, kk;
    double xi = gp[g * 2], eta = gp[g * 2 + 1], dJ, idJ;

    // shape functions
    shp[0] = 0.25 * (1 - xi) * (1 - eta);
    shp[3] = 0.25 * (1 + xi) * (1 - eta);
    shp[6] = 0.25 * (1 + xi) * (1 + eta);
    shp[9] = 0.25 * (1 - xi) * (1 + eta);

    // shape function derivatives
    double J[] = {0., 0., 0., 0.}, iJ[4], dshp[8];

    // shape function derivative with xi and eta
    dshp[0] = -0.25 * (1 - eta);
    dshp[1] = -0.25 * (1 - xi);
    dshp[2] = 0.25 * (1 - eta);
    dshp[3] = -0.25 * (1 + xi);
    dshp[4] = 0.25 * (1 + eta);
    dshp[5] = 0.25 * (1 + xi);
    dshp[6] = -0.25 * (1 + eta);
    dshp[7] = 0.25 * (1 - xi);

    // obtaining Jacobian matrix
    for (ii = 0; ii < 2; ii++)
        for (jj = 0; jj < 2; jj++)
            for (kk = 0; kk < 4; kk++)
                J[ii + jj + jj] += dshp[jj + kk + kk] * pl[ii + kk + kk];

    // determinant of the Jacobian
    dJ = J[0] * J[3] - J[1] * J[2];

    // inverse of the Jacobian
    idJ = 1. / dJ;

    iJ[0] = idJ * J[3];
    iJ[1] = -idJ * J[2];
    iJ[2] = -idJ * J[1];
    iJ[3] = idJ * J[0];

    // adding derivatives to shape function array
    for (ii = 0; ii < 2; ii++)
        for (jj = 0; jj < 4; jj++)
            shp[jj * nshp + ii + 1] =
                iJ[ii + ii] * dshp[jj + jj] + iJ[ii + ii + 1] * dshp[jj + jj + 1];

    dv = dJ * w[g];
}

inline void Element::shapeTri2D(int g)
{
    double xi = gp[g * 2], eta = gp[g * 2 + 1], dJ, a, x1 = pl[0], x2 = pl[2], x3 = pl[4],
           y1 = pl[1], y2 = pl[3], y3 = pl[5];

    dJ = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    a = 0.5 / dJ;

    // shape functions
    shp[0] = 1. - xi - eta;
    shp[3] = xi;
    shp[6] = eta;

    // shape function derivatives
    shp[1] = (y2 - y3) * a;
    shp[2] = (x3 - x2) * a;
    shp[4] = (y3 - y1) * a;
    shp[5] = (x1 - x3) * a;
    shp[7] = (y1 - y2) * a;
    shp[8] = (x2 - x1) * a;

    dv = dJ * w[g];
}

inline double Element::getElemSize(int e)
{
    if (ndm == 2) switch (nen)
        {
            // triangle
            case 3:
                return getElemSizeTri2D(e);  // TODO:

            // quad
            case 4:
                return getElemSizeQuad2D(e);

            default:
                throw std::runtime_error(
                    "Element::getElemSize: Unsupported number of element nodes (nen = " +
                    std::to_string(nen) + "). Only 3 (triangle) and 4 (quad) are supported.");
        }

    else
    {
        throw std::runtime_error("Element::getElemSize: 3D elements not yet implemented.");
    }
}

inline double Element::getElemSizeQuad2D(int e)
{
    // FIXME: need to call getLocalCoords before maybe make get Element pointer

    double area, *x = pl.data();

    area = getAreaTri2D(x, x + 2, x + 4) + getAreaTri2D(x + 2, x + 4, x);

    // area = abs(area);
    // cout << area << "\n"; getchar();

    if (area < 0)
    {
        throw std::runtime_error("Element::getElemSize: Negative area detected. Check element " +
                                 std::to_string(e) + " for invalid node ordering.");
    }

    return sqrt(area);
}

inline double Element::getElemSizeTri2D(int e)
{
    double area, *x = pl.data();

    area = getAreaTri2D(x, x + 2, x + 4);

    if (area < 0)
    {
        throw std::runtime_error("Element::getElemSize: Negative area detected. Check element " +
                                 std::to_string(e) + " for invalid node ordering.");
    }

    return sqrt(area);
}

inline double Element::getAreaTri2D(double *x0, double *x1, double *x2)
{
    return 0.5 * (x0[0] * (x1[1] - x2[1]) + x1[0] * (x2[1] - x0[1]) + x2[0] * (x0[1] - x1[1]));
}
