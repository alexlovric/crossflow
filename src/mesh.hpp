#pragma once

#include <algorithm>
#include <fstream>
#include <sstream>

#include "constants.hpp"
#include "fem.hpp"

/**
 * @class Mesh
 * @brief Handles mesh operations including reading, creating, and processing
 * mesh data
 *
 * Provides functionality for working with mesh files, creating structured
 * meshes, and querying mesh properties. Supports both custom formats and GMSH
 * formats.
 */
class Mesh : public virtual Fem
{
   public:
    Mesh(void) : Fem() {}

    virtual ~Mesh() {}

    /**
     * @brief Reads mesh data from a file
     * @param mesh_name Path to the mesh file
     * @param first Index offset for node numbering (0 or 1 based)
     */
    virtual void readMesh(const char *mesh_name, int first);

    /**
     * @brief Creates a structured 2D rectangular mesh
     * @param xmin Minimum x-coordinate
     * @param xmax Maximum x-coordinate
     * @param ymin Minimum y-coordinate
     * @param ymax Maximum y-coordinate
     * @param nx Number of divisions in x-direction
     * @param ny Number of divisions in y-direction
     */
    virtual void createBlockMesh2D(double xmin, double xmax, double ymin, double ymax, int nx,
                                   int ny);

    /**
     * @brief Identifies and stores boundary point connections for 2D meshes
     */
    virtual void getBoundaryPtsConnect2D(void);

    /**
     * @brief Gets edges connected to a specific element
     * @param edge Output vector of edge indices
     * @param e Index of the element to query
     */
    virtual void getEdges(vector<int> &edge, int e);

    /**
     * @brief Selects points based on spatial coordinates
     * @param select Output vector of selected point indices
     * @param xp X-coordinate threshold
     * @param yp Y-coordinate threshold
     * @param zp Z-coordinate threshold
     */
    virtual void selectPointsOn(vector<int> &select, double xp, double yp, double zp);

    /**
     * @brief Checks if angle between three points is acceptable
     * @param p0 Index of first point
     * @param p1 Index of second point (vertex)
     * @param p2 Index of third point
     * @param cAng Angle threshold in radians
     * @return true if angle is acceptable, false otherwise
     */
    virtual bool acceptableAngle(int p0, int p1, int p2, double cAng);

    /// Vector storing boundary point connections
    vector<vector<int>> boundaryPtsConnect;

#ifdef GMSH
    /**
     * @brief Reads a GMSH format mesh file
     * @param mesh_name Path to the GMSH mesh file
     */
    virtual void readGmsh(const char *mesh_name);

    /**
     * @brief Selects points on a specific boundary
     * @param select Output vector of selected point indices
     * @param boundaryNum Index of the boundary to query
     */
    virtual void selectPointsOnBoundary(vector<int> &select, int boundaryNum);

    /// Vector storing points with associated boundary indices for GMSH
    vector<int> pb;
#endif
};

/**
 * @brief Reads mesh data from a file
 *
 * Reads a mesh file and populates the mesh data structures.
 *
 * @param mesh_name Path to the mesh file
 * @param first Index offset for node numbering (0 or 1 based)
 */
inline void Mesh::readMesh(const char *mesh_name, int first)
{
    printf("Mesh::readMesh -> Reading mesh\n");

    nnp = 0;
    nel = 0;
    nen = 0;

    int ii;

    double jj;

    string row;

    ifstream meshfile(mesh_name);

    if (meshfile.is_open())
    {
        // reading up to nodes
        while (getline(meshfile, row) && row != "$Nodes")
            ;

        // reading nodes
        while (getline(meshfile, row) && row != "$EndNodes")
        {
            stringstream linestream(row);

            for (ii = 0; linestream >> jj; ii++)
            {
                if (ii > ndm)
                    break;
                else if (ii > 0)
                    p.push_back(jj);
            }

            nnp++;
        }

        // reading elements
        while (getline(meshfile, row) && row != "$Elements")
            ;

        getline(meshfile, row);

        stringstream linestream(row);

        for (ii = 0; linestream >> jj; ii++)
            if (ii > 0)
            {
                nen++;
                t.push_back(jj - first);
            }

        nel++;

        while (getline(meshfile, row) && row != "$EndElements")
        {
            stringstream linestream(row);

            for (ii = 0; linestream >> jj; ii++)
            {
                if (ii > 0) t.push_back(jj - first);
            }

            nel++;
        }

        meshfile.close();

        // get boundary points
        getBoundaryPtsConnect2D();

        // set global id and u
        initialiseVars();
    }

    else
    {
        printf(
            "ERROR: \t-> Mesh::readMesh\n"
            "\t-> Check if the file is in current directory.\n\n");
        exit(1);
    }
}

inline void Mesh::createBlockMesh2D(double xmin, double xmax, double ymin, double ymax, int nex,
                                     int ney)
{
    printf("Mesh::createBlockMesh2D : Creating 2D block mesh\n");

    // setting elements
    nen = 4;  // linear quads/blocks

    nel = nex * ney;

    t.resize(nen * nel);

    int npx = nex + 1, npy = ney + 1, ii, jj, kk;

    for (jj = 0; jj < ney; jj++)
    {
        for (ii = 0; ii < nex; ii++)
        {
            kk = (jj * nex + ii) * nen;

            t[kk] = jj * npx + ii;
            t[kk + 1] = t[kk] + 1;
            t[kk + 2] = t[kk + 1] + npx;
            t[kk + 3] = t[kk + 2] - 1;
        }
    }

    // setting nodes
    nnp = npx * npy;

    p.resize(nnp * ndm);

    double xc, yc, dx = (xmax - xmin) / nex, dy = (ymax - ymin) / ney;

    yc = ymin;
    for (jj = 0; jj < npy; jj++)
    {
        xc = xmin;
        for (ii = 0; ii < npx; ii++)
        {
            kk = jj * npx + ii;
            p[kk + kk] = xc;
            p[kk + kk + 1] = yc;
            xc += dx;
        }
        yc += dy;
    }

    // get boundary points
    getBoundaryPtsConnect2D();

    // set global id and u
    initialiseVars();
}

inline void Mesh::getEdges(vector<int> &edge, int e)
{
    int ii;

    switch (nen)
    {
            // triangle
        case 3:
        {
            edge.resize(6);

            edge[0] = t[e * 3];
            edge[1] = t[e * 3 + 1];
            edge[2] = t[e * 3 + 1];
            edge[3] = t[e * 3 + 2];
            edge[4] = t[e * 3 + 2];
            edge[5] = t[e * 3];

            break;
        }

            // quad
        case 4:
        {
            edge.resize(8);

            edge[0] = t[e * 4];
            edge[1] = t[e * 4 + 1];
            edge[2] = t[e * 4 + 1];
            edge[3] = t[e * 4 + 2];
            edge[4] = t[e * 4 + 2];
            edge[5] = t[e * 4 + 3];
            edge[6] = t[e * 4 + 3];
            edge[7] = t[e * 4];

            break;
        }

        default:
        {
            printf(
                "ERROR : Mesh::getEdges\n"
                "\t-> Not yet implemented for other nen.\n\n");
            exit(1);
        }
    }
}

inline void Mesh::getBoundaryPtsConnect2D(void)
{
    vector<vector<int>> tmp;
    tmp.resize(nnp);

    vector<int> edge;

    int ii, jj, p1, p2, c1, c2;

    boundaryPtsConnect.resize(nnp);

    for (ii = 0; ii < nel; ii++)
    {
        getEdges(edge, ii);

        for (jj = 0; jj < edge.size(); jj += 2)
        {
            p1 = edge[jj];
            p2 = edge[jj + 1];

            c1 = tmp[p2].size();

            tmp[p2].erase(remove(tmp[p2].begin(), tmp[p2].end(), p1), tmp[p2].end());

            c2 = tmp[p2].size();

            if (c2 - c1 == 0) tmp[p1].push_back(p2);
        }
    }

    for (ii = 0; ii < tmp.size(); ii++)
    {
        for (jj = 0; jj < tmp[ii].size(); jj++)
        {
            boundaryPtsConnect[ii].push_back(tmp[ii][jj]);
            boundaryPtsConnect[tmp[ii][jj]].push_back(ii);
        }
    }
}

inline void Mesh::selectPointsOn(vector<int> &select, double xp, double yp, double zp)
{
    int n, ii, jj, p0, p1, p2, m = select.size();
    double d, dmin, r[3] = {xp, yp, zp}, mxAng = 20., cAng = cos(mxAng / 180. * pi);

    // getting node closest to xp,yp,zp
    for (jj = 0; jj < ndm; jj++)
        dmin += (p[jj] - r[jj]) * (p[jj] - r[jj]);

    for (ii = 0; ii < nnp; ii++)
    {
        d = 0.;
        for (jj = 0; jj < ndm; jj++)
            d += (p[jj + ii * ndm] - r[jj]) * (p[jj + ii * ndm] - r[jj]);

        if (d < dmin)
        {
            dmin = d;
            n = ii;
        }
    }

    // finding points adjactent
    vector<int> queue;
    vector<bool> selected(nnp, false), inQueue(nnp, false);

    for (ii = 0; ii < m; ii++)
        selected[select[ii]] = true;

    if (!selected[n])
    {
        select.push_back(n);
        selected[n] = true;
    }

    if (ndm == 2)
    {
        queue.push_back(n);
        inQueue[n] = true;

        while (queue.size() > 0)
        {
            p0 = queue[0];

            for (auto kk = boundaryPtsConnect[p0].begin(); next(kk) != boundaryPtsConnect[p0].end();
                 ++kk)
            {
                p1 = *kk;
                p2 = *next(kk);

                if (acceptableAngle(p0, p1, p2, cAng))
                {
                    if (!selected[p1])
                    {
                        select.push_back(p1);
                        selected[p1] = true;
                    }
                    if (!selected[p2])
                    {
                        select.push_back(p2);
                        selected[p2] = true;
                    }

                    if (!inQueue[p1])
                    {
                        queue.push_back(p1);
                        inQueue[p1] = true;
                    }
                    if (!inQueue[p2])
                    {
                        queue.push_back(p2);
                        inQueue[p2] = true;
                    }
                }
            }
            queue.erase(queue.begin());
        }
    }

    else if (ndm == 3)
    {
        // TODO: do this
    }
}

inline bool Mesh::acceptableAngle(int p0, int p1, int p2, double cAng)
{
    double x0[2], x1[2], x2[2];

    for (int ii = 0; ii < 2; ii++)
    {
        x0[ii] = p[p1 * ndm + ii];
        x1[ii] = p[p0 * ndm + ii];
        x2[ii] = p[p2 * ndm + ii];
    }

    double dx1[2] = {x1[0] - x0[0], x1[1] - x0[1]}, dx2[2] = {x2[0] - x1[0], x2[1] - x1[1]},
           l1 = dx1[0] * dx1[0] + dx1[1] * dx1[1], l2 = dx2[0] * dx2[0] + dx2[1] * dx2[1];

    // cout << (dx1[0]*dx2[0]+dx1[1]*dx2[1]) / sqrt(l1 * l2);

    if ((dx1[0] * dx2[0] + dx1[1] * dx2[1]) / sqrt(l1 * l2) > cAng) return true;

    return false;
}

#ifdef GMSH

inline void Mesh::readGmsh(const char *mesh_name)
{
    printf("Mesh::readGmsh : Reading gmsh format\n");

    nen = 0;

    int ii, jj, e1 = 20, e2 = 0, epos = 3, tpos = 4, first = 1;

    bool physNames = false;

    double kk;

    string row;

    ifstream meshfile(mesh_name);

    if (ndm == 2)
    {
        if (meshfile.is_open())
        {
            // reading up to nodes
            while (getline(meshfile, row) && row != "$Nodes")
            {
                if (row == "$PhysicalNames")
                {
                    physNames = true;
                }
            }

            getline(meshfile, row);

            nnp = atoi(row.c_str());

            pb.resize(nnp);

            // reading nodes
            while (getline(meshfile, row) && row != "$EndNodes")
            {
                stringstream linestream(row);

                for (ii = 0; linestream >> kk; ii++)
                {
                    if (ii > ndm)
                        break;
                    else if (ii > 0)
                        p.push_back(kk);
                }
            }

            // reading up to elements
            while (getline(meshfile, row) && row != "$Elements")
                ;

            getline(meshfile, row);

            nel = atoi(row.c_str());

            if (physNames)
            {
                // reading edges
                while (getline(meshfile, row))
                {
                    stringstream linestream1(row), linestream2(row);

                    e2 = 0;
                    for (ii = 0; linestream1 >> kk; ii++)
                        e2++;

                    if (e2 > e1) break;
                    e1 = e2;

                    for (ii = 0; linestream2 >> kk; ii++)
                    {
                        if (ii == epos) jj = kk;
                        if (ii > tpos) pb[kk - first] = jj;
                    }

                    nel--;
                }
            }

            // reading elements
            stringstream linestream(row);

            for (ii = 0; linestream >> kk; ii++)
            {
                if (ii == epos) jj = kk;
                if (ii > tpos)
                {
                    nen++;
                    t.push_back(kk - first);
                }
            }

            // reading rest of elements
            while (getline(meshfile, row) && row != "$EndElements")
            {
                stringstream linestream(row);

                for (ii = 0; linestream >> kk; ii++)
                {
                    if (ii == epos) jj = kk;
                    if (ii > tpos)
                    {
                        t.push_back(kk - first);
                    }
                }
            }

            meshfile.close();

            // get boundary points
            getBoundaryPtsConnect2D();

            // set global id and u
            initialiseVars();
        }

        else
        {
            printf(
                "ERROR : Mesh::readGmsh\n"
                "\t-> Check if the file is in current directory.\n\n");
            exit(1);
        }
    }
    else
    {
        printf(
            "ERROR : Mesh::readGmsh\n"
            "\t-> 3D not yet implemented.\n\n");
        exit(1);
    }
}

inline void Mesh::selectPointsOnBoundary(vector<int> &select, int boundaryNum)
{
    int ii;

    for (ii = 0; ii < nnp; ii++)
    {
        if (pb[ii] == boundaryNum) select.push_back(ii);
    }
}

#endif
