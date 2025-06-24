#pragma once

#ifdef POSTPROCESS
#define vtkRenderingCore_AUTOINIT 3(vtkInteractionStyle, vtkRenderingFreeType, vtkRenderingOpenGL2)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL2)

#include "vtkAbstractArray.h"
#include "vtkActor.h"
#include "vtkCell.h"
#include "vtkDataArray.h"
#include "vtkDataSet.h"
#include "vtkDataSetMapper.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkFloatArray.h"
#include "vtkLine.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkProperty.h"
#include "vtkQuad.h"
#include "vtkSmartPointer.h"
#include "vtkSmartPointerBase.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkVertex.h"
#include "vtkXMLUnstructuredGridWriter.h"
#endif

#include "fem.hpp"

/**
 * @class Vtk
 * @brief VTK visualization interface for finite element analysis results.
 *
 * This class provides functionality to convert FEM mesh and solution data into
 * VTK format for visualization. Requires POSTPROCESS macro to be defined.
 */
class Vtk : public virtual Fem
{
   public:
#ifdef POSTPROCESS

    Vtk(void);

    virtual ~Vtk() {}

    /**
     * @brief Converts FEM mesh to VTK unstructured grid.
     *
     * Supports 2D linear triangles and quads. Nodes are mapped to VTK points,
     * elements to VTK cells.
     *
     * @throws std::runtime_error if unsupported element type or dimension
     */
    virtual void addVtkMesh(void);

    /**
     * @brief Writes selected nodes as points to VTK file.
     *
     * @param points Vector of node indices to visualize
     * @param filename Output VTK file path (.vtu)
     */
    virtual void writeVtkPoints(vector<int> &, const char *);

    /**
     * @brief Writes current VTK data to file in binary format.
     *
     * @param filename Output VTK file path (.vtu)
     */
    virtual void writeVtk(const char *);

    /**
     * @brief Adds solution field to VTK output.
     *
     * @param df Degree of freedom index (0-based)
     * @param type Field type (0=scalar, 1=vector)
     * @param name Field name in VTK output
     * @throws std::runtime_error for invalid type
     */
    virtual void addVtkSol(int df, int type, const char *name);

    /// VTK unstructured grid containing mesh and solution data
    vtkSmartPointer<vtkUnstructuredGrid> uGridVtk;

    /// VTK points container for node coordinates
    vtkSmartPointer<vtkPoints> pointsVtk;

    /// VTK XML writer for unstructured grid data
    vtkSmartPointer<vtkXMLUnstructuredGridWriter> writerVtk;

#endif
};

#ifdef POSTPROCESS

// #############################################################################
//  constructor
// #############################################################################
inline Vtk::Vtk(void) : Fem()
{
    uGridVtk = vtkSmartPointer<vtkUnstructuredGrid>::New();
    pointsVtk = vtkSmartPointer<vtkPoints>::New();
    writerVtk = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
}

inline void Vtk::addVtkMesh(void)
{
    printf("Vtk::addVtkMesh : Creating mesh vtk file\n");

    vtkSmartPointer<vtkQuad> quadVTK = vtkSmartPointer<vtkQuad>::New();

    vtkSmartPointer<vtkTriangle> triaVTK = vtkSmartPointer<vtkTriangle>::New();

    int ii, ll, ee;
    double xx, yy;

    vtkIdType pt[8];

    uGridVtk->Reset();
    pointsVtk->Reset();

    if (ndm == 2)
    {
        for (ii = 0; ii < nnp; ii++)
        {
            xx = p[ii * ndm];
            yy = p[ii * ndm + 1];

            pt[0] = pointsVtk->InsertNextPoint(xx, yy, 0.0);
        }

        if (nen == 3)
        {
            for (ee = 0; ee < nel; ee++)
            {
                for (ll = 0; ll < 3; ll++)
                    triaVTK->GetPointIds()->SetId(ll, t[ee * 3 + ll]);

                uGridVtk->InsertNextCell(triaVTK->GetCellType(), triaVTK->GetPointIds());
            }
        }
        else if (nen == 4)  // quad
        {
            for (ee = 0; ee < nel; ee++)
            {
                for (ll = 0; ll < 4; ll++)
                    quadVTK->GetPointIds()->SetId(ll, t[ee * 4 + ll]);

                uGridVtk->InsertNextCell(quadVTK->GetCellType(), quadVTK->GetPointIds());
            }
        }
        else
        {
            printf(
                "ERROR : Vtk::addVtkMesh\n"
                "\t-> Only linear triangles/quads implemented.\n\n");
            exit(1);
        }
    }
    else
    {
        printf(
            "ERROR : Vtk::addVtkMesh\n"
            "\t-> 1D/3D not implemented yet.\n\n");
        exit(1);
    }

    uGridVtk->SetPoints(pointsVtk);
};

inline void Vtk::writeVtkPoints(vector<int> &points, const char *filename)
{
    printf("Vtk::writeVtkPoints : Creating points vtk file\n");

    vtkSmartPointer<vtkVertex> vertexVTK = vtkSmartPointer<vtkVertex>::New();

    vtkIdType pt[8];

    uGridVtk->Reset();
    pointsVtk->Reset();

    int ii;

    if (ndm == 2)
    {
        for (ii = 0; ii < points.size(); ii++)
            pointsVtk->InsertNextPoint(p[points[ii] * ndm], p[points[ii] * ndm + 1], 0.);
    }

    for (ii = 0; ii < points.size(); ii++)
    {
        vertexVTK->GetPointIds()->SetId(0, ii);

        uGridVtk->InsertNextCell(vertexVTK->GetCellType(), vertexVTK->GetPointIds());
    }

    uGridVtk->SetPoints(pointsVtk);

    writerVtk->SetFileName(filename);
    writerVtk->SetInputData(uGridVtk);
    writerVtk->Write();
}

inline void Vtk::writeVtk(const char *filename)
{
    writerVtk->SetFileName(filename);
    writerVtk->SetInputData(uGridVtk);
    writerVtk->SetDataModeToBinary();
    writerVtk->Write();

    printf("\nWritten output file: %s\n", filename);
};

inline void Vtk::addVtkSol(int df, int type, const char *name)
{
    vtkSmartPointer<vtkDoubleArray> solVtk = vtkSmartPointer<vtkDoubleArray>::New();

    solVtk->SetName(name);

    int ii;

    switch (type)
    {
        // scalar
        case 0:
        {
            solVtk->SetNumberOfComponents(1);

            for (ii = 0; ii < nnp; ii++)
                solVtk->InsertNextValue(u[ii * ndf + df]);

            break;
        }

        // vector
        case 1:
        {
            int jj;
            double vec[3] = {0., 0., 0.};

            solVtk->SetNumberOfComponents(3);

            for (ii = 0; ii < nnp; ii++)
            {
                for (jj = 0; jj < ndm; jj++)
                    vec[jj] = u[ii * ndf + df + jj];

                solVtk->InsertNextTuple(vec);
            }

            break;
        }

        default:
        {
            printf(
                "ERROR : Vtk::addVtkSol\n"
                "\t-> Second entry must be 0 - scalar or 1 - vector.\n\n");
            exit(1);
        }
    }

    uGridVtk->GetPointData()->AddArray(solVtk);
}

#endif