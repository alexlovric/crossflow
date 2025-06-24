#pragma once

#include "element.hpp"
#include "mesh.hpp"
#include "navier_stokes_incomp.hpp"
#include "solver.hpp"
#include "time.hpp"
#include "vtk.hpp"

/**
 * @class Cflow
 * @brief Main problem class combining all FEM components
 */
class Cflow : public Mesh, public Time, public Element, public Solver, public Vtk, public NavierStokesIncomp
{
   public:
    Cflow() : Mesh(), Time(), Element(), Solver(), Vtk(), NavierStokesIncomp() {}

    virtual ~Cflow() = default;
};