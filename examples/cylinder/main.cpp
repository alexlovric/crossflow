#define POSTPROCESS

#include "cflow.hpp"

//==============================================================================
int main(int argc, char *argv[])
{
    clock_t time_init = clock();
    double time_tot;

    // ############################################
    //  preprocessing stage
    Cflow fem;

    fem.setDimension(2);
    fem.setDof(3);

    fem.readMesh("cylinder.msh", 1);
    // fem.createBlockMesh2D(0,1,0,1,64,64);

    fem.setNgp(3);
    fem.setProp("rho", 1.0);
    fem.setProp("mu", 0.01);
    fem.setProp("linear", false);

    fem.setTime("steady", false);
    fem.setTime("rhoInf", 0.5);
    fem.setTime("dt", 0.0001);
    fem.setTime("tFin", 200.);

    fem.setTime("adaptive", true);
    fem.setTime("dtMax", 0.2);
    fem.setTime("dtMin", 0.0001);
    fem.setTime("theta", 0.8);

    vector<int> points, points2;

    // slip
    fem.selectPointsOn(points, 0.0, -15., 0.0);  // bot
    fem.selectPointsOn(points, 0.0, 15., 0.0);   // top

    fem.prescribeDofs(points, 1, 0);

    points.clear();

    // velocity inlet
    fem.selectPointsOn(points, -10., 0.0, 0.0);  // left

    fem.prescribeDofs(points, 0, 1.);
    fem.prescribeDofs(points, 1, 0);

    points.clear();

    // pressure outlet
    fem.selectPointsOn(points, 20., 0.0, 0.0);  // right

    points2.push_back(52);

    fem.prescribeDofs(points2, 2, 0);

    points.clear();
    points2.clear();

    fem.selectPointsOn(points, 0.05, 0.0, 0.0);  // left

    fem.prescribeDofs(points, 0, 0);
    fem.prescribeDofs(points, 1, 0);

    points.clear();

    // ############################################
    //  processing stage

    fem.prepareSolver();

    char fname[100];

#ifdef POSTPROCESS
    fem.addVtkMesh();
    // fem.writeVtk(fname);
#endif

    // Time loop
    for (fem.it = 0; fem.time < fem.tFin; fem.it++)
    {
        fem.time += fem.dt;

        printf("\ntime: %1.5f", fem.time);

        fem.solve();

        fem.un = fem.u;
        fem.dun = fem.du;

        fem.tpr = fem.time;

#ifdef POSTPROCESS
        fem.addVtkSol(0, 1, "U");
        fem.addVtkSol(2, 0, "p");

        sprintf(fname, "out/file_%03d.vtu", fem.it);
        fem.writeVtk(fname);
#endif
    }

    time_tot = (clock() - time_init) / (double)CLOCKS_PER_SEC;
    printf("\nTime taken: %.6f\n\n", time_tot);
}
