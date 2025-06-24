#define POSTPROCESS

#include "poisson.hpp"  

int main(int argc, char *argv[])
{
    clock_t time_init = clock();
    double time_tot;

    Poisson fem;

    fem.setDimension(2);
    fem.setDof(1);
    
    fem.createBlockMesh2D(0, 1, 0, 1, 16, 16);

    fem.setNgp(4);

    fem.setProp("f", 1.);
    fem.setProp("mu", 1.);

    vector<int> points;

    fem.selectPointsOn(points, 0.0, 0.5, 0.0);  // left
    fem.selectPointsOn(points, 1.0, 0.5, 0.0);  // right
    fem.selectPointsOn(points, 0.5, 1.0, 0.0);  // top

    fem.prescribeDofs(points, 0, 0.);
    points.clear();

    fem.selectPointsOn(points, 0.5, 0.0, 0.0);  // bot
    fem.prescribeDofs(points, 0, 1.);
    points.clear();

    fem.prepareSolver();

    char fname[100];

#ifdef POSTPROCESS
    fem.addVtkMesh();
#endif

    fem.solve();

#ifdef POSTPROCESS
    fem.addVtkSol(0, 0, "phi");

    sprintf(fname, "out/poisson.vtu");
    // sprintf(fname, "out/file_%03d.vtu", fem.it);
    fem.writeVtk(fname);
#endif
}
