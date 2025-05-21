#include <iostream>
#include <petsc.h>

static char help[] = "PETSC time stepper example";

PetscErrorCode FormRHS(TS ts, double t, Vec y, Vec g, void * ctx) {
    const double *ay;
    double *ag;

    VecGetArrayRead(y, &ay);
    VecGetArray(g, &ag);

    ag[0] = ay[1];
    ag[1] = -ay[0] + t;

    VecRestoreArrayRead(y, &ay);
    VecRestoreArray(g, &ag);

    return 0;
}

PetscErrorCode SetExact(double t, Vec y) {
    double *ay;
    VecGetArray(y, &ay);
    ay[0] = t - sin(t);
    ay[1] = 1.0 - cos(t);
    VecRestoreArray(y, &ay);

    return 0;
}

int main(int argc, char **args) {
    PetscInitialize(&argc, &args, NULL, help);

    const int N = 2;
    int nsteps;
    double t0 = 0.0, tf = 10.0, dt = 0.1;
    double abserr;
    Vec y, yexact;
    TS ts;

    VecCreate(PETSC_COMM_WORLD, &y);
    VecSetSizes(y, PETSC_DECIDE, N);
    VecSetFromOptions(y);
    VecDuplicate(y, &yexact);

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetProblemType(ts, TS_NONLINEAR);
    TSSetRHSFunction(ts, NULL, FormRHS, NULL);
    TSSetType(ts, TSRK);

    TSSetTime(ts, t0);
    TSSetMaxTime(ts, tf);
    TSSetTimeStep(ts, dt);
    TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP);
    TSSetFromOptions(ts);

    TSGetTime(ts, &t0);
    SetExact(t0, y);
    TSSolve(ts, y);
    TSGetStepNumber(ts, &nsteps);
    TSGetTime(ts, &tf); // fetches exact time at the end of the solver
    SetExact(tf, yexact);
    VecAXPY(y, -1, yexact);
    VecNorm(y, NORM_INFINITY, &abserr);

    PetscPrintf(PETSC_COMM_WORLD, "finaltime: %1.4f\t n_steps: %d, error: %g\n", tf, nsteps, abserr);


    VecDestroy(&y), VecDestroy(&yexact), TSDestroy(&ts);
    PetscFinalize();
    return 0;
}
