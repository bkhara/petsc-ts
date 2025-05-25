#include <fstream>
#include <iostream>
#include <petsc.h>
#include <talyfem/talyfem.h>

static char help[] = "Van Der Pol time stepping";

typedef struct {
    double    mu, amp, freq; // Parameters of the problem
} AppCtx;

PetscErrorCode FormRHS(TS ts, double t, Vec y, Vec g, void * ctx) {
    const double *ay;
    double *ag;
    auto user = static_cast<AppCtx*>(ctx);

    VecGetArrayRead(y, &ay);
    VecGetArray(g, &ag);

    double u = ay[0];
    double v = ay[1];

    ag[0] = v;
    ag[1] = -user->mu * (1 - u*u)*v - u;

    VecRestoreArrayRead(y, &ay);
    VecRestoreArray(g, &ag);

    return 0;
}

PetscErrorCode SetExact(double t, Vec y) {
    double *ay;
    VecGetArray(y, &ay);
    ay[0] = 0.1;
    ay[1] = 0.0;
    VecRestoreArray(y, &ay);

    return 0;
}

PetscErrorCode TSMonitorCallback(TS ts, PetscInt steps, PetscReal time, Vec u, void *ctx) {
    const double *au;
    double *ag;

    VecGetArrayRead(u, &au);

    if (GetMPIRank() == 0) {
        std::ofstream outfile("solution.txt", std::ios::app);
        outfile << time << "," << au[0] << "," << au[1] << "\n";
    }

    return 0;
}

int main(int argc, char **args) {
    PetscInitialize(&argc, &args, NULL, help);

    std::ofstream outfile("solution.txt");

    const int N = 2;
    int nsteps;
    double t0 = 0.0, tf = 10.0, dt = 0.1;
    double abserr;
    Vec y, yexact;
    TS ts;

    AppCtx user;
    user.mu = 0.25;
    user.amp = 0.;
    user.freq = 0.2;

    VecCreate(PETSC_COMM_WORLD, &y);
    VecSetSizes(y, PETSC_DECIDE, N);
    VecSetFromOptions(y);
    VecDuplicate(y, &yexact);

    TSCreate(PETSC_COMM_WORLD, &ts);
    TSSetProblemType(ts, TS_NONLINEAR);
    TSSetRHSFunction(ts, NULL, FormRHS, &user);
    TSMonitorSet(ts, TSMonitorCallback, &user, NULL);
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
