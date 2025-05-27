//
// Created by khara on 5/25/25.
//

#include <petsc.h>
#include <talyfem/talyfem.h>

// this is a time-parallel solver for the van der Pol equation
// u'' - mu * (1-u^2) * u' + u = f(t)
// or written in first order format
// u' = v
// v' = mu * (1-u^2) * v - u + f(t)
// we pose a BVP in time
// with the initial condition specified at t = 0
// and no condition at t = T_final

static double Lx = 100.0;
typedef struct {
    double u, v;
} Field;

typedef struct {
    double D0;    // conductivity
    double mu;
} HeatCtx;

typedef struct {
    double l;
    double m;
    double n;
} wave_numbers;

static wave_numbers wn = {0.2, 3, 2};

static double u_exact(double t, double x, double y) {
    double l = wn.l;
    double m = wn.m;
    double n = wn.n;
    return exp(-l * t) * sin(m*PETSC_PI*x) * sin(n*PETSC_PI*y);
}

static double f_source(double t, double x, double y) {
    double l = wn.l;
    double m = wn.m;
    double n = wn.n;
    double u = u_exact(t, x, y);
    return (-l + (m*m + n*n) * PETSC_PI * PETSC_PI) * u;
}

static double gamma_neumann(double y) {
    return sin(4.0 * PETSC_PI * y);
}

PetscErrorCode Spacings(DMDALocalInfo *info, double *hx, double *hy) {
    if (hx)  *hx = Lx / (double)(info->mx-1);
    if (hy)  *hy = 1.0 / (double)(info->my-1);
    return 0;
}

PetscErrorCode FormRHSFunctionLocal(DMDALocalInfo *info, double t, Field *au, Field *aG, const HeatCtx *user) {
    PetscErrorCode ierr;
    int      i, j, mx = info->mx;
    double   hx, hy, y, uc, ul, ur, ud, uu, uxx, uyy;

    const double mu = user->mu;

    ierr = Spacings(info,&hx,&hy); CHKERRQ(ierr);
    for (i = info->xs; i < info->xs + info->xm; i++) {
        double const x = hx * i;

        double u = au[i].u;
        double v = au[i].v;

        if (i == 0) {
            aG[i].u = -(u - 0.2);
            aG[i].v = -(v - 0.1);
        } else {
            double u_dot = (au[i].u - au[i-1].u) / hx;
            double v_dot = (au[i].v - au[i-1].v) / hx;
            aG[i].u = -(u_dot - v);
            aG[i].v = -(v_dot - mu * (1 - u*u) * v + u);
        }
    }

    return 0;
}

PetscErrorCode FormRHSJacobianLocal(DMDALocalInfo *info, double t, Field *au, Mat J, Mat P, const HeatCtx *user) {
    PetscErrorCode ierr;
    constexpr int max_values = 4;
    int            i, j, ncols;
    const double   D = user->D0;
    double         hx, hy, hx2, hy2, v[max_values];
    MatStencil     col[max_values],row;

    ierr = Spacings(info,&hx,&hy); CHKERRQ(ierr);
    hx2 = hx * hx;  hy2 = hy * hy;
    const double mu = user->mu;
    for (i = info->xs; i < info->xs+info->xm; i++) {
        const double un = au[i].u;
        const double vn = au[i].v;
        if (i == 0) {
            ncols = 1;
            row.i = i; row.c = 0;
            col[0].i = i; col[0].c = 0;
            v[0] = -1;
            ierr = MatSetValuesStencil(P,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);

            ncols = 1;
            row.i = i; row.c = 1;
            col[0].i = i; col[0].c = 1;
            v[0] = -1;
            ierr = MatSetValuesStencil(P,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);
        } else {
            // column numbers: u{n-1} = 0, v{n-1} = 1, u{n} = 2, v{n} = 3
            ncols = 4;
            row.i = i; row.c = 0;
            col[0].i = i-1; col[0].c = 0;
            col[1].i = i-1; col[1].c = 1;
            col[2].i = i; col[2].c = 0;
            col[3].i = i; col[3].c = 1;
            v[0] = -(-1.0/hx);
            v[1] = -(0.0);
            v[2] = -(1.0/hx);
            v[3] = -(-1.0);
            ierr = MatSetValuesStencil(P,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);

            ncols = 4;
            row.i = i; row.c = 1;
            col[0].i = i-1; col[0].c = 0;
            col[1].i = i-1; col[1].c = 1;
            col[2].i = i; col[2].c = 0;
            col[3].i = i; col[3].c = 1;
            v[0] = -(0.0);
            v[1] = -(-1.0/hx);
            v[2] = -(2 * mu * un * vn + 1);
            v[3] = -(1/hx - mu *(1 - un*un));
            ierr = MatSetValuesStencil(P,1,&row,ncols,col,v,INSERT_VALUES); CHKERRQ(ierr);
        }
    }

    ierr = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    ierr = MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    if (J != P) {
        ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
        ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
    }
    return 0;
}

// PetscErrorCode ComputeInitialGuess(SNES snes, Vec u, HeatCtx* user) {
//     PetscErrorCode ierr;
//     DMDALocalInfo  info;
//     int            i,j;
//     double         sx,sy;
//     DMDACoor2d     *aC;
//     Field          *aY; // 2D array

//     ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);

//     ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
//     for (i = info.xs; i < info.xs+info.xm; i++)
//     {
//         aY[i].u = 0.2;
//         aY[i].v = 0.1;
//     }
//     ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
//     return 0;
// }

PetscErrorCode InitialState(DM da, Vec Y, HeatCtx* user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    int            i,j;
    double         sx,sy;
    DMDACoor2d     *aC;
    Field          *aY; // 2D array

    ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);

    ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
    for (i = info.xs; i < info.xs+info.xm; i++)
    {
        aY[i].u = 0.2;
        aY[i].v = 0.1;
    }
    ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode WriteSolution(DM da, Vec Y, HeatCtx* user) {
    PetscErrorCode ierr;
    DMDALocalInfo  info;
    int            i,j;
    double         sx,sy;
    DMDACoor2d     *aC;
    Field          *aY; // 2D array

    ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);

    double hx, hy;
    Spacings(&info, &hx, &hy);

    ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
    if (GetMPIRank() == 0) {
        std::ofstream outfile("solution.txt");
        for (i = info.xs; i < info.xs+info.xm; i++)
        {
            double x = i * hx;
            outfile << x << "," << aY[i].u << "," << aY[i].v << "\n";
        }
    }
    ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
    return 0;
}

PetscErrorCode PrintDMDALocalInfo(DMDALocalInfo* info) {
    PrintInfo("Printing DMDA local info:");
    PrintInfo("info->dim = ", info->dim);
    PrintInfo("info->dof = ", info->dof);
    PrintInfo("info->sw = ", info->sw);
    PrintInfo("info->mx, info->my, info->mz = ", info->mx, " ", info->my, " ", info->mz);
    PrintInfo("info->xs, info->ys, info->zs = ", info->xs, " ", info->ys, " ", info->zs);
    PrintInfo("info->xm, info->ym, info->zm = ", info->xm, " ", info->ym, " ", info->zm);
    // PetscSynchronizedPrintf(PETSC_COMM_WORLD, "info->xs = %" PetscInt_FMT ", info->ys = %" PetscInt_FMT ", info->zs = %" PetscInt_FMT "\n", info->xs, info->ys, info->zs);

    double hx, hy;
    Spacings(info, &hx, &hy);
    PrintInfo("hx = ", hx, ", hy = ", hy);
    PrintInfo("num proc = ", GetMPISize());
    return 0;
}


int main(int argc,char **argv)
{
    PetscErrorCode ierr;
    HeatCtx        user;
    TS             ts;
    SNES           snes;
    Mat            J;
    Vec            u, r;
    DM             da;
    DMDALocalInfo  info;

    KSP ksp;


    PetscInitialize(&argc,&argv,nullptr,"heat equation 2D");

    user.D0  = 1.0;
    user.mu = 0.25;

    PetscInt Nx = 250;
    PetscInt ndof = 2;
    PetscInt stencil_width = 1;

    ierr = DMDACreate1d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, Nx, ndof, stencil_width, nullptr, &da); CHKERRQ(ierr);
    ierr = DMSetFromOptions(da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,0,"u"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,1,"v"); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&u); CHKERRQ(ierr);
    // ierr = DMCreateMatrix(da, &J); CHKERRQ(ierr);
    // ierr = VecDuplicate(u, &r); CHKERRQ(ierr);

    DMDAGetLocalInfo(da, &info);
    PrintDMDALocalInfo(&info);

    ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
    ierr = TSSetProblemType(ts,TS_NONLINEAR); CHKERRQ(ierr);
    ierr = TSSetDM(ts,da); CHKERRQ(ierr);
    ierr = TSSetApplicationContext(ts,&user); CHKERRQ(ierr);
    ierr = DMDATSSetRHSFunctionLocal(da,INSERT_VALUES,
           reinterpret_cast<DMDATSRHSFunctionLocal>(FormRHSFunctionLocal),&user); CHKERRQ(ierr);
    ierr = DMDATSSetRHSJacobianLocal(da,
           reinterpret_cast<DMDATSRHSJacobianLocal>(FormRHSJacobianLocal),&user); CHKERRQ(ierr);

    ierr = TSSetType(ts,TSPSEUDO); CHKERRQ(ierr);
    ierr = TSSetTime(ts,0.0); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts,1e8); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts,0.001); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    ierr = TSPseudoSetTimeStep(ts, TSPseudoTimeStepDefault, 0); CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

    ierr = TSGetSNES(ts, &snes); CHKERRQ(ierr);


    ierr = InitialState(da,u,&user); CHKERRQ(ierr);
    ierr = TSSolve(ts,u); CHKERRQ(ierr);

    WriteSolution(da, u, &user);

    VecDestroy(&u);
    TSDestroy(&ts);
    DMDestroy(&da);

    return PetscFinalize();
}