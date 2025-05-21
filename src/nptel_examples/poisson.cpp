//
// Created by khara on 5/21/25.
//

#include <petsc.h>

PetscErrorCode u_exact(DM da, Vec uexact) {
    DMDALocalInfo info;
    int i, j;
    double hx, hy, x, y, **auexact;

    DMDAGetLocalInfo(da, &info);
    DMDAVecGetArray(da, uexact, &auexact);

    hx = 1.0 / (info.mx - 1);
    hy = 1.0 / (info.my - 1);
    for (j = info.ys; j < (info.ys + info.ym); j++) {
        y = j * hy;
        for (i = info.xs; i < (info.xs + info.xm); i++) {
            x = i * hx;
            auexact[j][i] = (x * x - x * x * x * x) * (y * y * y * y - y * y);
        }
    }

    DMDAVecRestoreArray(da, uexact, &auexact);
    return 0;
}


PetscErrorCode formMatrix(DM da, Mat A) {
    DMDALocalInfo info;
    int i, j, ncols;
    double v[5];
    MatStencil row, col[5];

    DMDAGetLocalInfo(da, &info);

    for (j = info.ys; j < (info.ys + info.ym); j++) {
        for (i = info.xs; i < (info.xs + info.xm); i++) {
            row.j = j;
            row.i = i;
            col[0].j = j;
            col[0].i = i;

            ncols = 1;

            if (i == 0 || i == info.mx - 1 || j == 0 || j == info.my - 1) {
                v[0] = 1.0;
            }
            else {
                v[0] = 4;
                if (i - 1 > 0) {
                    col[ncols].j = j;
                    col[ncols].i = i - 1;
                    v[ncols] = -1;
                    ncols = ncols + 1;
                }
                if (i + 1 < info.mx - 1) {
                    col[ncols].j = j;
                    col[ncols].i = i + 1;
                    v[ncols] = -1;
                    ncols = ncols + 1;
                }
                if (j - 1 > 0) {
                    col[ncols].j = j - 1;
                    col[ncols].i = i;
                    v[ncols] = -1;
                    ncols = ncols + 1;
                }
                if (j + 1 < info.my - 1) {
                    col[ncols].j = j + 1;
                    col[ncols].i = i;
                    v[ncols] = -1;
                    ncols = ncols + 1;
                }
            }


            MatSetValuesStencil(A, 1, &row, ncols, col, v, INSERT_VALUES);
        } // x loop
    } // y loop end

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    return 0;
}

PetscErrorCode formRHS(DM da, Vec b) {
    int i, j;
    double h, x, y, f, **ab;
    DMDALocalInfo info;

    DMDAGetLocalInfo(da, &info);
    h = 1.0 / (info.mx - 1);
    DMDAVecGetArray(da, b, &ab);
    for (j = info.ys; j < (info.ys + info.ym); j++) {
        y = j * h;
        for (i = info.xs; i < (info.xs + info.xm); i++) {
            x = i * h;
            f = 2 * ((1 - 6 * x * x) * y * y * (1 - y * y) + (1 - 6 * y * y) * x * x * (1 - x * x));
            ab[j][i] = f * h * h;
        }
    }

    DMDAVecRestoreArray(da, b, &ab);
    return 0;
}

int main(int argc, char** argv) {
    DM da;
    Mat A;
    Vec b, u, uexact;
    KSP ksp;
    double err;
    DMDALocalInfo info; // Data structure to hold information about the array

    PetscInitialize(&argc, &argv, NULL, "Solve Poisson equation in 2D");

    DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR, 9, 9, PETSC_DECIDE,
                 PETSC_DECIDE, 1, 1, NULL, NULL, &da);

    DMSetFromOptions(da);
    DMSetUp(da);
    DMCreateMatrix(da, &A);
    MatSetFromOptions(A);

    DMCreateGlobalVector(da, &b);
    VecDuplicate(b, &u);
    VecDuplicate(b, &uexact);

    // form Exact solution
    u_exact(da, uexact);
    // form the matrix A
    formMatrix(da, A);
    // form the RHS bar
    formRHS(da, b);

    KSPCreate(PETSC_COMM_WORLD, &ksp);
    KSPSetOperators(ksp, A, A);
    KSPSetFromOptions(ksp);
    KSPSolve(ksp, b, u);

    VecAXPY(u, -1, uexact); // axpy
    VecNorm(u, NORM_INFINITY, &err);
    DMDAGetLocalInfo(da, &info);
    PetscPrintf(PETSC_COMM_WORLD, "grid size: %d \t error: %e\n", info.mx, err);

    VecDestroy(&u);
    VecDestroy(&uexact);
    VecDestroy(&b);
    MatDestroy(&A);
    KSPDestroy(&ksp);
    DMDestroy(&da);

    return PetscFinalize();
}
