
#include <petsc.h>
#include <talyfem/talyfem.h>

// this program solves u_t = -laplacian u + f
// with all zero Dirichlet boundary conditions
// currently cannot handle nonzero boundary condition
// exact solution is found in u_exact function

typedef struct {
  double D0;    // conductivity
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
    return (1.0 - exp(-l * t)) * sin(m*PETSC_PI*x) * sin(n*PETSC_PI*y);
}

static double f_source(double t, double x, double y) {
    double l = wn.l;
    double m = wn.m;
    double n = wn.n;
    double pi = PETSC_PI;
    return (l * exp(-l*t) + (m*m + n*n)*pi*pi*(1.0 - exp(-l * t))) * sin(m*PETSC_PI*x) * sin(n*PETSC_PI*y);
}

static double gamma_neumann(double y) {
    return sin(4.0 * PETSC_PI * y);
}

// Function prototypes defined here
extern PetscErrorCode Spacings(DMDALocalInfo*, double*, double*);

extern PetscErrorCode FormRHSFunctionLocal(DMDALocalInfo*, double, double**,
                                           double**, HeatCtx*);
extern PetscErrorCode FormRHSJacobianLocal(DMDALocalInfo*, double, double**,
                                           Mat, Mat, HeatCtx*);

extern PetscErrorCode PrintDMDALocalInfo(DMDALocalInfo*);

int main(int argc,char **argv)
{
    PetscErrorCode ierr;
    HeatCtx        user;
    TS             ts;
    Vec            u;
    DM             da;
    DMDALocalInfo  info;


    PetscInitialize(&argc,&argv,NULL,"heat equation 2D");

    user.D0  = 1.0;

    PetscInt Nx = 5;
    PetscInt Ny = 5;
    PetscInt ndof = 1;
    PetscInt stencil_width = 1;

    ierr = DMDACreate2d(PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE, DM_BOUNDARY_NONE, DMDA_STENCIL_STAR,
      Nx,Ny,PETSC_DECIDE,PETSC_DECIDE,
      ndof,stencil_width,
      NULL,NULL,&da); CHKERRQ(ierr);
    ierr = DMSetFromOptions(da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,0,"u"); CHKERRQ(ierr);
    ierr = DMCreateGlobalVector(da,&u); CHKERRQ(ierr);

    DMDAGetLocalInfo(da, &info);
    PrintDMDALocalInfo(&info);

    ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
    ierr = TSSetProblemType(ts,TS_NONLINEAR); CHKERRQ(ierr);
    ierr = TSSetDM(ts,da); CHKERRQ(ierr);
    ierr = TSSetApplicationContext(ts,&user); CHKERRQ(ierr);
    ierr = DMDATSSetRHSFunctionLocal(da,INSERT_VALUES,
           (DMDATSRHSFunctionLocal)FormRHSFunctionLocal,&user); CHKERRQ(ierr);
    ierr = DMDATSSetRHSJacobianLocal(da,
           (DMDATSRHSJacobianLocal)FormRHSJacobianLocal,&user); CHKERRQ(ierr);

    ierr = TSSetType(ts,TSPSEUDO); CHKERRQ(ierr);
    ierr = TSSetTime(ts,0.0); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts,1e8); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts,0.001); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    ierr = TSPseudoSetTimeStep(ts, TSPseudoTimeStepDefault, 0); CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts);CHKERRQ(ierr);



    ierr = VecSet(u,0.0); CHKERRQ(ierr);
    ierr = TSSolve(ts,u); CHKERRQ(ierr);

    VecDestroy(&u);  TSDestroy(&ts);  DMDestroy(&da);
    return PetscFinalize();
}

PetscErrorCode Spacings(DMDALocalInfo *info, double *hx, double *hy) {
    if (hx)  *hx = 1.0 / (double)(info->mx-1);
    if (hy)  *hy = 1.0 / (double)(info->my-1);
    return 0;
}


PetscErrorCode FormRHSFunctionLocal(DMDALocalInfo *info,
                                    double t, double **au,
                                    double **aG, HeatCtx *user) {
  PetscErrorCode ierr;
  int      i, j, mx = info->mx;
  double   hx, hy, x, y, uc, ul, ur, ud, uu, uxx, uyy;

  ierr = Spacings(info,&hx,&hy); CHKERRQ(ierr);
  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j;
      for (i = info->xs; i < info->xs + info->xm; i++) {
          x = hx * i;

          double g = u_exact(t, x, y);

          if (i == 0 || i == info->mx - 1 || j == 0 || j == info->my - 1) {
              aG[j][i] = 0;
          } else {
              uc = au[j][i];
              ul = au[j][i-1];
              ur = au[j][i+1];
              ud = au[j-1][i];
              uu = au[j+1][i];

              uxx = (ul - 2.0 * uc + ur) / (hx*hx);
              uyy = (ud - 2.0 * uc + uu) / (hy*hy);
              aG[j][i] = user->D0 * (uxx + uyy) + f_source(t, x, y);
          }
      }
  }
  return 0;
}

PetscErrorCode FormRHSJacobianLocal(DMDALocalInfo *info,
                                    double t, double **au,
                                    Mat J, Mat P, HeatCtx *user) {
    PetscErrorCode ierr;
    int            i, j, ncols;
    const double   D = user->D0;
    double         hx, hy, hx2, hy2, v[5];
    MatStencil     col[5],row;

    ierr = Spacings(info,&hx,&hy); CHKERRQ(ierr);
    hx2 = hx * hx;  hy2 = hy * hy;
    for (j = info->ys; j < info->ys+info->ym; j++) {
        row.j = j;
        col[0].j = j;
        for (i = info->xs; i < info->xs+info->xm; i++) {
            if (i == 0 || i == info->mx - 1 || j == 0 || j == info->my - 1) {
                ncols = 1;
                row.i = i;
                col[0].i = i;
                v[0] = 1.;
            } else {
                ncols = 5;
                row.i = i;
                col[0].i = i;
                v[0] = - 2.0 * D * (1.0 / hx2 + 1.0 / hy2);

                col[1].j = j-1;  col[1].i = i;    v[1] = D / hy2;
                col[2].j = j+1;  col[2].i = i;    v[2] = D / hy2;
                col[3].j = j;    col[3].i = i-1;  v[3] = D / hx2;
                col[4].j = j;    col[4].i = i+1;  v[4] = D / hx2;
            }
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


