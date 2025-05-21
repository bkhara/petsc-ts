
#include <petsc.h>

typedef struct {
  double D0;    // conductivity
} HeatCtx;

static double f_source(double x, double y) {
    return exp( -(x-0.5) * (x-0.5)) * sin(2.0*PETSC_PI*y);
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



  ierr = DMDACreate2d(PETSC_COMM_WORLD,
      DM_BOUNDARY_NONE, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_STAR,
      5,4,PETSC_DECIDE,PETSC_DECIDE,  
      1,1,                            
      NULL,NULL,&da); CHKERRQ(ierr);
  ierr = DMSetFromOptions(da); CHKERRQ(ierr);
  ierr = DMSetUp(da); CHKERRQ(ierr);
  ierr = DMCreateGlobalVector(da,&u); CHKERRQ(ierr);

  ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
  ierr = TSSetProblemType(ts,TS_NONLINEAR); CHKERRQ(ierr);
  ierr = TSSetDM(ts,da); CHKERRQ(ierr);
  ierr = TSSetApplicationContext(ts,&user); CHKERRQ(ierr);
  ierr = DMDATSSetRHSFunctionLocal(da,INSERT_VALUES,
           (DMDATSRHSFunctionLocal)FormRHSFunctionLocal,&user); CHKERRQ(ierr);
  ierr = DMDATSSetRHSJacobianLocal(da,
           (DMDATSRHSJacobianLocal)FormRHSJacobianLocal,&user); CHKERRQ(ierr);

  ierr = TSSetType(ts,TSBDF); CHKERRQ(ierr);
  ierr = TSSetTime(ts,0.0); CHKERRQ(ierr);
  ierr = TSSetMaxTime(ts,0.1); CHKERRQ(ierr);
  ierr = TSSetTimeStep(ts,0.001); CHKERRQ(ierr);
  ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
  ierr = TSSetFromOptions(ts);CHKERRQ(ierr);



  ierr = VecSet(u,0.0); CHKERRQ(ierr);   
  ierr = TSSolve(ts,u); CHKERRQ(ierr);

  VecDestroy(&u);  TSDestroy(&ts);  DMDestroy(&da);
  return PetscFinalize();
}

PetscErrorCode Spacings(DMDALocalInfo *info, double *hx, double *hy) {
    if (hx)  *hx = 1.0 / (double)(info->mx-1);
    if (hy)  *hy = 1.0 / (double)(info->my);  
    return 0;
}




PetscErrorCode FormRHSFunctionLocal(DMDALocalInfo *info,
                                    double t, double **au,
                                    double **aG, HeatCtx *user) {
  PetscErrorCode ierr;
  int      i, j, mx = info->mx;
  double   hx, hy, x, y, ul, ur, uxx, uyy;

  ierr = Spacings(info,&hx,&hy); CHKERRQ(ierr);
  for (j = info->ys; j < info->ys + info->ym; j++) {
      y = hy * j;
      for (i = info->xs; i < info->xs + info->xm; i++) {
          x = hx * i;

          ul = (i == 0) ? au[j][i+1] + 2.0 * hx * gamma_neumann(y)
                        : au[j][i-1];
          ur = (i == mx-1) ? au[j][i-1] : au[j][i+1];
          uxx = (ul - 2.0 * au[j][i]+ ur) / (hx*hx);

          uyy = (au[j-1][i] - 2.0 * au[j][i]+ au[j+1][i]) / (hy*hy);
          aG[j][i] = user->D0 * (uxx + uyy) + f_source(x,y);
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
        row.j = j;  col[0].j = j;
        for (i = info->xs; i < info->xs+info->xm; i++) {

            row.i = i;
            col[0].i = i;
            v[0] = - 2.0 * D * (1.0 / hx2 + 1.0 / hy2);
            col[1].j = j-1;  col[1].i = i;    v[1] = D / hy2;
            col[2].j = j+1;  col[2].i = i;    v[2] = D / hy2;
            col[3].j = j;    col[3].i = i-1;  v[3] = D / hx2;
            col[4].j = j;    col[4].i = i+1;  v[4] = D / hx2;
            ncols = 5;

            if (i == 0) {
                ncols = 4;
                col[3].j = j;  col[3].i = i+1;  v[3] = 2.0 * D / hx2;
            } else if (i == info->mx-1) {
                ncols = 4;
                col[3].j = j;  col[3].i = i-1;  v[3] = 2.0 * D / hx2;
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


