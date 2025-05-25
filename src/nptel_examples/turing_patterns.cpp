#include <petsc.h>

typedef struct {
  double u, v;
} Field;

typedef struct {
  double    L, Du, Dv, phi, kappa; // Parameters of the problem
} AppCtx;

PetscErrorCode InitialState(DM da, Vec Y, double noiselevel, AppCtx* user) {
  PetscErrorCode ierr;
  DMDALocalInfo  info;
  int            i,j;
  double         sx,sy;
  const double   ledge = (user->L - 0.5) / 2.0, redge = user->L - ledge;       // declare some analytical function on top the random variables between x = 1.0, to x = 1.5
  DMDACoor2d     **aC;
  Field          **aY; // 2D array

  ierr = VecSet(Y,0.0); CHKERRQ(ierr);
  ierr = VecSetRandom(Y,NULL); CHKERRQ(ierr);
  ierr = VecScale(Y,noiselevel); CHKERRQ(ierr);
	
  ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);
  ierr = DMDAGetCoordinateArray(da,&aC); CHKERRQ(ierr); // ac[j][i].x or ac[j][i].y
  ierr = DMDAVecGetArray(da,Y,&aY); CHKERRQ(ierr);
	
  for (j = info.ys; j < info.ys+info.ym; j++) 
	{
    for (i = info.xs; i < info.xs+info.xm; i++) 
		{
      if ((aC[j][i].x >= ledge) && (aC[j][i].x <= redge) && (aC[j][i].y >= ledge) && (aC[j][i].y <= redge)) 
			{
          sx = sin(4.0 * PETSC_PI * aC[j][i].x);
          sy = sin(4.0 * PETSC_PI * aC[j][i].y);
          aY[j][i].v += 0.5 * sx * sx * sy * sy;
      }
      aY[j][i].u += 1.0 - 2.0 * aY[j][i].v;
    }
  }
  ierr = DMDAVecRestoreArray(da,Y,&aY); CHKERRQ(ierr);
  ierr = DMDARestoreCoordinateArray(da,&aC); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode FormRHSFunctionLocal(DMDALocalInfo *info, double t, Field **aY, Field **aG, AppCtx *user) 
{
  int            i, j;
  double         uv2;

  for (j = info->ys; j < info->ys + info->ym; j++) 
	{
      for (i = info->xs; i < info->xs + info->xm; i++) 
			{
          uv2 = aY[j][i].u * aY[j][i].v * aY[j][i].v;
          aG[j][i].u = - uv2 + user->phi * (1.0 - aY[j][i].u);
          aG[j][i].v = + uv2 - (user->phi + user->kappa) * aY[j][i].v;
      }
  }
  return 0;
}


PetscErrorCode FormRHSJacobianLocal(DMDALocalInfo *info, double t, Field **aY, Mat J, Mat P, AppCtx *user) 
{
    PetscErrorCode ierr;
    int            i, j;
    double         v[2], uv, v2;
    MatStencil     col[2],row;

    for (j = info->ys; j < info->ys+info->ym; j++) 
		{
        row.j = j;  col[0].j = j;  col[1].j = j;
        for (i = info->xs; i < info->xs+info->xm; i++) 
				{
            row.i = i;  col[0].i = i;  col[1].i = i;
            uv = aY[j][i].u * aY[j][i].v;
            v2 = aY[j][i].v * aY[j][i].v;

            row.c = 0;  col[0].c = 0;  col[1].c = 1;
            v[0] = - v2 - user->phi;
            v[1] = - 2.0 * uv;
            ierr = MatSetValuesStencil(P,1,&row,2,col,v,INSERT_VALUES); CHKERRQ(ierr);

            row.c = 1;  col[0].c = 0;  col[1].c = 1;
            v[0] = v2;
            v[1] = 2.0 * uv - (user->phi + user->kappa);
            ierr = MatSetValuesStencil(P,1,&row,2,col,v,INSERT_VALUES); CHKERRQ(ierr);
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


PetscErrorCode FormIFunctionLocal(DMDALocalInfo *info, double t, Field **aY, Field **aYdot, Field **aF, AppCtx *user) 
{
  int            i, j;
  const double   h = user->L / (double)(info->mx),
                 Cu = user->Du / (6.0 * h * h),
                 Cv = user->Dv / (6.0 * h * h);
  double         u, v, lapu, lapv;

  for (j = info->ys; j < info->ys + info->ym; j++) 
	{
      for (i = info->xs; i < info->xs + info->xm; i++) 
			{
          u = aY[j][i].u;
          v = aY[j][i].v;
          lapu =     aY[j+1][i-1].u + 4.0*aY[j+1][i].u +   aY[j+1][i+1].u
                 + 4.0*aY[j][i-1].u -    20.0*u        + 4.0*aY[j][i+1].u
                 +   aY[j-1][i-1].u + 4.0*aY[j-1][i].u +   aY[j-1][i+1].u;
          lapv =     aY[j+1][i-1].v + 4.0*aY[j+1][i].v +   aY[j+1][i+1].v
                 + 4.0*aY[j][i-1].v -    20.0*v        + 4.0*aY[j][i+1].v
                 +   aY[j-1][i-1].v + 4.0*aY[j-1][i].v +   aY[j-1][i+1].v;
          aF[j][i].u = aYdot[j][i].u - Cu * lapu;
          aF[j][i].v = aYdot[j][i].v - Cv * lapv;
      }
  }
  return 0;
}

//     J = (shift) dF/d(dot Y) + dF/dY

PetscErrorCode FormIJacobianLocal(DMDALocalInfo *info,
                   double t, Field **aY, Field **aYdot, double shift,
                   Mat J, Mat P, AppCtx *user) {
    PetscErrorCode ierr;
    int            i, j, s, c;
    const double   h = user->L / (double)(info->mx),
                   Cu = user->Du / (6.0 * h * h),
                   Cv = user->Dv / (6.0 * h * h);
    double         val[9], CC;
    MatStencil     col[9], row;

    for (j = info->ys; j < info->ys + info->ym; j++) {
        row.j = j;
        for (i = info->xs; i < info->xs + info->xm; i++) {
            row.i = i;
            for (c = 0; c < 2; c++) { // u,v equations are c=0,1
                row.c = c; 
                CC = (c == 0) ? Cu : Cv;
                for (s = 0; s < 9; s++)
                    col[s].c = c;
                col[0].i = i;   col[0].j = j;
                val[0] = shift + 20.0 * CC;
                col[1].i = i-1; col[1].j = j;    val[1] = - 4.0 * CC;
                col[2].i = i+1; col[2].j = j;    val[2] = - 4.0 * CC;
                col[3].i = i;   col[3].j = j-1;  val[3] = - 4.0 * CC;
                col[4].i = i;   col[4].j = j+1;  val[4] = - 4.0 * CC;
                col[5].i = i-1; col[5].j = j-1;  val[5] = - CC;
                col[6].i = i-1; col[6].j = j+1;  val[6] = - CC;
                col[7].i = i+1; col[7].j = j-1;  val[7] = - CC;
                col[8].i = i+1; col[8].j = j+1;  val[8] = - CC;
                ierr = MatSetValuesStencil(P,1,&row,9,col,val,INSERT_VALUES); CHKERRQ(ierr);
            }
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


int main(int argc,char **argv)
{
    PetscErrorCode ierr;
    AppCtx     		 user;
    TS             ts;
    Vec            x;
    DM             da;
    DMDALocalInfo  info;
    double         noiselevel = 0.15;

    PetscInitialize(&argc,&argv,NULL,"Solve coupled PDE");

    user.L      = 2.5;
    user.Du     = 8.0e-5;
    user.Dv     = 4.0e-5;
    user.phi    = 0.05;
    user.kappa  = 0.063;

    PetscInt Nx = 3;
    PetscInt Ny = 3;
    PetscInt ndof = 2;
    PetscInt stencil_width = 1;

    ierr = DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, DMDA_STENCIL_BOX, Nx,Ny,PETSC_DECIDE,PETSC_DECIDE, ndof, stencil_width, NULL,NULL,&da); CHKERRQ(ierr);
    ierr = DMSetFromOptions(da); CHKERRQ(ierr);
    ierr = DMSetUp(da); CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,0,"u"); CHKERRQ(ierr);
    ierr = DMDASetFieldName(da,1,"v"); CHKERRQ(ierr);
    ierr = DMDAGetLocalInfo(da,&info); CHKERRQ(ierr);

    ierr = DMDASetUniformCoordinates(da, 0.0, user.L, 0.0, user.L, -1.0, -1.0); CHKERRQ(ierr);

    ierr = TSCreate(PETSC_COMM_WORLD,&ts); CHKERRQ(ierr);
    ierr = TSSetDM(ts,da); CHKERRQ(ierr); // Link the time-stepper with the DMDA
    ierr = TSSetApplicationContext(ts,&user); CHKERRQ(ierr);
    ierr = DMDATSSetRHSFunctionLocal(da,INSERT_VALUES, (DMDATSRHSFunctionLocal)FormRHSFunctionLocal,&user); CHKERRQ(ierr);
    ierr = DMDATSSetRHSJacobianLocal(da,(DMDATSRHSJacobianLocal)FormRHSJacobianLocal,&user); CHKERRQ(ierr);
    ierr = DMDATSSetIFunctionLocal(da,INSERT_VALUES,(DMDATSIFunctionLocal)FormIFunctionLocal,&user); CHKERRQ(ierr);
    ierr = DMDATSSetIJacobianLocal(da,(DMDATSIJacobianLocal)FormIJacobianLocal,&user); CHKERRQ(ierr);

    ierr = TSSetType(ts,TSARKIMEX); CHKERRQ(ierr);
    ierr = TSSetTime(ts,0.0); CHKERRQ(ierr);
    ierr = TSSetMaxTime(ts,15000.0); CHKERRQ(ierr);
    ierr = TSSetTimeStep(ts,5.0); CHKERRQ(ierr);
    ierr = TSSetExactFinalTime(ts,TS_EXACTFINALTIME_MATCHSTEP); CHKERRQ(ierr);
    ierr = TSSetFromOptions(ts);CHKERRQ(ierr);

    ierr = DMCreateGlobalVector(da,&x); CHKERRQ(ierr);
    ierr = InitialState(da,x,noiselevel,&user); CHKERRQ(ierr);
    ierr = TSSolve(ts,x); CHKERRQ(ierr);

    VecDestroy(&x);  TSDestroy(&ts);  DMDestroy(&da);
    return PetscFinalize();
}
