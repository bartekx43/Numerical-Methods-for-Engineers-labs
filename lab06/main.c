#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mgmres.c"

////////////////// DEFINITIONS ////////////////////
#define itr_max 500
#define mr 500
#define tol_abs 1e-8
#define tol_rel 1e-8
#define N (nx+1)*(ny+1)
#define delta 0.1
#define x_max delta*nx
#define y_max delta*ny
#define sigma x_max/10.0

#define nx 4
#define ny 4
#define eps1 1.0
#define eps2 1.0
#define v1 10.0
#define v2 -10.0
#define v3 10.0
#define v4 -10.0
////////////////// FUNCTIONS ////////////////////

double rho1(double x, double y){
  return exp(-pow(x-0.25*x_max,2)/pow(sigma,2)-pow(y-0.5*y_max,2)/pow(sigma,2));
}

double rho2(double x, double y){
  return -exp(-pow(x-0.75*x_max,2)/pow(sigma,2)-pow(y-0.5*y_max,2)/pow(sigma,2));
}

int get_j(int l) { return floor(l/(nx+1)); }
int get_i(int l) { return l - get_j(l)*(nx+1); }
	
double eps(int l){
	if(get_i(l)<=nx/2) return eps1;
	else return eps2;
}

void fill(double* a, int* ja, int* ia, double* b, double* V){

  int k = -1;
  for(int l = 0 ; l < N ; l ++){
    int brzeg = 0;
    double vb = 0;

    int j = get_j(l);
    int i = get_i(l);

    if(i==0){
      brzeg = 1;
      vb = v1;
    }
    if(j==ny){
      brzeg = 1;
      vb = v2;
    }
    if(i==nx){
      brzeg = 1;
      vb = v3;
    }
    if(j==0){
      brzeg = 1;
      vb = v4;
    }
    ///////////////////////////// TU TRZEBA ZMIENI RHO /////////////////
    b[l] = 0.0; // -(rho1(delta*i, delta*j) + rho2(delta*i, delta*j));

    if(brzeg==1) b[l] = vb;
    ia[l] = -1;
    //LEWA SKRAJNA PRZEKATNA
    if(l-nx-1>=0 && brzeg==0){
      k++;
      if(ia[l]<0) ia[l]=k;
      a[k] = eps(l)/pow(delta,2); 
      ja[k] = l-nx-1;
    }
    //PODDIAGONALA
    if(l-1>=0 && brzeg==0){
      k++;
      if(ia[l]<0) ia[l]=k;
      a[k] = eps(l)/pow(delta,2);
      ja[k] = l-1;
    }
    //DIAGONALA
    k++;
    if(ia[l]<0) ia[l]=k;
    if(brzeg==0) a[k] = -(2*eps(l)+eps(l+1)+eps(l+nx+1))/pow(delta,2);
    else a[k] = -1;
    ja[k] = l;
    
    //NADDIAGONALA
    if(l<N && brzeg==0){
      k++;
      a[k] = eps(l+1)/pow(delta,2);
      ja[k] = l+1;
    }
    //PRAWA SKRAJNA PRZEKATNA
    if(l<N+nx+1 && brzeg==0){
      k++;
      a[k] = eps(l+nx+1)/pow(delta,2);
      ja[k] = l+nx+1;
    }
  //END FOR
  }
  int nz_num = k+1;
  ia[N] = nz_num;
}

////////////////// MAIN ////////////////////
int main(){
  
  double a[5*N], b[N], V[N];
  int ia[N+1], ja[5*N];
  fill(a,ja,ia,b,V);

  pmgmres_ilu_cr(N, ia[N], ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel);

////////////////// PIERWSZA CZESC ////////////////////
/*
  FILE* file1 = fopen("A.dat", "a+");
  FILE* file2 = fopen("b.dat", "a+");
  for(int l = 0 ; l < ia[N] ; l ++){
    fprintf(file1, "l:%d  i:%d  j:%d  a[l]:%f\n", l, get_i(l), get_j(l), a[l]);
  }
  for(int l = 0 ; l < N ; l ++){
  fprintf(file2, "l:%d  i:%d  j:%d  b[l]:%f\n", l, get_i(l), get_j(l), b[l]);
  }
  fclose(file1);
  fclose(file2);
*/

////////////////// DRUGA CZESC ////////////////////
  FILE* file = fopen("V6a.dat", "a+");
  for(int i = 0 ; i < nx+1 ; i ++){
    for(int j = 0 ; j < ny+1 ; j ++){
      fprintf(file, "%d %d %f\n", j, i, V[j+i*(nx+1)]);
    }
    fprintf(file, "\n");
  }
  fclose(file);

return 0;
}
