#include <stdio.h>
#include <stdlib.h>
#include <math.h>
////////////////// definicje ////////////////////
#define eps 1
#define delta 0.1
#define nx 150
#define ny 100
#define V1 10.0
#define V2 0.0
#define x_max delta*nx
#define y_max delta*ny
#define sx 0.1*x_max
#define sy 0.1*y_max
#define TOL 1e-8
////////////////// zmienny parametr ////////////////////
#define omegaG 1.0
#define omegaL 1.9
////////////////// funkcje ////////////////////
double get_rho1(double x, double y){
  return exp(-pow(x-0.35*x_max,2)/pow(sx,2) - pow(y-0.5*y_max,2)/pow(sy,2));
}

double get_rho2(double x, double y){
  return -exp(-pow(x-0.65*x_max,2)/pow(sx,2) - pow(y-0.5*y_max,2)/pow(sy,2));
}

double get_rho(double x, double y){ 
  return get_rho1(x,y)+get_rho2(x,y);
}

void fill_WB(double V[nx+1][ny+1]){
  for(int i = 0 ; i < nx+1 ; i ++){
    V[i][0] = V1;
    V[i][ny] = V2;
  }
}

void global_step(double Vn[nx+1][ny+1], double Vs[nx+1][ny+1], double rho[nx+1][ny+1]){

  for(int i = 1 ; i <= nx-1 ; i ++){
    for(int j = 1 ; j <= ny-1 ; j ++){
      Vn[i][j] = (Vs[i+1][j]+Vs[i-1][j]+Vs[i][j+1]+Vs[i][j-1]+pow(delta,2)*rho[i][j]/eps)/4.0;
    }
  }
  for(int i = 1 ; i <= ny-1 ; i ++){
    Vn[0][i] = Vn[1][i];
    Vn[nx][i] = Vn[nx-1][i];
  }
  for(int i = 0 ; i < nx+1 ; i ++){
    for(int j = 0 ; j < ny+1 ; j ++){
      Vs[i][j] = (1-omegaG)*Vs[i][j] + omegaG*Vn[i][j];
    }
  }
}

void local_step(double V[nx+1][ny+1], double rho[nx+1][ny+1]){
 
  for(int i = 1 ; i <= nx-1 ; i ++){
    for(int j = 1 ; j <= ny-1 ; j ++){
      V[i][j] = (1-omegaL)*V[i][j]+omegaL/4.0*(V[i+1][j]+V[i-1][j]+V[i][j+1]+V[i][j-1]+pow(delta,2)*rho[i][j]/eps);
    }
  }
  for(int i = 1 ; i <= ny-1 ; i ++){
    V[0][i] = V[1][i];
    V[nx][i] = V[nx-1][i];
  }
}

double get_d2V(double V[nx+1][ny+1], int i, int j){
  return (V[i+1][j]-2*V[i][j]+V[i-1][j]+V[i][j+1]-2*V[i][j]+V[i][j-1])/pow(delta,2);
}

double get_s(double V[nx+1][ny+1], double rho[nx+1][ny+1]){
  double s = 0.0;
  for(int i = 0 ; i < nx ; i ++){
    for(int j = 0 ; j < ny ; j ++){
      s += pow((V[i+1][j]-V[i][j])/delta,2)/2.0 + pow((V[i][j+1]-V[i][j])/delta,2)/2.0 + rho[i][j]*V[i][j];
    }
  }
  s *= pow(delta,2);
  return s;
}

////////////////// main ////////////////////
int main(){

double V[nx+1][ny+1] = {0};
fill_WB(V);
double Vn[nx+1][ny+1] = {0};
fill_WB(Vn);

double rho[nx+1][ny+1] = {0};
for(int i = 0 ; i < nx+1 ; i ++){
  for(int j = 0 ; j < ny+1 ; j ++){
    rho[i][j] = get_rho(delta*i, delta*j);
  }
}

FILE* file = fopen("sG1.dat", "a+");
//FILE* fg = fopen("global.dat", "a+");
//FILE* ferr = fopen("err.dat", "a+");

double s_old, s_new;
s_old = s_new = 0.0;
do{
  global_step(Vn, V, rho);
  s_old = s_new;
  s_new = get_s(V, rho);
  fprintf(file, "%f\n", s_new);
} while(fabs((s_new - s_old)/s_old) >= TOL);
/*
for(int i = 0 ; i < nx+1 ; i ++){
  for(int j = 0 ; j < ny+1 ; j ++){
    fprintf(fg, "%f %f %f\n", i*delta, j*delta, V[i][j]);  
    fprintf(ferr, "%f %f %f\n", i*delta, j*delta, get_d2V(V, i, j)+rho[i][j]/eps);  
  }
  fprintf(fg, "\n");
  fprintf(ferr, "\n");
}

fclose(fg);
fclose(ferr);
*/
fclose(file);

return 0;
}
