#include <stdio.h>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <string.h>
////////////////// definicje ////////////////////
#define delta 0.2
#define nx 128
#define ny 128
#define x_max delta*nx
#define y_max delta*ny
#define TOL 1e-8
#define PI 3.14159265
////////////////// funkcje ////////////////////

void fill_WB(double V[nx+1][ny+1]){
  for(int i = 0 ; i < nx+1 ; i ++){
    V[0][i] = sin(1.0*PI*i/nx);
    V[i][ny] = -sin(2.0*PI*i/nx);
    V[nx][i] = sin(1.0*PI*i/nx);
    V[i][0] = sin(2.0*PI*i/nx);
  }
}

void iter(double V[nx+1][ny+1], int k){ 
  for(int i = k ; i < nx ; i += k){
    for(int j = k ; j < ny ; j += k){
      V[i][j] = (V[i+k][j]+V[i-k][j]+V[i][j+k]+V[i][j-k])/4.0;
    }
  }
}

double get_s(double V[nx+1][ny+1], int k){
  double s = 0.0;
  for(int i = 0 ; i < nx ; i += k){
    for(int j = 0 ; j < ny ; j += k){
      s += (pow((V[i+k][j]-V[i][j]+V[i+k][j+k]-V[i][j+k])/(2*k*delta),2)+pow((V[i][j+k]-V[i][j]+V[i+k][j+k]-V[i+k][j])/(2*k*delta),2)); 
    }
  }
  s *= (pow(k*delta,2)/2.0);
  return s;
}

void relax(double V[nx+1][ny+1], int k){

  char f1[10] = {0};
  char f2[10] = {0};
  sprintf(f1, "%d", k);
  sprintf(f2, "%d", k);
  FILE* file = fopen(strcat(f1, "k.dat"), "a+");
  FILE* fg = fopen(strcat(f2, "map.dat"), "a+");
  double s_old, s_new;
  s_old = s_new = 0.0;
  do{
    iter(V, k);
    s_old = s_new;
    s_new = get_s(V, k);
    fprintf(file, "%.30f\n", s_new);
  } while(fabs((s_new - s_old)/s_old) >= TOL);

  for(int i = 0 ; i <= nx ; i += k){
    for(int j = 0 ; j <= ny ; j += k){
      fprintf(fg, "%.1f %.1f %.20f\n", i*delta, j*delta, V[i][j]);  
    }
    fprintf(fg, "\n");
  }
  fclose(fg);
  fclose(file);
}

void narrow(double V[nx+1][ny+1], int k){
  for(int i = 0 ; i < nx ; i += k){
    for(int j = 0 ; j < ny ; j += k){
      V[i+k/2][j+k/2] = (V[i][j]+V[i+k][j]+V[i][j+k]+V[i+k][j+k])/4.0;
      V[i+k][j+k/2] = (V[i+k][j]+V[i+k][j+k])/2.0;
      V[i+k/2][j+k] = (V[i][j+k]+V[i+k][j+k])/2.0;
      V[i+k/2][j] = (V[i][j]+V[i+k][j])/2.0;
      V[i][j+k/2] = (V[i][j]+V[i][j+k])/2.0;
    }
  }
  fill_WB(V);
}

////////////////// main ////////////////////
int main(){

double V[nx+1][ny+1] = {0};
fill_WB(V);

for(int k = 16 ; k >= 1 ; k /= 2){
  printf("relax%d\n", k);
  relax(V, k);
  printf("relaxed%d\n", k);
  narrow(V, k);
  printf("narrowed%d\n", k);
}

return 0;
}
