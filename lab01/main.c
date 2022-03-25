#include <stdio.h>
#include <stdlib.h>
#include<math.h>
#define R 100
#define L 0.1
#define C 0.001
#define DT 0.0001
#define PI 3.14159265
#define om0 1.0/sqrt(0.0001)
//Wartosc ta byla zmieniana przy kolejnych wykonaniach programu, dla drugiej czesci zadania
#define OMV 1.2*om0


double an_val(double lambda, double t){
  return pow(exp(1), lambda*t);
}

void euler(double* arr, int n, double dt, int lambda){
  for(int i = 1; i < n ; i ++){
    arr[i] = arr[i-1] + arr[i-1]*dt*lambda;
  }
}

void rk2(double* arr, int n, double dt, int lambda){
  double k1, k2;
  for(int i = 1; i < n ; i ++){
    k1 = lambda*arr[i-1];
    k2 = lambda*(arr[i-1] + dt*k1);
    arr[i] = arr[i-1] + dt*(k1+k2)/2.0;
  }
}

void rk4(double* arr, int n, double dt, int lambda){
  double k1, k2, k3, k4;
  for(int i = 1; i < n ; i ++){
    k1 = lambda*arr[i-1];
    k2 = lambda*(arr[i-1] + dt*k1/2.0);
    k3 = lambda*(arr[i-1] + dt*k2/2.0);
    k4 = lambda*(arr[i-1] + dt*k3);
    arr[i] = arr[i-1] + dt*(k1+2*k2+2*k3+k4)/6.0;
  }
}

void err(double* arr, int n, double dt){
  for(int i = 0 ; i < n ; i ++){
    arr[i] = fabs(arr[i] - an_val(-1, dt*i));
  }
}

void create_file(double* arr, int n, double dt, char* name, int _i){
  FILE* file = fopen(name, "a+");
  for(int i = _i ; i < n ; i ++){
    fprintf(file, "%f %.13f\n", i*dt, arr[i]);
  }
  fclose(file);
}

////////////////// 2 func ////////////////////

double v(double omv, double t){
  return 10*sin(omv*t);
}

double f(double t, double Q, double I){
  return I;
}

double g(double t, double Q, double I){
  return (v(OMV, t)/L - Q/(L*C) - R*I/L);
}

int main(){

double y0 = 1;
double lambda = -1;
double y_a[501];
double y_b[51];
double y_c[6];
double y_an[501];
double t = 0.0;

double d_a[501];
double d_b[51];
double d_c[6];
///////////////////// 1.1 ////////////////////

//Ustawienie wartosci poczatkowych
y_a[0]=y_b[0]=y_c[0] = y0;
y_an[0] = an_val(0,0);

//Wypelnienie tablicy wartosci an alitycznych
for(int i = 1; i < 501 ; i ++){
  y_an[i] = an_val(lambda, i*0.01);
}
create_file(y_an, 501, 0.01, "an.dat", 0);

//Wypelnienie tablic wartosciami

euler(y_a, 501, 0.01, lambda);
euler(y_b, 51, 0.1, lambda);
euler(y_c, 6, 1, lambda);

create_file(y_a, 501, 0.01, "1_a.dat", 0);
create_file(y_b, 51, 0.1, "1_b.dat", 0);
create_file(y_c, 6, 1, "1_c.dat", 0);

err(y_a, 501, 0.01);
err(y_b, 51, 0.1);
err(y_c, 6, 1);

create_file(y_a, 501, 0.01, "err1_a.dat", 1);
create_file(y_b, 51, 0.1, "err1_b.dat", 1);
create_file(y_c, 6, 1, "err1_c.dat", 1);

//////////////////// 1.2 /////////////////////

rk2(y_a, 501, 0.01, lambda);
rk2(y_b, 51, 0.1, lambda);
rk2(y_c, 6, 1, lambda);

create_file(y_a, 501, 0.01, "2_a.dat", 0);
create_file(y_b, 51, 0.1, "2_b.dat", 0);
create_file(y_c, 6, 1, "2_c.dat", 0);

err(y_a, 501, 0.01);
err(y_b, 51, 0.1);
err(y_c, 6, 1);

create_file(y_a, 501, 0.01, "err2_a.dat", 1);
create_file(y_b, 51, 0.1, "err2_b.dat", 1);
create_file(y_c, 6, 1, "err2_c.dat", 1);

/////////////////// 1.3 /////////////////////

rk4(y_a, 501, 0.01, lambda);
rk4(y_b, 51, 0.1, lambda);
rk4(y_c, 6, 1, lambda);

create_file(y_a, 501, 0.01, "3_a.dat",0);
create_file(y_b, 51, 0.1, "3_b.dat",0);
create_file(y_c, 6, 1, "3_c.dat",0);

err(y_a, 501, 0.01);
err(y_b, 51, 0.1);
err(y_c, 6, 1);

create_file(y_a, 501, 0.01, "err3_a.dat", 1);
create_file(y_b, 51, 0.1, "err3_b.dat", 1);
create_file(y_c, 6, 1, "err3_c.dat", 1);

////////////////// 2 ////////////////////

double T0 = 2*3.14*sqrt(0.1*0.001);
int n = 4*T0/DT + 1;
double I;
double Q;
double k1q, k1i, k2q, k2i, k3q, k3i, k4q, k4i;
I = Q = 0.0;

FILE* file = fopen("2d.dat", "a+");


for(int i = 0 ; i < n ; i ++){
  double t = i*DT;
  k1q = f(t, Q, I);
  k1i = g(t, Q, I);
  
  k2q = f(t+DT/2.0, Q+DT*k1q/2.0, I+DT*k1i/2.0);
  k2i = g(t+DT/2.0, Q+DT*k1q/2.0, I+DT*k1i/2.0);

  k3q = f(t+DT/2.0, Q+DT*k2q/2.0, I+DT*k2i/2.0);
  k3i = g(t+DT/2.0, Q+DT*k2q/2.0, I+DT*k2i/2.0);

  k4q = f(t+DT/2.0, Q+DT*k3q, I+DT*k3i);
  k4i = g(t+DT/2.0, Q+DT*k3q, I+DT*k3i);
 
  fprintf(file, "%f %f %f\n", DT*i, Q, I);

  Q = Q + DT*(k1q + 2*k2q + 2*k3q + k4q)/6.0;
  I = I + DT*(k1i + 2*k2i + 2*k3i + k4i)/6.0;
}

fclose(file);

return 0;
}
