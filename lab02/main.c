#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 500
#define t_max 100
#define TOL 0.000001
#define dT 0.1
#define gamma 0.1
#define beta 0.001
#define c1 0.5-sqrt(3)/6.0
#define c2 0.5+sqrt(3)/6.0
#define a11 0.25
#define a22 0.25
#define a12 0.25-sqrt(3)/6.0
#define a21 0.25+sqrt(3)/6.0
#define b1 0.5
#define b2 0.5

double f(double u){
  return (beta*N-gamma)*u - beta*u*u;
}

////////////////// 1 ////////////////////

double picard(double u_n){
  double u_n_new = u_n;
  double u_n_mi = u_n;
  double alpha = beta*N - gamma;
  int mi = 0;

  while(mi <= 20 && fabs(u_n_new-u_n_mi)<TOL){
    u_n_mi = u_n_new;
    u_n_new = u_n + dT*(f(u_n)+f(u_n_mi))/2.0;
    mi ++;
  }
  return u_n_new;
}

double newton(double u_n){
  double u_n_new = u_n;
  double u_n_mi = u_n;
  double alpha = beta*N - gamma;
  int mi = 0;

  while(mi <= 20 && fabs(u_n_new-u_n_mi)<TOL){
    u_n_mi = u_n_new;
    u_n_new = u_n_mi - (u_n_mi-u_n-dT*(f(u_n)+f(u_n_mi))/2.0)/(1-dT*(alpha-2*beta*u_n_mi)/2.0);
    mi ++;
  }
  return u_n_new;
}

////////////////// 2 ////////////////////

double F1(double U1, double U2, double u){
  double alpha = beta*N - gamma;
  return U1 - u - dT*(a11*(alpha*U1-beta*pow(U1,2))+a12*(alpha*U2-beta*pow(U2,2)));
}

double F2(double U1, double U2, double u){
  double alpha = beta*N - gamma;
  return U2 - u - dT*(a21*(alpha*U1-beta*pow(U1,2))+a22*(alpha*U2-beta*pow(U2,2)));
}

double next_u(double u){
  double alpha = beta*N - gamma;
  double u_n;
  double dU1, dU2;
  dU1 = dU2 = 0.0;
  double U1, U2;
  u_n = U1 = U2 = u;
  double mi = 0;
  while(mi <= 20 && fabs(dU1)<TOL && fabs(dU2)<TOL){    
    double m11 = 1 - dT*a11*(alpha-2*beta*U1); 
    double m12 = -dT*a12*(alpha-2*beta*U2);
    double m21 = -dT*a21*(alpha-2*beta*U1);
    double m22 = 1 - dT*a22*(alpha-2*beta*U2);
    
    dU1 = (F2(U1,U2,u_n)*m12 - F1(U1,U2,u_n)*m22)/(m11*m22-m12*m21);
    dU2 = (F1(U1,U2,u_n)*m21 - F2(U1,U2,u_n)*m11)/(m11*m22-m12*m21);
    U1 += dU1;
    U2 += dU2;
    u_n = u_n + dT*(b1*f(U1)+b2*f(U2));
    mi ++;
  }
  return u_n;
}

////////////////// main ////////////////////
int main(){

  double u0 = 1;

  ////////////////// 1 ////////////////////

  FILE* file1 = fopen("picard.dat", "a+");
  FILE* file2 = fopen("newton.dat", "a+");
  double u_picard = u0;
  double u_newton = u0;
  double t = 0.0;
  while(t <= 100.0){
    fprintf(file1, "%f %f\n", t, u_picard);
    fprintf(file2, "%f %f\n", t, u_newton);
    u_picard = picard(u_picard);
    u_newton = newton(u_newton);
    t += dT;
  }
  fclose(file1);
  fclose(file2);
  
  ////////////////// 2 ////////////////////

  t = 0.0;
  FILE* file = fopen("rk2.dat", "a+");
  double u = u0;
  while(t <= 100.0){
    fprintf(file,"%f %f\n", t, u);
    u = next_u(u);
    t += dT;
  }
  fclose(file);

  return 0;
}
