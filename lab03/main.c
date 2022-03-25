#include <stdio.h>
#include <stdlib.h>
#include <math.h>

////////////////// parametry stale ////////////////////

#define x0 0.01
#define v0 0.0
#define dt0 1.0
#define S 0.75
#define p 2.0
#define alpha 5.0
#define t_max 40.0

//Zaleznie od cwiczenia
#define TOL 0.01

////////////////// funkcje ////////////////////

typedef struct XV{
  double x;
  double v;
} XV;

double f(double x, double v){ return v; }
double g(double x, double v){ return alpha*(1-pow(x,2))*v - x; }

XV update_rk2(double x, double v, double dt){
  double k1x = v;
  double k1v = g(x, v);
  double k2x = v + dt *k1v;
  double k2v = g(x + dt*k1x, v + dt*k1v);
  x += dt*(k1x+k2x)/2.0;
  v += dt*(k1v+k2v)/2.0;

  XV new_xv;
  new_xv.x = x;
  new_xv.v = v;

  return new_xv;
}

XV update_trapez (double x, double v, double dt){
  double new_x = x;
  double new_v = v;
  double dx, dv;
  double delta = pow(10, -10);
  do{

    double F = new_x - x - dt*(f(x, v) + f(new_x, new_v))/2.0;
    double G = new_v - v - dt*(g(x, v) + g(new_x, new_v))/2.0;
    double a11 = 1.0;
    double a12 = -dt/2.0;
    double a21 = -dt*(-2*alpha*new_x*new_v - 1)/2.0;
    double a22 = 1 - dt*alpha*(1 - pow(new_x, 2))/2.0;

    dx = (-F*a22 - (-G)*a12)/(a11*a22 - a12*a21);
    dv = (a11*(-G) - a21*(-F))/(a11*a22 - a12*a21);
    
    x += dx;
    v += dv;
  } while(fabs(dx)<delta && fabs(dv)<delta);

  XV new_xv;
  new_xv.x = x;
  new_xv.v = v;

  return new_xv;
}

int main(){

//zamiennie w zaleznosci od wywolania
//FILE* f1 = fopen("rk2_tol1.dat", "a+");
FILE* f1 = fopen("trapez_tol1.dat", "a+");
//FILE* f1 = fopen("rk2_tol2.dat", "a+");
//FILE* f1 = fopen("trapez_tol2.dat", "a+");

double t = 0.0;
double dt = dt0;
double x = x0;
double v = v0;

while(t < t_max){
  XV new_11, new_12, new_2;
  
  //new_11 = update_rk2(x, v, dt);
  //new_12 = update_rk2(new_11.x, new_11.v, dt);
  //new_2 = update_rk2(x, v, 2*dt);
  
  new_11 = update_trapez(x, v, dt);
  new_12 = update_trapez(new_11.x, new_11.v, dt);
  new_2 = update_trapez(x, v, 2*dt);
  
  double ex = (new_12.x - new_2.x) / (pow(2,p) - 1);
  double ev = (new_12.v - new_2.v) / (pow(2,p) - 1);
  double e_max = fabs(ex) > fabs(ev) ? fabs(ex) : fabs(ev);
  if (e_max < TOL){
    t += 2*dt;
    x = new_12.x;
    v = new_12.v;
    fprintf(f1, "%f %f %f %f\n", t, dt, x, v);
  }

  dt *= pow(S*TOL/e_max, (1.0/(p + 1.0)));

}

fclose(f1);


return 0;
}
