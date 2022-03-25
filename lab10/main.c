#include <stdio.h>
#include <stdlib.h>
#include <math.h>

////////////////// DEFINICJE ////////////////////

#define nx 150
#define nt 1000
#define delta 0.1
#define dt 0.05
#define xa 7.5
#define sigma 0.5

////////////////// FUNKCJE ////////////////////

void initial(double U[nx+1], double V[nx+1]){
    U[0] = U[nx] = 0.0;
    V[0] = V[nx] = 0.0;
    for(int i = 1 ; i < nx ; i ++){
        U[i] = exp(-(pow(delta*i-xa, 2)/(2*pow(sigma,2))));
        V[i] = 0.0;
    }
}

void copy(double U0[nx+1], double U[nx+1]){
    for(int i = 1; i < nx ; i++){
        U0[i] = U[i];
    }
}

void fillA(double A[nx+1], double U[nx+1], double U0[nx+1], double alpha, double beta, int it){
    for(int i = 1; i < nx ; i++){
        ///////////////////////////////////////////////// TU TRZEBA DOPISAĆ COŚ PO ALFIE!!!!!!!!!!
        A[i] = (U[i+1] - 2*U[i] + U[i-1])/pow(sigma,2) - beta*(U[i] - U0[i])/dt + alpha*cos(50*it/nt);
    }
}

void iter(double U0[nx+1], double U[nx+1], double V[nx+1], double Vp[nx+1], double A[nx+1], double alpha, double beta, int it){
    for(int i = 1 ; i < nx ; i ++) Vp[i] = V[i] + dt*A[i]/2.0;
    for(int i = 1 ; i < nx ; i ++) copy(U0, U);
    for(int i = 1 ; i < nx ; i ++) U[i] = U[i] + dt*Vp[i];
    for(int i = 1 ; i < nx ; i ++) fillA(A, U, U0, alpha, beta, it);
    for(int i = 1 ; i < nx ; i ++) V[i] = Vp[i] + dt*A[i]/2.0;
}

////////////////// MAIN ////////////////////

int main(){

    double U0[nx+1], U[nx+1], V[nx+1], Vp[nx+1], A[nx+1];
    double alpha = 0.0;
    double beta = 0.0;
    FILE* f1 = fopen("E0.dat", "w+");
    FILE* f2 = fopen("U0.dat", "w+");

    initial(U, V);
    copy(U0, U);
    fillA(A, U, U0, alpha, beta, 1);
    for(int it = 1 ; it <= nt ; it++){
        iter(U0, U, V, Vp, A, alpha, beta, it);
        double E = delta/4.0*(pow((U[1]-U[0])/delta,2)+pow((U[nx]-U[nx-1])/delta,2));
        for(int i = 1 ; i < nx ; i ++){
            E += delta/2.0*(pow(V[i],2)+pow((U[i+1]-U[i-1])/(2*delta),2));
        }
        fprintf(f1, "%f %f\n", it*dt, E);
        for(int i = 0; i <= nx ; i ++){
            fprintf(f2, "%f %f %f\n", it*dt, i*delta, U[i]);
        }
        fprintf(f2, "\n");
    }

    fclose(f1);
    fclose(f2);

    return 0;
}