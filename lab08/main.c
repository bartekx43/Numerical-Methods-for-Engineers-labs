#include <stdio.h>
#include <stdlib.h>
#include <math.h>

////////////////// DEFINICJE ////////////////////

#define nx 400
#define ny 90
#define i1 200
#define i2 210
#define j1 50
#define delta 0.01
#define sigma 10*delta
#define xA 0.45
#define yA 0.45
#define IT_MAX 9500

////////////////// FUNKCJE ////////////////////

void read_psi(float psi[nx+1][ny+1]){
  int a, b;
  float tmp;
  FILE* file = fopen("psi.dat", "r");
  for(int i = 0 ; i <= nx ; i ++){
    for(int j = 0 ; j <= ny ; j ++){
      fscanf(file, "%d", &a);
      fscanf(file, "%d", &b);
      fscanf(file, "%f", &tmp);
      psi[i][j]=tmp;
    }
  }
  fclose(file);
}

void get_V(float Vx[nx+1][ny+1], float Vy[nx+1][ny+1], float psi[nx+1][ny+1]){
  for(int i = 1 ; i <= nx-1 ; i ++){
    for(int j = 1 ; j <= ny-1 ; j ++){
      Vx[i][j] = (psi[i][j+1]-psi[i][j-1])/(2*delta);
      Vy[i][j] = -(psi[i+1][j]-psi[i-1][j])/(2*delta);
    }
  }

  for(int i = i1 ; i <= i2 ; i ++){
    for(int j = 0 ; j <= j1 ; j ++){
      Vx[i][j] = Vy[i][j] = 0.0;
    }
  }

  for(int i = 1 ; i <= nx-1  ; i ++){
    Vx[i][0] = Vy[i][ny] = 0.0;
  }

  for(int i = 0 ; i <= ny ; i ++){
    Vx[0][i] = Vx[1][i];
    Vx[nx][i] = Vx[nx-1][i];
  }
}

float get_max(float Vx[nx+1][ny+1], float Vy[nx+1][ny+1]){
  float v_max = 0.0;
  for(int i = 1 ; i <= nx ; i ++){
    for(int j = 0 ; j <= ny ; j ++){
      float tmp = sqrt(pow(Vx[i][j],2)+pow(Vy[i][j],2));
      if (tmp > v_max) v_max = tmp;
    }    
  }
  return v_max;
}

void initialize_U(float U[nx+1][ny+1]){
  for(int i = 0 ; i <= nx ; i ++){
    for(int j = 0 ; j <= ny ; j ++){
      U[i][j] = (1.0/(2*M_PI*sigma*sigma))*exp(-(pow(delta*i - xA,2) + pow(delta*j - yA,2))/(2*sigma*sigma));
    }
  }
}

void write_to_file(float arr[nx+1][ny+1], char* filename){
  FILE* f1 = fopen(filename, "w+");
  for(int i = 0 ; i <= nx ; i ++){
    for(int j = 0 ; j <= ny ; j ++){
      fprintf(f1, "%d %d %f\n", i, j, arr[i][j]);
    }
    fprintf(f1, "\n");
  }
  fclose(f1);
}

void new_U(float U[nx+1][ny+1], float Un[nx+1][ny+1], float Vx[nx+1][ny+1], float Vy[nx+1][ny+1], float dt, float D){
  for(int i = 0 ; i <= nx ; i ++){
    for(int j = 1 ; j <= ny-1 ; j ++){
      if (i >= i1 && i <= i2 && j <= j1){}
      else if(i==0){
        U[i][j] = (1.0/(1.0+2*D*dt/delta/delta)) * (Un[i][j]-dt/2.0*Vx[i][j]*((Un[i+1][j]-Un[nx][j])/(2*delta)+(U[i+1][j]-U[nx][j])/(2*delta)) - dt/2.0*Vy[i][j]*((Un[i][j+1]-Un[i][j-1])/(2*delta)+(U[i][j+1]-U[i][j-1])/(2*delta)) + dt/2.0*D*((Un[i+1][j]+Un[nx][j]+Un[i][j+1]+Un[i][j-1]-4*Un[i][j])/pow(delta,2)+(U[i+1][j]+U[nx][j]+U[i][j+1]+U[i][j-1])/pow(delta,2)));
      }
      else if(i==nx){
        U[i][j] = (1.0/(1.0+2*D*dt/delta/delta)) * (Un[i][j]-dt/2.0*Vx[i][j]*((Un[0][j]-Un[i-1][j])/(2*delta)+(U[0][j]-U[i-1][j])/(2*delta)) - dt/2.0*Vy[i][j]*((Un[i][j+1]-Un[i][j-1])/(2*delta)+(U[i][j+1]-U[i][j-1])/(2*delta)) + dt/2.0*D*((Un[0][j]+Un[i-1][j]+Un[i][j+1]+Un[i][j-1]-4*Un[i][j])/pow(delta,2)+(U[0][j]+U[i-1][j]+U[i][j+1]+U[i][j-1])/pow(delta,2)));
      }
      else{
        U[i][j] = (1.0/(1.0+2*D*dt/delta/delta)) * (Un[i][j]-dt/2.0*Vx[i][j]*((Un[i+1][j]-Un[i-1][j])/(2*delta)+(U[i+1][j]-U[i-1][j])/(2*delta)) - dt/2.0*Vy[i][j]*((Un[i][j+1]-Un[i][j-1])/(2*delta)+(U[i][j+1]-U[i][j-1])/(2*delta)) + dt/2.0*D*((Un[i+1][j]+Un[i-1][j]+Un[i][j+1]+Un[i][j-1]-4*Un[i][j])/pow(delta,2)+(U[i+1][j]+U[i-1][j]+U[i][j+1]+U[i][j-1])/pow(delta,2)));
      }
    }
  }    
}


////////////////// MAIN ////////////////////

int main(){
  
  float psi[nx+1][ny+1];
  read_psi(psi);

  float Vx[nx+1][ny+1] = {{0.0}};
  float Vy[nx+1][ny+1] = {{0.0}};
  get_V(Vx, Vy, psi);
  
  float v_max = get_max(Vx, Vy);
  float dt = delta/(4.0*v_max);

  float U[nx+1][ny+1] = {{0.0}};
  initialize_U(U);
  float Un[nx+1][ny+1];
  initialize_U(Un);

  write_to_file(Vx, "Vx.dat");
  write_to_file(Vy, "Vy.dat");

  FILE* file = fopen("dane_a.dat", "w+");
  for(int it = 0 ; it <= IT_MAX ; it ++){  
    
    for(int i = 0 ; i <= nx ; i ++){
      for(int j = 0 ; j <= ny ; j ++){
        U[i][j] = Un[i][j];
      }
    }

    for(int i = 0 ; i < 20 ; i ++){
      new_U(U, Un, Vx, Vy, dt, 0.0);
    }

    double c = 0.0;
    double x_sr = 0.0;
    for(int i = 0 ; i <= nx ; i ++){
      for(int j = 0 ; j <= ny ; j ++){
        Un[i][j] = U[i][j];
        c += U[i][j]*delta*delta;
        x_sr += delta*i*U[i][j]*delta*delta;
      }
    }
    fprintf(file, "%f %f %f\n", it*delta, c, x_sr);
    if (it == 0){
      write_to_file(U, "U0_a.dat");
    }
    if(it == IT_MAX/5.0){
      write_to_file(U, "U1_a.dat");
    }
    if(it == 2*IT_MAX/5.0){
      write_to_file(U, "U2_a.dat");
    }
    if(it == 3*IT_MAX/5.0){
      write_to_file(U, "U3_a.dat");
    }
    if(it == 4*IT_MAX/5.0){
      write_to_file(U, "U4_a.dat");
    }
    if(it == IT_MAX-1){
      write_to_file(U, "U5_a.dat");
    }
    if(it%100 == 0){
      printf("%d %f\n", it, c);
    }
  }
  fclose(file);

  return 0;
}




















































































