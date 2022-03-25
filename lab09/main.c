#include <stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permute_vector.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

#define delta 1.0
#define dt 1.0
#define TA 40.0
#define TB 0.0
#define TC 30.0
#define TD 0.0
#define KB 0.1
#define KD 0.6
#define nx 40
#define ny 40
#define N (nx*ny+nx+ny+1)
#define IT_MAX 2000

void write_to_file(gsl_vector* v, char* filename){
    FILE* f1 = fopen(filename, "w+");
    for(int i = 0 ; i <= nx ; i ++){
        for(int j = 0 ; j <= ny ; j ++){
            int l = i + j*(nx+1);
            fprintf(f1, "%d %d %f\n", i, j, gsl_vector_get(v, l));
        }
    fprintf(f1, "\n");
    }
    fclose(f1);
}

int main(){
    gsl_matrix* A = gsl_matrix_alloc(N,N);
    gsl_matrix* B = gsl_matrix_alloc(N,N);
    gsl_vector* c = gsl_vector_alloc(N);
    for (int i = 0 ; i <= nx ; i++){
        for (int j = 0 ; j <= ny ; j++){
            gsl_matrix_set(A, i, j, 0.0);
            gsl_matrix_set(B, i, j, 0.0);
        }
        gsl_vector_set(c, i, 0.0);
    }
    for (int i = 1 ; i <= nx-1 ; i++){
        for (int j = 1 ; j <= ny-1 ; j++){
            int l = i + j*(nx+1);
            gsl_matrix_set(A, l, l-nx-1, dt/(2.0*delta*delta));
            gsl_matrix_set(A, l, l-1, dt/(2.0*delta*delta));
            gsl_matrix_set(A, l, l+1, dt/(2.0*delta*delta));
            gsl_matrix_set(A, l, l+nx+1, dt/(2.0*delta*delta));
            gsl_matrix_set(A, l, l, -2.0*dt/(delta*delta)-1.0);
            gsl_matrix_set(B, l, l-nx-1, -dt/(2.0*delta*delta));
            gsl_matrix_set(B, l, l-1, -dt/(2.0*delta*delta));
            gsl_matrix_set(B, l, l+1, -dt/(2.0*delta*delta));
            gsl_matrix_set(B, l, l+nx+1, -dt/(2.0*delta*delta));
            gsl_matrix_set(B, l, l, 2.0*dt/(delta*delta)-1.0);
        }
    }
    for(int i = 0 ; i <= nx ; i += nx){
        for (int j = 0 ; j <= ny ; j++){
            int l = i + j*(nx+1);
            gsl_matrix_set(A, l, l, 1.0);
            gsl_matrix_set(B, l, l, 1.0);
            gsl_vector_set(c, l, 0.0);
        }
    }
     for(int i = 1 ; i <= nx-1 ; i ++){
        int l = i + ny*(nx+1);
        gsl_matrix_set(A, l, l-nx-1, -1.0/(KB*delta));
        gsl_matrix_set(A, l, l, 1.0+1.0/(KB*delta));
        gsl_vector_set(c, l, TB);
        for(int k = 0 ; k < N; k++) gsl_matrix_set(B, l, k, 0.0);
        
        l = i + 0*(nx+1);
        gsl_matrix_set(A, l, l+nx+1, -1.0/(KD*delta));
        gsl_matrix_set(A, l, l, 1.0+1.0/(KD*delta));
        gsl_vector_set(c, l, TD);
        for(int k = 0 ; k < N; k++) gsl_matrix_set(B, l, k, 0.0);
    }

    gsl_vector* T = gsl_vector_alloc(N);
    for (int i = 0 ; i <= nx ; i++){
        for (int j = 0 ; j <= ny ; j++){
            if (i == 0) gsl_vector_set(T, i+j*(nx+1), TA);
            else if (i == nx) gsl_vector_set(T, i+j*(nx+1), TC);
            else gsl_vector_set(T, i+j*(nx+1), 0.0);
        }
    }
    int s = 0;
    gsl_permutation* p = gsl_permutation_alloc(N);
    gsl_linalg_LU_decomp(A, p, &s);
    gsl_vector* d = gsl_vector_alloc(N);
    gsl_vector* T_old = gsl_vector_alloc(N);
    for(int it = 0 ; it <= IT_MAX ; it++){
        gsl_blas_dgemv(CblasNoTrans, 1.0, B, T, 0.0, d);
        gsl_blas_daxpy(1.0, c, d);
        if (it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000) {
            for(int i = 0 ; i < N ; i ++) gsl_vector_set(T_old, i, gsl_vector_get(T,i));
            gsl_linalg_LU_solve(A, p, d, T);
            char it_str[5];
            char filename[10] = {0};
            filename[0] = 'T';
            sprintf(it_str, "%d", it);
            strcat(filename, it_str);
            strcat(filename, ".dat");
            write_to_file(T, filename);
            filename[0] = 'E';
            FILE* f = fopen(filename, "w+");
            for(int i = 0 ; i <= nx ; i ++){
                for(int j = 0 ; j <= ny ; j ++){
                    int l = i + j*(nx+1);
                    fprintf(f, "%d %d %.10f\n", i, j, (gsl_vector_get(T, l)-gsl_vector_get(T_old, l))/dt);
                }
                fprintf(f, "\n");
            }
            fclose(f);
        }
        else gsl_linalg_LU_solve(A, p, d, T);
    }
    return 0;
}