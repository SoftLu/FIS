# include <string.h>
# include <math.h>
# include <stdio.h>
# include <stdlib.h>
# include "mmio.h"
# include "definitions.h"  // Definitions of all functions

// All the vectors I'm creating start with entry 0 = 0, then that entry has no meaningful
// information and then the entries correspond to nz entry values (like in matlab)

int main(int argc, char *argv[])
{
    int ret_code;
    MM_typecode matcode;
    FILE * f;
    int M, N, nz, e, u, ind_x, c = 0, i1, i2;         // The matrix has a size of MxN with nz (number of nonzero values)
    int i, *I, *J, k;
    double *val;
    
 
    f = fopen("mysym.mtx", "r");
    
    if(f == NULL)
    {
        printf("Exiting...\n");
        exit(1);
    }

    ret_code = mm_read_mtx_crd_size(f, &M, &N, &nz);
    I = (int *) malloc(nz * sizeof(int));
    J = (int *) malloc(nz * sizeof(int));
    val = (double *) malloc(nz * sizeof(double));
    
    for(i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %lg\n", &I[i], &J[i], &val[i]);
        //I[i]--;  /* adjust from 1-based to 0-based */
        //J[i]--;
    }

    mm_write_banner(stdout, matcode);
    mm_write_mtx_crd_size(stdout, M, N, nz);
    
    for(i = 0; i < nz; i++)
        printf("I %d  J %d \n", I[i], J[i]);
        
    //int num_rows_cols = sizeof(val); // should be equal to nz but it is not...
    //printf("%d\n", N);
    
    //printf("%d\n", N);
    //for(i = 0; i < 10; i++)
      //  fprintf(stdout, "%d ", I[i]);  //There are no information contained in the first entry of I and J, the info starts in I[1] and J[1]
    //printf("%d ",I[90500]);
    // Creo una matriz en la que guardaré el número de veces (de nz values) que hay en cada row y lo guardo en la entry correspondiente a esa row
    //int *num_ind = (int *) calloc((N + 2), sizeof(int)); // Saving the number of times each row is repeated
    int *I_CSC = (int *) malloc((N + 2)* sizeof(int));  // The new vector I in CSR format
    int *J_CSC = (int *) malloc((nz + 1)* sizeof(int));  // The new vector J in CSR format
    double *V_CSC = (double *) malloc((nz + 1)* sizeof(double));  // The new vector V in CSR format
    // int *I_CSR = (int *) malloc((N + 2)* sizeof(int));  // The new vector I in CSR format
    // int *J_CSR = (int *) malloc((nz + 1)* sizeof(int));  // The new vector J in CSR format
    // double *V_CSR = (double *) malloc((nz + 1)* sizeof(double));
    double *x = (double *) malloc((N + 1) * sizeof(double));
    //double y = (double *) malloc((N + 2) * sizeof(int));
    double y[N + 1], V_CSR[nz + 1]; 
    int I_CSR[N + 2], J_CSR[nz + 1];// I_CSC[N + 2], J_CSC[nz + 1];
    //printf("%d\n", N);

    // Todos los componentes del nuevo vector a 0:
    // for (i = 0; i < N + 2; i++) 
    // {
    //     num_ind[i] = 0;
    //     I_CSR[i] = 0;
    // }
    for (i = 0; i < N + 2; i++) 
    {
        x[i] = 1;
        //y[i] = 0;
        //printf("x[%d] = %lf\n", i, x[i]);
    }
    for (i = 0; i < nz + 1; i++) 
    {
        J_CSR[i] = 0;
        J_CSC[i] = 0;
        V_CSR[i] = 0;
        V_CSC[i] = 0;
    }

    // double V_CSC[nz + 1];
    // int I_CSC[N + 2], J_CSC[nz + 1];
    // for (i = 0; i < nz + 1; i++) 
    // {
    //     //J_CSR[i] = 0;
    //     J_CSC[i] = 0;
    //     //V_CSR[i] = 0;
    //     V_CSC[i] = 0;
    // }


    //for(i = 0; i < N + 2; i++)    
        //printf("%d %d\n", num_ind[i], I_CSR[i]);

    // Going through all values (all rows) seeing how many nz there are in each row
    // for (i = 1; i < N + 2; i++) 
    // {
    //     for (e = 0; e < nz; e++) // hasta nz
    //     {
    //         if (I[e] == i)
    //         {
    //             c++;
    //             num_ind [i] = num_ind [i] + 1;
    //             J_CSR[c] = J[e];    // My new J has all values ordered in ascending corresponding row number
    //             V_CSR[c] = val[e];   // V saves all values in the same order
    //         }
    //     }
        
    // }
    // // for(i = 0; i < nz + 1; i++)    //hasta N+1
    // //     printf("J_CSR %lf \n", &J_CSR[i]);

    // c = 0;
    // for (i = 1; i < N + 2; i++)   // Creating the new I for the CSR format
    // {
    //         I_CSR [i] = i + c;
    //         c = c + num_ind[i] - 1;
    // }
    COO_to_CSR(I_CSR, J_CSR, V_CSR, I, J, val, N, nz);
    COO_to_CSC(I_CSC, J_CSC, V_CSC, I, J, val, N, nz);
    // for(i = 0; i < N + 2; i++)    //hasta N+1
    //     printf("num_ind %d I_CSR %d x %d\n", num_ind[i], I_CSR[i], x[i]);
    
    for(i = 1; i < N + 2; i++)    //hasta N+1
        printf("I_CSR %d \n", I_CSR[i]);   

    for(i = 1; i < nz + 1; i++)    //hasta N+1
        printf("J_CSR %d V_CSR %lf \n", J_CSR[i], V_CSR[i]);

    for(i = 1; i < N + 2; i++)    //hasta N+1
        printf("I_CSC %d x: %lf \n", I_CSC[i], x[i]);   

    for(i = 1; i < nz + 1; i++)    //hasta N+1
        printf("J_CSC %d V_CSC %lf \n", J_CSC[i], V_CSC[i]);

    //Trying mat-vec-sym:
    // mat_vec_mul_sym(y, x, I_CSR, J_CSR, V_CSR, I_CSC, J_CSC, V_CSC, N);
    
    // Trying mat-vec: it works
    // double yfin[N + 2];
    // // for(i = 1; i < N + 1; i++)    //hasta N+1
    // //         printf("fin_y %lf x %lf \n", y[i], x[i]); 
    // mat_vec_mul_CSR(y, x, I_CSR, J_CSR, V_CSR, N);

    // for(i = 1; i < N + 1; i++)    //hasta N+1
    //         printf("fin_y %lf x %lf \n", y[i], x[i]); 

    // // Trying vec-vec: it works
    // double num_fin;
    // num_fin = vec_vec_mul(y, x, N);
    // printf("fin vec-vec %lf \n", num_fin);

    // //Trying modulus: it works 
    // double mod;
    // mod = modulus(x, N);
    // printf("fin mod x %lf \n", mod);

    // Implementation of GMRES:
    //GramSchmidt implementation:

    int m, j, p;
    double mod_r, rho, tol = 1e-8, temp_h, temp_h_1;
    // Declaration of vectors and matrices all intialized as zero and these
    // are reinitialized as m changes
    double *x_star = (double*) calloc(N + 1, sizeof(double));
    for (i = 0; i < N + 1; i++){x_star[i] = 1;} // Intializing the input data all with 1
    
    double *r_z_in = (double*) calloc(N + 1, sizeof(double));
    double *x_z = (double*) calloc(N + 1, sizeof(double));
    double *v_for_r = (double*) calloc(N + 1, sizeof(double));
    double *b = (double*) calloc(N + 1, sizeof(double));
    // Calculating b: b = A*x_star
    mat_vec_mul_CSR(b, x_star, I_CSR, J_CSR, V_CSR, N);
    for(i = 0; i < N + 1; i++)    //hasta N+1
            printf("b %lf \n", b[i]);
    
    // Calculating r_z_in: r_z = b - A*x_z
    // First v_for_r is calculated: v_for_r = A*x_z
    mat_vec_mul_CSR(v_for_r, x_z, I_CSR, J_CSR, V_CSR, N);
    // Finally r_z_in is computed
    for(i = 0; i < N + 1; i++) // Equal r_z to b - A*X_z
    {
        r_z_in[i] = b[i] - v_for_r[i];
    }

    // Calculus of rho and tol
    // Calculating modulus of r_z = r
    mod_r = modulus(r_z_in, N);
    //printf("mod %lf \n", mod_r);
    rho = mod_r;
    tol = (modulus(r_z_in, N))/(modulus(r_z_in, N));
    printf("TOL: %.18lf \n",tol);
    

   
    m = 4;  // The value of m should be updated

    while (j < m + 1)//(rho > tol) //(j < m + 1)//(rho > tol) //(j < m + 1)
    {
        // Declaration of vectors and matrices all intialized as zero and these
        // are reinitialized as m changes
        double **H = (double**)calloc(m + 1, sizeof(double));
        for(i = 0; i < m + 2; i++) 
        { 
            H[i] = (double*)calloc(m + 2, sizeof(double));
        }
        // Now I can access H entries as H[i][j]
        double **v_gk = (double**) calloc(N + 1, sizeof(double));
        for(i = 0; i < m + 2; i++) 
        { 
            v_gk[i] = (double*)calloc(m + 2, sizeof(double));
        }
        double *r_z = (double*) calloc(N + 1, sizeof(double));
        double *v_for_r = (double*) calloc(N + 1, sizeof(double));
        double *v_R = (double*) calloc(m + 1, sizeof(double)); // lo que me sale de multiplicar Vm(R^-1*g)
        double *w = (double*) calloc(N + 1, sizeof(double));
        //double *x_z = (double*) calloc(N + 1, sizeof(double));
        double *C = (double*) calloc(m + 2, sizeof(double));
        double *S = (double*) calloc(m + 2, sizeof(double));
        double *g = (double*) calloc(m + 2, sizeof(double));
        double *R_g = (double*) calloc(m + 1, sizeof(double)); // el resultado de R^-1*g
        double *e1 = (double*) calloc(m + 2, sizeof(double));
        e1[1] = 1;
        
        // Calculating r_z: r_z = b - A*x_z
        // First v_for_r is calculated: v_for_r = A*x_z
        mat_vec_mul_CSR(v_for_r, x_z, I_CSR, J_CSR, V_CSR, N);
        // Finally r_z is computed
        for(i = 0; i < N + 1; i++) // Equal r_z to b - A*X_z
        {
            r_z[i] = b[i] - v_for_r[i];
        }

        // for(i = 0; i < N + 2; i++)    //hasta N+1
        //         printf("r_z %lf \n", r_z[i]);
        // Calculating modulus of r_z = r
        mod_r = modulus(r_z, N);
        printf("RES: %.18lf \n", mod_r);
        rho = mod_r;
        tol = (modulus(r_z, N))/(modulus(r_z_in, N));
        //tol = (modulus(r_z, N))/(modulus(r_z, N));

        // Here the while should start with while rho > tol
        // Calculating g and v1
        for(i = 0; i < N + 2; i++)
        {
            g[i] = mod_r * e1[i];
            v_gk[1][i] = r_z[i] / mod_r;
        }
        // for(i = 0; i < N + 2; i++)
        //     printf("v_gk %lf \n", v_gk[1][i]);
        
        for (j = 1; j < m + 1; j++) // Para que al menos haga una iteración cuando m=1
        {   
            printf("j = %d \n", j);
            printf("m = %d \n", m);
            for(i = 1; i < N + 1; i++){
                printf("v_gk %.18lf \n", v_gk[j][i]);
                w[i] = 0;
            }

            mat_vec_mul_CSR(w, v_gk[j], I_CSR, J_CSR, V_CSR, N);

            for(i = 1; i < N + 1; i++){
                printf("w 2nd %.18lf \n", w[i]);
            }
                
            G_K(v_gk, H, I_CSR, J_CSR, V_CSR, w, m, j, N);

            for (i = 1; i < N + 1; i++){
                printf("GK results V[j]:  %.18lf \n", v_gk[j][i]);
            }
                
            // for (i = 0; i < m + 2; i++)
            // {
            //     for (e = 0; e < m + 1; e++)
            //     {
            //         printf(" %d % d H %lf \n", i, e, H[i][e]);
            //     }
            // }

            // printf("Cj %lf \n", C[j]);
            // printf("Hjj %lf \n", H[j][j]);
            // printf("Hj+1j %lf \n", H[j+1][j]);
            // printf("j %d \n",j);
            // printf("temp_h_1 %lf \n", H[1][j]);
            for (k = 2; k <= j; k++) // Overwritting the values
             {   
                printf("entra k: %d\t j: %d\n", k, j);
                temp_h_1 = H[k - 1][j];
                printf("%.18lf \n", temp_h_1);
                temp_h = H[k][j];
                printf("%.18lf \n", temp_h);
                H[k - 1][j] = (C[k - 1] * temp_h_1) + (S[k - 1] * temp_h);
                H[k][j] = -(S[k - 1] * temp_h_1) + (C[k - 1] * temp_h);
             }
            
            printf(" Givens rotation \n");
            C[j] = (H[j][j]) / (sqrt(pow(H[j][j], 2) + pow(H[j + 1][j], 2)));
            S[j] = (H[j + 1][j]) / (sqrt(pow(H[j][j], 2) + pow(H[j + 1][j], 2)));
            printf("Cj %.18lf \n", C[j]);
            printf("Sj %.18lf \n", S[j]);
            // printf("Hj+1j %lf \n", H[j+1][j]);

            H[j][j] = (C[j] * H[j][j]) + (S[j] * H[j + 1][j]);
            printf("g[j] %.18lf \n", g[j]);
            g[j + 1] = - S[j] * g[j];
            printf("g[j+1] %.18lf \n", g[j+1]);
            g[j] = C[j] * g[j];
            printf("g[j] %.18lf \n", g[j]);

            for (i = 1; i < m + 2; i++)
            {
                for (e = 1; e < m + 1; e++)
                {
                    printf(" %d % d H %.18lf \n", i, e, H[i][e]);
                }
            }

        }
        printf("%d\n", j);
        for(i = 1; i < j; i++){
                printf("g 2nd %.18lf \n", g[i]);}
        inverse_mul(R_g, H, g, m); // para este que la m es fixed
        for (e = 1; e < m + 1; e++)
            printf("R_g[j] %.18lf \n", R_g[e]);
        mat_vec(v_R, v_gk, R_g, m); //tamaño de mi matriz R_g
        for (e = 1; e < m + 1; e++)
            printf("v_R[j] %.18lf \n", v_R[e]);
        for (i = 1; i < N + 1; i++)
        {
            x_z[i] = x_z[i] + v_R[i];
        }
        for (e = 1; e < m + 1; e++)
            printf("x_z[j] %.18lf \n", x_z[e]);
        
        //m = m + 1;
    }

        
    // Printing outputs from GK
    

// Deallocate memory
  // free(num_ind);
  // free(I_CSR);
  // free(J_CSR);
  // free(V_CSR);
  // free(v_gk);
  // free(H);
  // free(x);

	return 0;
}

// FUNCTIONS DEFINITIONS: //

// Conversion from COO to CSR format:
void COO_to_CSR(
    int* I_CSR,
    int* J_CSR,
    double* V_CSR,
    int* I,
    int* J,
    double* val,
    int N,
    int nz
    ){

    int i, e, c = 0 ;
    int *num_ind = (int *) calloc((N + 2), sizeof(int));
    // Going through all values (all rows) seeing how many nz there are in each row
    for (i = 1; i < N + 2; i++) 
    {
        for (e = 0; e < nz + 1; e++) // hasta nz
        {
            if (I[e] == i)
            {
                c++;
                //printf("c: %d\n", c);
                num_ind [i] = num_ind[i] + 1;
                // printf(" i: %d e: %d I: %d \n", i, e, I[e]);
                J_CSR[c] = J[e];    // My new J has all values ordered in ascending corresponding row number
                V_CSR[c] = val[e];   // V saves all values in the same order
            }
        }
        
    }
    // for(i = 0; i < nz + 1; i++)    //hasta N+1
    //     printf("J_CSR %lf \n", &J_CSR[i]);

    c = 0;
    for (i = 1; i < N + 2; i++)   // Creating the new I for the CSR format
    {
            I_CSR [i] = i + c;
            c = c + num_ind[i] - 1;
    }

}

// Conversion from COO to CSC format:
void COO_to_CSC(
    int* I_CSC,
    int* J_CSC,
    double* V_CSC,
    int* I,
    int* J,
    double* val,
    int N,
    int nz
    ){

    int i, e, c = 0;
    int *num_ind = (int *) calloc((N + 2), sizeof(int));
    // Going through all values (all rows) seeing how many nz there are in each row
    for (i = 1; i < N + 2; i++) 
    {
        for (e = 0; e < nz + 1; e++) // hasta nz
        {
            if (J[e] == i)
            {
                c++;
                num_ind [i] = num_ind [i] + 1;
                J_CSC[c] = I[e];    // My new J has all values ordered in ascending corresponding row number
                V_CSC[c] = val[e];   // V saves all values in the same order
                // printf(" i: %d e: %d c: %d V_CSC: %lf \n", i, e, c, V_CSC[e]);
            }
        }
        
    }
    // for(i = 0; i < nz + 1; i++)    //hasta N+1
    //     printf("J_CSR %lf \n", &J_CSR[i]);

    c = 0;
    for (i = 1; i < N + 2; i++)   // Creating the new I for the CSR format
    {
            I_CSC [i] = i + c;
            c = c + num_ind[i] - 1;
    }

}

// Matrix-vector multiplication for symmetric matrices:
void mat_vec_mul_sym(
    double* y,                  // Output: the resulting vector
    const double* const x,      // Input: the vector that is multiplying
    const int* const I_CSR,     // Input: the I_CSR vector of the matrix
    const int* const J_CSR,     // Input: the J_CSR vector of the matrix
    const double* const V_CSR,  // Input: the V_CSR vector of the matrix
    const int* const I_CSC,     // Input: the I_CSR vector of the matrix
    const int* const J_CSC,     // Input: the J_CSR vector of the matrix
    const double* const V_CSC,  // Input: the V_CSR vector of the matrix
    int N 
    ){

    int i, e, ir1, ir2, ic1, ic2;
    // double y[N + 1];

        for (i = 1; i < N + 2; i++)
        {   
            printf("i: %d\n", i);
            ir1 = I_CSR[i];
            ir2 = I_CSR[i + 1] - 1;
            ic1 = I_CSC[i];
            ic2 = I_CSC[i + 1] - 1;
        printf("%d %d \n", ir1, ir2);
        printf("%d %d \n", ic1, ic2);
            //if (ir1 == ir2){
            for (e = ic1; e < ic2 + 1; e++)
            {
                y[i] = y[i] + (V_CSC[e] * x[J_CSC[e]]);
                // printf("V_CSC %lf \n",V_CSC[e]);
                // printf("y: %lf \n", y[i]);
            }
            //}
            //if (ic1 == ic2){
            for (e = ir1; e < ir2; e++)
            {
                y[i] = y[i] + (V_CSR[e] * x[J_CSR[e]]);
                // printf("V_CSR %lf \n",V_CSR[e]);
                // printf("y: %lf \n", y[i]);
            }
            //}
            
        }
         // for(i = 0; i < N + 2; i++)    //hasta N+1
         //     printf("fin_y %lf \n", y[i]); 

}
// Matrix-vector multiplication for CSR matrices:
void mat_vec_mul_CSR(
        double* y,                  // Output: the resulting vector
        const double* const x,      // Input: the vector that is multiplying
        const int* const I_CSR,     // Input: the I_CSR vector of the matrix
        const int* const J_CSR,     // Input: the J_CSR vector of the matrix
        const double* const V_CSR,  // Input: the V_CSR vector of the matrix
        int N                       // Input: the size of the matrix
        ){

        int i, e, i1, i2;
        for (i = 1; i < N + 2; i++)
        {
            i1 = I_CSR[i];
            i2 = I_CSR[i + 1] - 1;
        //printf("%d %d \n", i1, i2);
            for (e = i1; e < i2 + 1; e++)
            {
                y[i] = y[i] + (V_CSR[e] * x[J_CSR[e]]);
            }
        }
        for(i = 0; i < N + 2; i++)    //hasta N+1
            printf("fin_y %lf \n", y[i]); 
    }

// Matrix-vector multiplication:

        
// Vector-vector multiplication:
double vec_vec_mul(
    double* vec1,        // Input: One of the vectors
    double* vec2,        // Input: Other vector
    int N
    ){

    int i;
    double sum_mul = 0;
    //printf("N %d \n", N);
    //N = 3;
    for (i = 1; i < N + 1; i++)
    {
        sum_mul = sum_mul + (vec1[i] * vec2[i]);
        //printf("sum %lf \n",sum);
    }
   
    return sum_mul;
}

// Modulus:
double modulus(
    const double* const vector,
    const int N
    ){

    int i;
    double sum_mod = 0, mod;

    for (i = 1; i < N + 1; i++)
    {
        sum_mod = sum_mod + (pow(vector[i],2));
        //printf("bla");
    }

    mod = sqrt(sum_mod);
    return mod;
}

//Get Krylov method:
void G_K(
    double** v_gk,                   // Output: the resulting vector
    double** h_gk,                  // Output: the resulting new values of h(matrix)
    const int* const I_CSR,         // Input: the I_CSR vector of the matrix
    const int* const J_CSR,         // Input: the J_CSR vector of the matrix
    const double* const V_CSR,      // Input: the V_CSR vector of the matrix              // Input: the first V vector
    double* w,
    int m,                          // Input: the size of the matrix A
    int j,
    int N                           // Input: j number of times the outer for loop will iterate
    ){

    double h;
    int i, e, f, u;
    

    for (i = 1; i < j + 1; i++)
    {   
        // Computation of hij = Vi*w
        h_gk[i][j] = vec_vec_mul(v_gk[i], w, m);
        printf("h_gk %d %d obtenido %.18lf \n", i, j, h_gk[i][j]);

         for(e = 1; e < N + 1; e++)
         {  
            //printf("h_gk obtenido %lf \n", h_gk[i][j]);
            // Updating the value of w inside GK
            w[e] = w[e] - (h_gk[i][j] * v_gk[i][e]);
            printf("%d w %.18lf \n", e, w[e]);
            // printf("h_gk %lf \n", h_gk[i][j]);
            // printf("v_pre_gk %lf \n", v_pre_gk[e]);
         }
    }

    h_gk[j + 1][j] = modulus(w, N);
    printf("h[j+1,j] %.18lf inside GK \n", h_gk[j + 1][j]);

    for(f = 1; f < N + 1; f++)
    {
        v_gk[j + 1][f] = w[f] / h_gk[j + 1][j];
    } 
}

// Inverse of a matrix:
void inverse_mul(
    double* inv_m,
    double **M,
    double* g,
    int m
    ){

    int i, u;

    //inv_m[m] = g[m] / M[m][m];
    printf("%lf\n",g[m]);
    printf("%lf\n",M[m][m]);
    printf("%lf\n",inv_m[m]);

    for (i = m - 1 ; i > - 1; i--)
    {
        inv_m[m] = g[m] / M[m][m];
        //printf("%lf\n",inv_m[m]);
        for (u = m; u > i; u--)
        {
            g[i] = g[i] - (M[i][u]*inv_m[u]);
        }
        inv_m[i] = g[i] / M[i][i];
    }
}

void mat_vec(
    double* vec_fin,
    double** mat,
    double* vec,
    int m                 // debe ir hasta el tamaño de mi matrix
    ){

double sum_bef = 0;
int c, r;

for (r = 1; r < m + 1 ; r++)
    {
        for (c = 1; c < m + 1; c++)
        {
            sum_bef = mat[r][c] * vec[c];
            vec_fin[r] = vec_fin[r] + sum_bef;
        }
    }

}

// Conjugate Gradient function:
void CG(
    double* r_z,      // Input: r = b - A*x
    double* x_z,      // Input: x
    int N,
    ){

    int i, m;
    double alpha, num_alpha, den_alpha, beta, num_beta, den_beta;
    double *p = (double*) calloc(N + 1, sizeof(double));
    // double *alpha = (double*) calloc(N + 1, sizeof(double));
    double *for_den_alpha = (double*) calloc(N + 1, sizeof(double));

    for (i = 0; i < N + 2; i++)
    {
        p[i] = r_z[i];
    }
    for (m = 1; m < N + 2; m++)
    {
        num_alpha = vec_vec_mul(r_z,r_z, N);
        mat_vec_mul_sym(for_den_alpha, p,I_CSR, J_CSR, V_CSR, I_CSC, J_CSC, V_CSC, N);
        den_alpha = vec_vec_mul(for_den_alpha, p, N);
        alpha = (num_alpha) / (den_alpha);
        for (e = 1; e < N + 2; e++)
        {
            x_z[e] = x_z[e] + alpha * p[e];
            r_z[e] = r_z[e] - alpha * p[e];
        }
        num_beta = vec_vec_mul(r_z, r_z, N);
        den_beta = vec_vec_mul(p, p, N);  // This is the way I use the previous r
        beta = (num_beta) / (den_alpha);
        for (e = 1; e < N + 2; e++)
        {
            p[e] = r_z[e] + beta * p[e];
        }
    }

}




        