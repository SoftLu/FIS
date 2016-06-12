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
    );

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
    );

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
    );

// Matrix-vector multiplication:
void mat_vec_mul_CSR(
        double* y,                  // Output: the resulting vector
        const double* const x,      // Input: the vector that is multiplying
        const int* const I_CSR,     // Input: the I_CSR vector of the matrix
        const int* const J_CSR,     // Input: the J_CSR vector of the matrix
        const double* const V_CSR,  // Input: the V_CSR vector of the matrix
        int N                       // Input: the size of the matrix
        );


double vec_vec_mul(
    double* vec1,        // Input: One of the vectors
    double* vec2,        // Input: Other vector
    int N
    );

double modulus(
    const double* const vector,
    const int N
    );

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
    );

// Inverse of a matrix:
void inverse_mul(
    double* inv_m,
    double **M,
    double* g,
    int m
    );

void mat_vec(
    double* vec_fin,
    double** mat,
    double* vec,
    int N                 // debe ir hasta el tamaño de mi matrix
    );

// Conjugate Gradient function:
void CG(
    double* r_z,      // Input: r = b - A*x
    double* x_z,      // Input: x
    int N,
    );