const char *dgemm_desc = "Naive, three-loop dgemm.";
void square_dgemm(int n, double *A, double *B, double *C) {
    // Initialize C to zero
    for (int i = 0; i < n * n; ++i) {
        C[i] = 0.0;
    }
    // Perform matrix multiplication
    for (int i = 0; i < n; ++i) {
        for (int k = 0; k < n; ++k) {
            double a_ik = A[i + k * n]; // Accessing A(i,k)
            for (int j = 0; j < n; ++j) {
                C[i + j * n] += a_ik * B[k + j * n]; // Updating C(i,j)
            }
        }
    }
}

