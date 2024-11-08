// #include <stdio.h>
// #include <stdlib.h>
// #include <math.h>
// #include <Accelerate/Accelerate.h>

// void print_matrix(const char* desc, int m, int n, double* a, int lda) {
//     int i, j;
//     printf("\n %s\n", desc);
//     for (i = 0; i < m; i++) {
//         for (j = 0; j < n; j++) printf(" %6.2f", a[i*lda+j]);
//         printf("\n");
//     }
// }

// int main() {
//     int m = 1, n = 3;
//     double A[3] = {1.0, 1.0, 1.0};
//     double S[1], U[1], VT[9];
//     double superb[1];
//     int info;
//     int lwork = -1;
//     double wkopt;
//     double* work;

//     // Query and allocate the optimal workspace
//     dgesvd_("A", "A", &m, &n, A, &n, S, U, &m, VT, &n, &wkopt, &lwork, &info);
//     lwork = (int)wkopt;
//     work = (double*)malloc(lwork * sizeof(double));

//     // Perform SVD
//     dgesvd_("A", "A", &m, &n, A, &n, S, U, &m, VT, &n, work, &lwork, &info);

//     if (info > 0) {
//         printf("The algorithm computing SVD failed to converge.\n");
//         exit(1);
//     }

//     // Print the singular values
//     print_matrix("Singular values", 1, 1, S, 1);

//     // Print the right singular vectors (stored row-wise)
//     print_matrix("Right singular vectors (stored row-wise)", n, n, VT, n);

//     // The nullspace is the last n-r columns of VT, where r is the rank of A
//     // In this case, the rank of A is 1, so the nullspace is the last 2 columns of VT
//     printf("\nNullspace:\n");
//     for (int i = 0; i < n; i++) {
//         for (int j = 1; j < n; j++) {
//             printf(" %6.2f", VT[i*n + j]);
//         }
//         printf("\n");
//     }

//     free(work);
//     return 0;
// }
#include <stdio.h>
#include <stdlib.h>
#include <Accelerate/Accelerate.h>

void print_matrix(const char* desc, int m, int n, double* a, int lda) {
    int i, j;
    printf("\n %s\n", desc);
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) printf(" %6.5f", a[i*lda+j]);
        printf("\n");
    }
}

void compute_nullspace(int n) {
    int m = 1;
    double* A = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        A[i] = 1.0;
    }
    double* S = (double*)malloc(n * sizeof(double));
    double* U = (double*)malloc(m * m * sizeof(double));
    double* VT = (double*)malloc(n * n * sizeof(double));
    int info;
    int lwork = -1;
    double wkopt;
    double* work;

    // Query and allocate the optimal workspace
    dgesvd_("A", "A", &m, &n, A, &n, S, U, &m, VT, &n, &wkopt, &lwork, &info);
    lwork = (int)wkopt;
    work = (double*)malloc(lwork * sizeof(double));

    // Perform SVD
    dgesvd_("A", "A", &m, &n, A, &n, S, U, &m, VT, &n, work, &lwork, &info);

    if (info > 0) {
        printf("The algorithm computing SVD failed to converge.\n");
        exit(1);
    }

    // Print the singular values
    print_matrix("Singular values", 1, n, S, 1);

    // Print the right singular vectors (stored row-wise)
    print_matrix("Right singular vectors (stored row-wise)", n, n, VT, n);

    // The nullspace is the last n-1 columns of VT
    printf("\nNullspace:\n");
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < n; j++) {
            printf(" %6.5f", VT[i*n + j]);
        }
        printf("\n");
    }

    free(A);
    free(S);
    free(U);
    free(VT);
    free(work);
}

int main() {
    int n;
    printf("Enter the value of n: ");
    scanf("%d", &n);

    compute_nullspace(n);

    return 0;
}