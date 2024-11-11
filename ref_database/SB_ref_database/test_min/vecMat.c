#include <stdio.h>
#include <stdlib.h>

// gcc -o vector_matrix_multiplication vecMat.c -lm

// Function to perform vector-matrix multiplication
void vector_matrix_multiplication(double* v, double** M, double* result, int vector_size, int matrix_cols) {
    for (int j = 0; j < matrix_cols; j++) {
        result[j] = 0.0;
        for (int i = 0; i < vector_size; i++) {
            result[j] += v[i] * M[i][j];
        }
    }
}
// Function to perform matrix-vector multiplication
void matrix_vector_multiplication(double** M, double* v, double* result, int matrix_rows, int matrix_cols) {
    for (int i = 0; i < matrix_rows; i++) {
        result[i] = 0.0;
        for (int j = 0; j < matrix_cols; j++) {
            result[i] += M[i][j] * v[j];
        }
    }
}
int main() {
    int vector_size = 4;
    int matrix_cols = 3;

    // Define a vector
    double v[] = {1.0, 2.0, 3.0,4.0};

    // Define a matrix
    double** M = (double**)malloc(vector_size * sizeof(double*));
    for (int i = 0; i < vector_size; i++) {
        M[i] = (double*)malloc(matrix_cols * sizeof(double));
    }

    // Initialize the matrix
    M[0][0] = 1.0; M[0][1] = 2.0; M[0][2] = 3.0;
    M[1][0] = 4.0; M[1][1] = 5.0; M[1][2] = 6.0;
    M[2][0] = 7.0; M[2][1] = 8.0; M[2][2] = 9.0;
    M[3][0] = 3.0; M[3][1] = 3.0; M[3][2] = 3.0;
    // Allocate memory for the result
    double* result = (double*)malloc(matrix_cols * sizeof(double));
    vector_matrix_multiplication(v, M, result, vector_size, matrix_cols);

    // Allocate memory for the result of matrix-vector multiplication
    double* result2 = (double*)malloc(vector_size * sizeof(double));
    matrix_vector_multiplication(M, result, result2, vector_size, matrix_cols);



    // Print the result
    printf("Result of v' * M:\n");
    for (int j = 0; j < matrix_cols; j++) {
        printf("%6.2f ", result[j]);
    }
    printf("\n");

    // Print the result of matrix-vector multiplication
    printf("Result of M * result:\n");
    for (int i = 0; i < vector_size; i++) {
        printf("%6.2f ", result2[i]);
    }
    printf("\n");

    // Free allocated memory
    for (int i = 0; i < vector_size; i++) {
        free(M[i]);
    }
    free(M);
    free(result);
    free(result2);
    
    return 0;
}

// # Define the vector
// v = [1.0, 2.0, 3.0,4.0]

// # Define the matrix
// M = [1.0 2.0 3.0;
//      4.0 5.0 6.0;
//      7.0 8.0 9.0;
//      3.0 3.0 3.0]

// # Perform vector-matrix multiplication
// result = v' * M
