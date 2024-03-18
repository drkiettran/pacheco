#include <iostream>
#include <string>
#include <mpi.h>

// const int m_val = 4, n_val = 5;
const int m_val = 7, n_val = 5;


// double raw_A[] = {2, 1, 3, 4, 0, 5, -1, 2, -2, 4, 0, 3, 4, 1, 2, 2, 3, 1, -3, 0};
// double raw_A_1[] = {2, 1, 3, 4, 0, 5, -1, 2, -2, 4, 0, 3, 4, 1, 2, 2, 1, 3, 4, 0, 5, -1, 2, -2, 4, 0, 3, 4, 1, 2, 2, 3, 1, -3, 0};
double raw_x[] = {3, 1, 4, 0, 3};

// double new_A[4][5] = {{2, 1, 3, 4, 0}, 
//                       {5, -1, 2, -2, 4}, 
//                       {0, 3, 4, 1, 2}, 
//                       {2, 3, 1, -3, 0}
//                      };

double raw_A[7][5] = {{2, 1, 3, 4, 0}, 
                      {5, -1, 2, -2, 4}, 
                      {0, 3, 4, 1, 2}, 
                      {2, 1, 3, 4, 0}, 
                      {5, -1, 2, -2, 4}, 
                      {0, 3, 4, 1, 2}, 
                      {2, 3, 1, -3, 0}
                     };
void get_array(double **A, int* m, int*n) {
    *m = m_val;
    *n = n_val;
    *A = &raw_A[0][0];
    return;
}

void get_vector(double** x, int m) {
    *x = &raw_x[0];
}

void print_array(int rank, double *A, int m, int n) {
    std::cout << std::endl << "Array A for " << rank << ":" << std::endl;
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << A[i*n + j] << ' ';
        }
        std::cout << std::endl;
    }
}

void print_vector(int rank, double *x, int m) {
    std::cout << std::endl << "Vector x for " << rank << ":" << std::endl;
    for (int i = 0; i < m; i++) {
        std::cout << x[i] << ' ';
    }
    std::cout << std::endl;
}

void get_input(int rank, double** A, double** x, int *m, int *n) {
    get_array(A, m, n);
    get_vector(x, *n);
    print_array(rank, *A, *m, *n);
    print_vector(rank, *x, *n);
}

double* mv_mult(double* A, double* x, int m, int n) {
    std::cout << "Do matrix-vector multiplication" << std::endl;
    double* B = (double*) malloc(m*sizeof(double));
    
    for (int row = 0; row < m; row++) {
        B[row] = 0;
        for (int j = 0; j < n; j++) {
            B[row] += A[row*n + j] * x[j];
        }
    }
    return B;
}

int mv_mult(int argc, char* argv[]) {
    MPI_Comm comm;
    int size, rank;
    double* A, *x;
    int m, n;
    const int root = 0;
    int count, rem;


    MPI_Init(&argc, &argv);

    comm = MPI_COMM_WORLD;

    MPI_Comm_size(comm, &size);
    MPI_Comm_rank(comm, &rank);

    if (rank == 0) {
        get_input(rank, &A, &x, &m, &n);
        if (size > 1) {
            count = m / (size-1);
            rem = m % (size-1);
        } else {
            count = m;
            rem = m;
        }
    }

    // broadcast m, n and x
    MPI_Bcast(&count, 1, MPI_INT, 0, comm);
    MPI_Bcast(&n, 1, MPI_INT, 0, comm);
    

    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            std::cout << "Sending to " << i << " " << count << " row(s)" << std::endl;
            MPI_Send(x, n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
            MPI_Send(&A[(i-1)*count*n], count*n, MPI_DOUBLE, i, 1, MPI_COMM_WORLD);
        }
        std::cout << "to process from " << rem << " rows" << std::endl;
        double* B = mv_mult(&A[count*(size-1)*n], x, rem, n);
        
        if (size > 1) {
            double* C = (double*)malloc(count*sizeof(double));
            for (int i = 1; i < size; i++) {
                MPI_Recv(C, count, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                for (int row = 0; row < count; row++) {
                    std::cout << "rank " << i << ":" << C[row] << std::endl;
                }
            }
        }
        for (int row = 0; row < rem; row++) {
            std::cout << "rank " << rank << ":" << B[row] << std::endl;
        }


    } else {
        std::cout << "alloc " << count*n*sizeof(MPI_DOUBLE) << " bytes" << std::endl;
        x = (double*) malloc(n*sizeof(MPI_DOUBLE));
        A = (double*) malloc(count*n*sizeof(MPI_DOUBLE));
        MPI_Recv(x, n, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(A, count*n, MPI_DOUBLE, root, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        
        print_vector(rank, x, n);
        print_array(rank, A, count, n);

        double* B = mv_mult(A, x, count, n);        
        // for (int row = 0; row < count; row++) {
        //     std::cout << "rank " << rank << " result:" << B[row] << std::endl;
        // }

        MPI_Send(B, count, MPI_DOUBLE, root, 1, MPI_COMM_WORLD);
        std::cout << "rank " << rank << " comletes!" << std::endl;
    }

    MPI_Finalize();
    return 0;
}

int main(int argc, char* argv[]) {
    mv_mult(argc, argv);
    return 0;
}