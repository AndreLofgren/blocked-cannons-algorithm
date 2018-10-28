#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

void fillMatrixRow(double *M, int n);
void mul(double *C, double *A, double *B, int n);
int isEqual(double *A, double *B, int n);
int isEqualBlock(double *A, double *B, int n, int nx);
void printMatrixRow(double *M, int n);
void printMatrixBlock(double *M, int n, int nx);
double timer();

static const int root = 0;

int main(int argc, char *argv[]) {
    if (argc < 2)  {
        fprintf(stderr, "usage: mpirun -np p ./matrixMul n\n");
        exit(EXIT_FAILURE);
    }
    int n = atoi(argv[1]);

    int i, j, k;
    int nproc, rank, grid_rank, row_rank, col_rank, px, py;
    int coords[2], pos[2], reorder = 1, ndim = 2;
    int dims[2] = {0, 0}, rowdims[2] = {0, 1}, coldims[2] = {1, 0}, periods[2] = {1, 1};

    MPI_Init(&argc, &argv);
    MPI_Status status[4];
    MPI_Request request[4];
    MPI_Comm GRID, ROW, COL;

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    dims[0] = sqrt(nproc); px = dims[0];
    dims[1] = dims[0]; py = dims[1];

    int nx = n/px; int ny = n/py;

    /* Check if sqrt(p) is integer */
    if (px*py != nproc) {
        if (rank == root) {
            fprintf(stderr, "sqrt(p) is not an integer\n");
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }

    /* Create new data type representing a submatrix */
    MPI_Datatype Block;
    int count = nx, blocklen = ny, stride = n;
    MPI_Type_vector(count, blocklen, stride, MPI_DOUBLE, &Block);
    MPI_Type_commit(&Block);

    /* Split communicator into a Cartesian topology */
    MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, reorder, &GRID);
    MPI_Comm_rank(GRID, &grid_rank);
    MPI_Cart_coords(GRID, grid_rank, ndim, coords);
    MPI_Cart_sub(GRID, rowdims, &ROW);
    MPI_Cart_sub(GRID, coldims, &COL);
    MPI_Comm_rank(ROW, &row_rank);
    MPI_Comm_rank(COL, &col_rank);

    /* Allocate matrices */
    double *A, *B, *C, *C_ref;
    double time;
    if (grid_rank == root) {
        A = malloc(n*n*sizeof(double));
        B = malloc(n*n*sizeof(double));
        C = calloc(n*n, sizeof(double));
        C_ref = calloc(n*n, sizeof(double)); /* reference */
        srand48(10);
        fillMatrixRow(A, n);
        srand48(20);
        fillMatrixRow(B, n);
        time = timer();
        mul(C_ref, A, B, n);
        printf("serial time %1.4lf sec\n", timer() - time);
    }

    /* Start timer */
    if (grid_rank == root) {
        time = timer();
    }

    /* Allocate block matrices*/
    double *subMatA = malloc(nx*ny*sizeof(double));
    double *subMatB = malloc(nx*nx*sizeof(double));
    double *subMatC = calloc(nx*ny, sizeof(double));
    double *bufferA = malloc(nx*ny*sizeof(double));
    double *bufferB = malloc(nx*ny*sizeof(double));

    /* Distribute blocks to processors */
    int sendA, sendB, recvA, recvB;
    int coordA[2], coordB[2];
    int index;
    if (grid_rank == root) {
        for (i = 0; i < py; i++) {
            for (j = 0; j < px; j++) {
                coordA[0] = i;
                coordA[1] = j-i;
                coordB[0] = i-j;
                coordB[1] = j;
                MPI_Cart_rank(GRID, coordA, &sendA);
                MPI_Cart_rank(GRID, coordB, &sendB);
                index = i*n*ny + j*nx;
                MPI_Isend(&A[index], 1, Block, sendA, 0, GRID, &request[0]);
                MPI_Isend(&B[index], 1, Block, sendB, 1, GRID, &request[1]);
            }
        }
    }
    MPI_Irecv(&subMatA[0], nx*ny, MPI_DOUBLE, root, 0, GRID, &request[0]);
    MPI_Irecv(&subMatB[0], nx*ny, MPI_DOUBLE, root, 1, GRID, &request[1]);
    MPI_Wait(&request[0], MPI_STATUSES_IGNORE);
    MPI_Wait(&request[1], MPI_STATUSES_IGNORE);

    /* Don't need A and B anymore */
    MPI_Barrier(GRID);
    if (grid_rank == root) {
        free(A);
        free(B);
    }

    /* Cannon's algorithm */
    for (k = 0; k < px; k++) {
        int i, j;

        sendA = row_rank-1 > -1 ? row_rank-1 : px + row_rank-1;
        recvA = (row_rank+1) % px;
        sendB = col_rank-1 > -1 ? col_rank-1 : py + col_rank-1;
        recvB = (col_rank+1) % py;

        MPI_Isend(subMatA, nx*ny, MPI_DOUBLE, sendA, 0, ROW, &request[0]);
        MPI_Isend(subMatB, nx*ny, MPI_DOUBLE, sendB, 1, COL, &request[1]);
        MPI_Irecv(bufferA, nx*ny, MPI_DOUBLE, recvA, 0, ROW, &request[2]);
        MPI_Irecv(bufferB, nx*ny, MPI_DOUBLE, recvB, 1, COL, &request[3]);

        mul(subMatC, subMatA, subMatB, nx);

        MPI_Wait(&request[0], MPI_STATUSES_IGNORE);
        MPI_Wait(&request[1], MPI_STATUSES_IGNORE);
        MPI_Wait(&request[2], MPI_STATUSES_IGNORE);
        MPI_Wait(&request[3], MPI_STATUSES_IGNORE);

        /* switch buffers */
        double *tempA = bufferA;
        double *tempB = bufferB;
        bufferA = subMatA;
        bufferB = subMatB;
        subMatA = tempA;
        subMatB = tempB;
    }

    /* Collect Blocks */
    MPI_Isend(subMatC, nx*ny, MPI_DOUBLE, root, grid_rank, GRID, &request[0]);
    if (grid_rank == root) {
        int index, coord[2];
        for (k = 0; k < nproc; k++) {
            MPI_Cart_coords(GRID, k, ndim, coord);
            i = coord[0]; j = coord[1];
            index = i*n*ny + j*nx;
            MPI_Recv(&C[index], 1, Block, k, k, GRID, &status[0]);
        }
    }
    MPI_Wait(&request[0], &status[0]);

    MPI_Barrier(GRID);
    if (grid_rank == root) {
        printf("parallel time %1.4lf sec\n", timer() - time);
        printf("correct: %s\n", isEqual(C, C_ref, n) ? "yes" : "no");
        free(C); free(C_ref);
    }
    free(subMatA); free(subMatB); free(subMatC); free(bufferA); free(bufferB);

    MPI_Finalize();
    return 0;
}

void fillMatrixRow(double *M, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            M[i*n+j] = drand48();
        }
    }
}

void mul(double *C, double *A, double *B, int n) {
    int i, j, k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                C[i*n+j] += A[i*n+k]*B[k*n+j];
            }
        }
    }
}

int isEqual(double *A, double *B, int n) {
    int i, j;
    double tol = 1e-8;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (fabs(A[i*n + j] - B[i*n + j]) > tol) {
                return 0;
            }
        }
    }
    return 1;
}

int isEqualBlock(double *A, double *B, int n, int nx) {
    int i, j;
    double tol = 1e-8;
    for (i = 0; i < nx; i++) {
        for (j = 0; j < nx; j++) {
            if (fabs(A[i*nx + j] - B[i*n + j]) > tol) {
                return 0;
            }
        }
    }
    return 1;
}

void printMatrixRow(double *M, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        printf("[");
        for (j = 0; j < n; j++) {
            if (j == n-1) {
                printf("%1.2lf]\n", M[i*n + j]);
            } else {
                printf("%1.2lf, ", M[i*n + j]);
            }
        }
    }
}

void printMatrixBlock(double *M, int n, int nx) {
    int i, j;
    for (i = 0; i < nx; i++) {
        printf("[");
        for (j = 0; j < nx; j++) {
            if (j == nx-1) {
                printf("%1.2lf]\n", M[i*n + j]);
            } else {
                printf("%1.2lf, ", M[i*n + j]);
            }
        }
    }
}

double timer() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}
