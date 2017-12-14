/* filename: matMultiplyWithMPI.cpp
 * parallel matrix multiplication with MPI
 * C(m,n) = A(m,p) * B(p,n)
 * input: three parameters - m, p, n
 * @copyright: fengfu-chris
 */
#include<iostream>
#include<mpi.h>
#include<math.h>
#include<stdlib.h>
#include "padm/padm.h"
using namespace std;

void initMatrixWithRV(float *A, int rows, int cols);
void matMultiplyWithSingleThread(float *A, float *B, float *matResult, int m, int p, int n);

int main(int argc, char** argv)
{
    int m = 20; //atoi(argv[1]);
    int p = 20; //atoi(argv[2]);
    int n = 20; //atoi(argv[3]);
    
    float *A, *B, *C;
    float *bA, *bC;  

    int myrank, numprocs;

    MPI_Status status;
  
    MPI_Init(&argc, &argv);  // 并行开始
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs); 
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank); 
    
    int bm = m / numprocs;

    bA = new float[bm * p];
    B  = new float[p * n];
    bC = new float[bm * n];

    if(myrank == 0){
        A = new float[m * p];
        C = new float[m * n];
        
        initMatrixWithRV(A, m, p);
        initMatrixWithRV(B, p, n);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    //communication
	dim _A0 = data(m*p);
	dim _Ab = proc(numprocs) * data(bm*p);
	copyto(_A0, A, _Ab, bA, sizeof(float));
	//MPI_Scatter(A, bm * p, MPI_FLOAT, bA, bm *p, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
	dim _B0 = multi(data(p*n), numprocs);
   	dim _Bb = proc(numprocs) * data(p*n);
	copyto(_B0, B, _Bb, B, sizeof(float));	
	//MPI_Bcast(B, p * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
	//compute on each process
	matMultiplyWithSingleThread(bA, B, bC, bm, p, n);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //collect the results
	dim _Cb = proc(numprocs) * data(bm*n);
	dim _C0 = data(m*n);
	copyto(_Cb, bC, _C0, C, sizeof(float));
	//MPI_Gather(bC, bm * n, MPI_FLOAT, C, bm * n, MPI_FLOAT, 0, MPI_COMM_WORLD);
  
    //the left behind blocks
	int remainRowsStartId = bm * numprocs;
    if(myrank == 0 && remainRowsStartId < m){
        int remainRows = m - remainRowsStartId;
        matMultiplyWithSingleThread(A + remainRowsStartId * p, B, C + remainRowsStartId * n, remainRows, p, n);
    }

    delete[] bA;
    delete[] B;
    delete[] bC;
  	
    if(myrank == 0){
        delete[] A;
        delete[] C;
    }
    
	cout << "over." << endl;

    MPI_Finalize(); // 并行结束

    return 0;
}

void initMatrixWithRV(float *A, int rows, int cols)
{
    srand((unsigned)time(NULL));
    for(int i = 0; i < rows*cols; i++){
        A[i] = (float)rand() / RAND_MAX;
    }
}
void matMultiplyWithSingleThread(float *A, float *B, float *matResult, int m, int p, int n)
{
    for(int i=0; i<m; i++){
        for(int j=0; j<n; j++){
            float temp = 0;
            for(int k=0; k<p; k++){
                temp += A[i*p+k] * B[k*n + j];
            }
            matResult[i*n+j] = temp;
        }
    }
}
