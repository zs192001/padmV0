#ifndef trans_h
#define trans_h
#include <mpi.h>
#include "dim.h"
#include "func.h"
#include "scanComm.h"
#include "bufComm.h"
#include "local.h"
double copyto(dim A, void *s, dim B, void *r, int elem_size){
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	int choice = 0;	struct timespec start, stop; double accum; 
	MPI_Barrier(MPI_COMM_WORLD);
	//size match?
	if(rank==0 && _prt_flag){ 
		printf("from:"); A.print(); printf("to  :"); B.print();
		if(A.size!=B.size){ printf("size not match: %d, %d.\n", A.size, B.size); return 0.0; }
		printf("size is %lld.\n", A.size);
	}
	//set args
	setArgs(A, B, elem_size, pnum); if(rank==0 && _prt_flag) prtArgs();
	MPI_Barrier(MPI_COMM_WORLD); clock_gettime(CLOCK_REALTIME, &start);
	//local
	if(local(A, B)){
		if(rank==0) printf("local data copy.\n");
		localComm(A, s, B, r);	
	}	
	//comm
	else if(_len>=_gran){ 
		scanComm(rank, A, s, B, r);
	}
	else{ choice = 1; 
		bufComm(rank, A, s, B, r);
	}
	MPI_Barrier(MPI_COMM_WORLD); clock_gettime(CLOCK_REALTIME, &stop);
 	accum = (stop.tv_sec-start.tv_sec)+(stop.tv_nsec-start.tv_nsec)/1e9; 
    if(rank==0 && _prt_flag){ 
		printf("----copyto, ");
		if(choice) printf("bufComm:"); else printf("scanComm:");
		printf("----%.3f s, \n", accum);
	}
	return accum;
}
double scanCopy(dim A, void *s, dim B, void *r, int elem_size){
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	struct timespec start, stop; double accum; 
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && _prt_flag){ 
		printf("from:"); A.print(); printf("to  :"); B.print();
		if(A.size!=B.size){ printf("size not match: %d, %d.\n", A.size, B.size); return 0.0; }
		printf("size is %lld.\n", A.size);
	}
	setArgs(A, B, elem_size, pnum); if(rank==0 && _prt_flag) prtArgs();
	MPI_Barrier(MPI_COMM_WORLD); clock_gettime(CLOCK_REALTIME, &start);
	scanComm(rank, A, s, B, r);
	MPI_Barrier(MPI_COMM_WORLD); clock_gettime(CLOCK_REALTIME, &stop);
 	accum = (stop.tv_sec-start.tv_sec)+(stop.tv_nsec-start.tv_nsec)/1e9; 
    if(rank==0 && _prt_flag) printf("----scanTrans, %.3f. \n", accum);
	return accum;
}
double bufTrans(dim A, void *s, dim B, void *r, int elem_size){
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	struct timespec start, stop; double accum; 
	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0 && _prt_flag){ 
		printf("from:"); A.print(); printf("to  :"); B.print();
		if(A.size!=B.size){ printf("size not match: %d, %d.\n", A.size, B.size); return 0.0; }
		printf("size is %lld.\n", A.size);
	}
	setArgs(A, B, elem_size, pnum); if(rank==0 && _prt_flag) prtArgs();
	MPI_Barrier(MPI_COMM_WORLD); clock_gettime(CLOCK_REALTIME, &start);
	bufComm(rank, A, s, B, r);
	MPI_Barrier(MPI_COMM_WORLD); clock_gettime(CLOCK_REALTIME, &stop);
 	accum = (stop.tv_sec-start.tv_sec)+(stop.tv_nsec-start.tv_nsec)/1e9; 
    if(rank==0 && _prt_flag) printf("----bufTrans, %.3f. \n", accum);
	return accum;
}
#endif
//end of trans.h
