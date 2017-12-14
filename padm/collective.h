#ifndef collective_h
#define collective_h
#include "trans.h"
#define BCAST_LONG_MSG (1024*1024*8)
void PA_Alltoall(void *s, void *r, int count, int elem_size){
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	_PA_Dim P = proc(pnum);
	_PA_Dim A = array( P * data(pnum*count));
	_PA_Dim B = array( data(pnum) * P * data(count));
	copyto(A, s, B, r, elem_size);
}
void PA_Bcast_2(void *s, int count, int root, int elem_size){
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	_PA_Dim P1 = _PA_Node(root);
	_PA_Dim A = array( P1 * data(count));
	_PA_Dim B = array( proc(pnum) * data(count/pnum));
	void* buf; dimMalloc(B, (void**)&buf, elem_size);
	copyto(A, s, B, buf, elem_size);

	_PA_Dim P2 = array(dim(pnum), proc(1))*proc(pnum);
	_PA_Dim C = P2 * data(count/pnum);
	_PA_Dim D = proc(pnum)*data(count);
	copyto(C, buf, D, s, elem_size);
}
void PA_Bcast_1(void* s, int count, int root, int elem_size){
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	_PA_Dim P = _PA_Node(root, pnum);
	_PA_Dim A = array( P * data(count));
	_PA_Dim B = array( proc(pnum) * data(count));
	copyto(A, s, B, s, elem_size);	
}
void PA_Bcast(void* s, int count, int root, int elem_size){
	if((count*elem_size)>=BCAST_LONG_MSG){
		PA_Bcast_2(s, count, root, elem_size);
	}
	else PA_Bcast_1(s, count, root, elem_size);
	return;
}
void PA_Gather(void* s, void* r, int count, int root, int elem_size){
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	_PA_Dim P = _PA_Node(root);
	_PA_Dim A = array( proc(pnum) * data(count));
	_PA_Dim B = array( P * data(pnum*count));
	copyto(A, s, B, r, elem_size);	
}
void PA_Scatter(void* s, void* r, int count, int root, int elem_size){
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	_PA_Dim P = _PA_Node(root);
	_PA_Dim A = array( P * data(pnum*count));
	_PA_Dim B = array( proc(pnum) * data(count));
	copyto(A, s, B, r, elem_size);	
}
void PA_Allgather(void* s, void* r, int count, int elem_size){
	int pnum; MPI_Comm_size(MPI_COMM_WORLD, &pnum);
	_PA_Dim A = multi(proc(pnum)*data(count), pnum);	
	_PA_Dim B = proc(pnum) * data(pnum*count);
	copyto(A, s, B, r, elem_size);
}
#endif
//end of collective
