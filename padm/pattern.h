#ifndef pattern_h
#define pattern_h
#include "trans.h"
#define LONG (2)
int clpNum(_PA_Dim A){
	if(A.tag == TAG_DIMS) return clpNum(*A.left)+clpNum(*A.right);
	if(A.is_ref){
		if(A.size>A.base->size) return 1;
	}
	return 0;	
}
int clpStep(_PA_Dim A, int k=1){
	if(A.tag == TAG_DIMS){
		return clpStep(*A.left, k*(A.right->size)) *clpStep(*A.right, k);
	}	
	if(A.is_ref) return (k*(A.base->size)); 
	return 1;
}
//split proc vs data //without REF_DIMS
void split(_PA_Dim A, _PA_Dim *P, _PA_Dim *D){
	if(A.tag == TAG_DIMS){ split(*A.left, P, D); split(*A.right, P, D); return;}
	if(A.tag==TAG_PROC){ P->append(A); return; }
	if(A.tag==TAG_DATA){ D->append(A); return; }
	return;
}
bool chunk(_PA_Dim A, int choice){
	int an = A.dims(); _PA_Dim *arr = (_PA_Dim*)malloc(sizeof(_PA_Dim)*an);
	toArr(A, arr); int count=an; _PA_Dim TD;
	for(int i=0; i<count;){
		if(arr[i].tag == choice && !arr[i].is_ref) i++; //no ref dim //debug
		else { count--; for(int j=i; j<count; j++) arr[j]=arr[j+1];}
	}
	//sort asd
	for(int i=0; i<count-1; i++) if(arr[i].step>arr[i+1].step){ TD=arr[i]; arr[i]=arr[i+1]; arr[i+1]=TD;}
	int step = arr[0].step; if(step!=1){ free(arr); return false;} 
	int scope = arr[0].size*step;
	for(int i=1; i<count; i++){ 
		if(scope!=arr[i].step){ free(arr); return false;}
		scope *= arr[i].size;
	}	
	return true;
}
bool chunk_data(_PA_Dim A){ return chunk(A, TAG_DATA);}
bool chunk_proc(_PA_Dim A){ return chunk(A, TAG_PROC);}
//situation: A.dsize=B.dsize; A.psize<B.psize;
//A, B chunk_data; A.psize < B.psize
//only one collapsed dim
bool chk_bcast(_PA_Dim A, _PA_Dim B){
	if(!match(A, B)) return false;
	if(A.psize()<B.psize() && clpNum(A)==1 && chunk_data(A) && chunk_data(B)) return true;
	return false;
}
void BcastGen(_PA_Dim AP, _PA_Dim AD, void *s, _PA_Dim BP, _PA_Dim BD, void *r, int elem_size){
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int dummy = BP.size; int root = AP.root(); 
	_PA_Dim Ra = node(root);
	_PA_Dim Rb = array(dim(1), cyclic(BP, root));
	copyto(Ra*AD, s, Rb*BD, r, elem_size);
	int size = 1; int step = dummy; int mid = root+step/2;
	_PA_Dim TA; _PA_Dim TB;
	for(int i=0; i<log(dummy); i++){
		TA = array(dim(size, step), cyclic(BP, root));
		TB = array(dim(size, step), cyclic(BP, mid));
		copyto(TA*AD, r, TB*BD, r, elem_size);
		size *= 2; step/= 2; mid = root + step/2;
	}
}
void BcastGrpGen(_PA_Dim AP, _PA_Dim AD, void *s, _PA_Dim BP, _PA_Dim BD, void *r, int elem_size){
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0) printf("group general bcast\n");
	int sp = clpStep(AP); int loop = AP.psize(); 
	int mul = AP.size/loop; int dp = loop/sp;
	int k1,k2;
	for(int i=0; i<loop; i++){
		k1 = AP.nodes(i*dp); k2 = BP.nodes(i*dp);
		_PA_Dim A = array(dim(mul,sp), AP+k1);
		_PA_Dim B = array(dim(mul,sp), BP+k2);
		BcastGen(A, AD, s, B, BD, r, elem_size);
	}	
}
void BcastOpt(_PA_Dim AP, _PA_Dim AD, void *s, _PA_Dim BP, _PA_Dim BD, void *r, int elem_size){
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0) printf("opt bcast:\n");
	int dummy = BP.size; int size = BD.size/dummy; 
	_PA_Dim A = (proc(1)+AP.root()) * AD;
	_PA_Dim B = BP * data(size);
	void *buf; dimMalloc(B, (void**)&buf, elem_size); ///////
	copyto(A, s, B, buf, elem_size);
	_PA_Dim C = multi(B, dummy);
	_PA_Dim D = BP * BD;
	copyto(C, buf, D, r, elem_size);
}
void BcastGrpOpt(_PA_Dim AP, _PA_Dim AD, void *s, _PA_Dim BP, _PA_Dim BD, void *r, int elem_size){
	int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	if(rank==0) printf("group opt bcast\n");
	int sp = clpStep(AP); int loop = AP.psize(); 
	int mul = AP.size/loop; int dp = loop/sp;
	int k1,k2;
	for(int i=0; i<loop; i++){
		k1 = AP.nodes(i*dp); k2 = BP.nodes(i*dp);
		_PA_Dim A = array(dim(mul,sp), AP+k1);
		_PA_Dim B = array(dim(mul,sp), BP+k2);
		BcastOpt(A, AD, s, B, BD, r, elem_size);
	}	
}
void pattern_copyto(_PA_Dim A, void *s, _PA_Dim B, void *r, int elem_size){
	bool isb = chk_bcast(A, B);
	_PA_Dim AP, AD, BP, BD;
	split(A, &AP, &AD); split(B, &BP, &BD);	
	if(isb && A.dsize()>LONG){ 
		if(AP.psize()==1) BcastOpt(AP, AD, s, BP, BD, r, elem_size); 
		else BcastGrpOpt(AP, AD, s, BP, BD, r, elem_size);
	}
	else{ 
		if(AP.psize()==1) BcastGen(AP, AD, s, BP, BD, r, elem_size);
		else BcastGrpGen(AP, AD, s, BP, BD, r, elem_size);
	}
}
#endif
