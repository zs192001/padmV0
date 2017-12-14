#ifndef local_copy_h
#define local_copy_h
#include <string.h>
#include "dim.h"
#include "func_c.h"
void copy(int tid, _PA_Dim A, void *s, _PA_Dim B, void *r){ 
	int start = tid*_scope; int end = start+_scope; if(end>A.size) end = A.size;
	int sl, rl; int msize = _len*_elem_size;
	for(int c=start; c<end; c+=_len){
        sl = A.offset(c); rl = B.offset(c);
		memcpy(r+rl*_elem_size, s+sl*_elem_size, msize);
    }   
    __sync_fetch_and_add(&packs_count, 1); 
    if(packs_count == _tnum) { 
		pack_over=true; 
		//printf("copy over: tid %d\n", tid);
	}
	return;
}
void localComm(_PA_Dim A, void *s, _PA_Dim B, void *r){
	int scansize = A.size; int len=_scope; 
	launchWork();
    #pragma omp parallel num_threads(1+_tnum)
    {   int ttid = omp_get_thread_num();
        if(ttid==0){}
        else{ int tid=ttid-1;
            copy(tid, A, s, B, r);
        }
    }
	clearWork();
}
double localcopy(_PA_Dim A, void *s, _PA_Dim B, void *r, int elem_size){
	struct timespec start, stop; double accum; 
	if(A.size!=B.size) return 0.0;
	if(_prt_flag){ printf("from:"); A.print(); printf("to:"); B.print();}
	setArgs(A, B, elem_size); if(_prt_flag) prtArgs();
	clock_gettime(CLOCK_REALTIME, &start);
	localComm(A, s, B, r);
	clock_gettime(CLOCK_REALTIME, &stop);
	accum = (stop.tv_sec-start.tv_sec)+(stop.tv_nsec-start.tv_nsec)/1e9;
	if(_prt_flag) printf("----localTrans, %.3f s \n", accum);
	return accum;
}
#endif
//end of local.h
