#ifndef dim_support_h
#define dim_support_h
#include "dim_class.h"
//column and row
dim col(dim A){ return *(A.left); }
dim row(dim A){ return *(A.right); }

//get part
dim get_mid(dim A, int len){
	int disp = (A.size-len)/2;
	return array(dim(len)+disp, A);
}
dim get_low(dim A, int len){
	return array(dim(len), A);
}
dim get_high(dim A, int len){
	int disp = A.size-len;
	return array(dim(len)+disp, A);
}
//cut part
dim cut_low(dim A, int k){
	return array(dim(A.size-k)+k, A);
}
dim cut_high(dim A, int k){
	return array(dim(A.size-k), A);
}
//cyclic
_PA_Dim cyclic(_PA_Dim A, int k){
	int disp = (A.size+k)%A.size;
	return array(dim(A.size)+disp, A);
}
dim cyclicCol(dim A, int k){
	int size = A.left->size;
	int disp = (size + k)%size;
	dim C = array(dim(size)+disp, col(A));
	return C * row(A);
}
dim cyclicRow(dim A, int k){
	int size = A.right->size;
	int disp = (size+k)%size;
	dim C = array(dim(size)+disp, row(A));
	return col(A)*C;
}
//multi
dim multi(dim A, int m){
	return array(dim(m*A.size), A);
}
dim divide(dim A, int step){
	return array(dim(A.size/step, step), A);
}
dim split(dim A, int len){
	int size = A.size/len;
	dim C = array(dim(size,len), A);
	dim R = get_low(A, len);
	return (C*R);
}
//major
_PA_Dim block_major(_PA_Dim A, int k){//A is 2D matrix, k is the edge of sub-block, divisible 
	if(A.tag != TAG_DIMS) return NULL;
	dim sub = _PA_Array(dim(k)*dim(k), A);
	int x = A.left->size/k; int y = A.right->size/k;
	dim g = _PA_Array(dim(x,k)*dim(y,k), A);
	return g*sub;		
}
_PA_Dim col_major(_PA_Dim A){//A is 2D matrix //transpose 2 dims
	if(A.tag != TAG_DIMS) return NULL;
	return (A.rd()*A.ld());
}
_PA_Dim centerWin(dim A, int edge=1){
	dim W; if(A.tag!=TAG_DIMS) return W;
	int k = edge; 
	int lsize = A.left->size-(k*2);
	int rsize = A.right->size-(k*2);
	dim L = array(dim(lsize)+k, *A.left);
	dim R = array(dim(rsize)+k, *A.right);
	W = L * R;
	return W;
}
//dimMalloc
void dimMalloc(_PA_Dim D, void **a, int elem_size){
	int size = D.dsize()*elem_size;
	*a = (void*)malloc(size);	
}
void dimMalloc(_PA_Dim D, void **a, void **b, int elem_size){
	int size = D.dsize()*elem_size;
	*a = (void*)malloc(size);	
	*b = (void*)malloc(size);
}
void dimMalloc(_PA_Dim D, void ** a, void **b, void **c, int elem_size){
	int size = D.dsize()*elem_size;
	*a = (void*)malloc(size);	
	*b = (void*)malloc(size);
	*c = (void*)malloc(size);
}
void dimMalloc(_PA_Dim D, void ** a, void **b, void **c, void **d, int elem_size){
	int size = D.dsize()*elem_size;
	*a = (void*)malloc(size);	
	*b = (void*)malloc(size);
	*c = (void*)malloc(size);
	*d = (void*)malloc(size);
}

//print value
template <typename T>
void prtVal(_PA_Dim D, T* A){
	for(int i=0; i<D.size; i++){
		if(D.tag==TAG_DIMS && i%D.right->size==0) cout << endl;
		cout << A[D.offset(i)] << ", ";
	}
	cout << endl;
}
//set value
template <typename T>
void setValIdx(_PA_Dim D, T* a, T val){
	for(int i=0; i<D.size; i++) a[D.offset(i)] = val+i;
}
template <typename T>
void setVal(_PA_Dim D, T* a, T val){
	for(int i=0; i<D.size; i++) a[D.offset(i)] = val;
}
template <typename T>
void setVal(_PA_Dim D, T* a, T* b, T val){
	int idx;
	for(int i=0; i<D.size; i++){
	   	idx = D.offset(i);
		a[idx] = b[idx] = val;
	}
}
template <typename T>
void setVal(_PA_Dim D, T* a, T* b, T* c, T val){
	int idx;
	for(int i=0; i<D.size; i++){ 
		idx = D.offset(i);
		a[idx] = b[idx] = c[idx] = val;
	}
}
template <typename T>
void setVal(_PA_Dim D, T* a, T* b, T* c, T* d, T val){
	int idx;
	for(int i=0; i<D.size; i++){ 
		idx = D.offset(i);
		a[idx] = b[idx] = c[idx] = d[idx] = val;
	}
}
//check pattern
bool hasRoot(_PA_Dim A, int root){
	for(int i=0; i<A.size; i++) if(A.nodes(i)==root) return true;
	return false;
}
void toArr(_PA_Dim A, _PA_Dim *res, int k=0){
	if(A.tag==TAG_DIMS){
		toArr(*A.left, res, k); k+=A.left->dims();
		toArr(*A.right, res, k); return;
	}
	res[k] = A; return;
}
//proc-proc, data-data, with same size => match
bool match(_PA_Dim A, _PA_Dim B){
	if(A.size!=B.size) return false;
	int an = A.dims(); _PA_Dim *arr = (_PA_Dim*)malloc(sizeof(_PA_Dim)*an);
	int bn = B.dims(); _PA_Dim *brr = (_PA_Dim*)malloc(sizeof(_PA_Dim)*bn);
	toArr(A, arr); toArr(B, brr); int ta, tb, sa, sb; int i=0, j=0;
	while(i<an && j<bn){
		sa = arr[i].size; sb = brr[j].size;	 if(sa==1){ i++; continue; } if(sb==1){ j++; continue; }
		ta = arr[i].tag; tb = brr[j].tag;  if(ta != tb || ta>1 || tb>1){ free(arr); free(brr); return false;}
		while(sa!=sb){
			if(sa>sb){ if(tb!= brr[j+1].tag){ free(arr); free(brr); return false;} sb *= brr[j+1].size; j++;}
			if(sa<sb){ if(ta!= arr[i+1].tag){ free(arr); free(brr); return false;} sa *= arr[i+1].size; i++;}
		}
		i++; j++;
	}
	free(arr); free(brr); return true;
}
void procDims(_PA_Dim A, _PA_Dim *P){
	if(A.tag==TAG_DIMS){procDims(*A.left, P); procDims(*A.right, P); return;}
	if(A.tag==TAG_PROC){ P->append(A); return;}
	return;
}
bool sameProcSeq(_PA_Dim A, _PA_Dim B){
	if(A.size!=B.size || A.psize()!=B.psize()) return false;
	for(int i=0; i<A.size; i++) if(A.nodes(i)!=B.nodes(i)) return false;
	return true;
}
//if data movement only ocur in local
bool local(_PA_Dim A, _PA_Dim B){
	if(match(A, B)){
		_PA_Dim P1, P2; procDims(A, &P1); procDims(B, &P2);
		if(sameProcSeq(P1, P2)) return true;
	}
	return false;	
}

#endif
//end of dim_support.h
