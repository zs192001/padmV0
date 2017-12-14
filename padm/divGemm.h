//blk mul
void gemm(dim A, float *a, float *b, float *c){ 
	int n = A.right->size;
	int blk_len = A.left->size;    
    for (int i=0; i<blk_len; i++)
        for(int j=0; j<blk_len; j++)
            for(int k=0; k<n; k++)
                c[i*blk_len+j] += a[i*n+k]*b[k*blk_len+j];
}

//A B different
void gemm(dim A, float *a, dim B, float *b, float *c){ 
	int Ah = A.height(); int Aw = A.width(); int Bw = B.width();
    for (int i=0; i<Ah; i++)
        for(int j=0; j<Bw; j++)
            for(int k=0; k<Aw; k++)
                c[i*Bw+j] += a[i*Aw+k]*b[k*Bw+j];
}

//3d Gemm
void sumGemm(int rank, int blk_len, _PA_Dim A, float* a, _PA_Dim B, float* b, _PA_Dim C, float* c){
	int len = blk_len;
	int m = A.left->size/len;
	int k = A.right->size/len;
	int l = B.right->size/len;
	double trans_time = 0.0; int elem_size = sizeof(float);
	//dims 
	dim Mm = proc(m, k*l); dim Mk = proc(k, l); dim Ml = proc(l, 1);
	dim Blk = array( data(len)*data(len) );
	float *ina, *inb, *amb; 
	dimMalloc(Blk, (void**)&ina, (void**)&inb, (void**)&amb, sizeof(float)); 
	setVal(Blk, ina, inb, amb, (float)0.0);
	
	dim A_Blk = multi(block_major(A, len), l); 
	dim B_Blk = multi(block_major(B, len), m); 
	dim INA = Ml*Mm*Mk*Blk;
	dim INB = Mm*Mk*Ml*Blk;
	trans_time += copyto(A_Blk, a, INA, ina, elem_size);
	trans_time += copyto(B_Blk, b, INB, inb, elem_size);
	
	MPI_Barrier(MPI_COMM_WORLD);	
	if(CHK_FLAG) gemm(Blk, ina, inb, amb);
	MPI_Barrier(MPI_COMM_WORLD);	
	
	dim AMB = Mk*Mm*Ml*Blk;
	dim CP = data(k, C.dsize());
	dim CSUM = CP * block_major(C,len);
	float *csum; dimMalloc(CSUM, (void**)&csum, sizeof(float)); //setVal(CSUM, csum, (float)0.0);
	trans_time += copyto(AMB, amb, CSUM, csum, elem_size);
	if(CHK_FLAG) for(int i=0; i<CP.dsize(); i++) for(int j=0; j<C.dsize(); j++) c[j] += csum[CP.offset(i)+j];
	MPI_Barrier(MPI_COMM_WORLD);	
	free(ina);free(inb);free(amb);free(csum);	
	if(rank==0) printf("thread num %d, all trans time is %.3f \n", _tnum, trans_time);
}
//sum gemm for square matrix
void sumGemm(int rank, int blk_len, _PA_Dim A, float* a, float* b, float* c){
	sumGemm(rank, blk_len, A, a, A, b, A, c);
}
//sum gemm on 64 nodes
void sumGemm_64(int rank, _PA_Dim A, float* a, float* b, float* c){
	int r = 4; int p = 8; int size = A.right->size; int elem_size = sizeof(float);
	int len = (size+r-1)/r; int n = len*r;
	if(size == n){ sumGemm(rank, len, A, a, b, c); return; }

	if(rank==0) printf("patch matrix with zero.\n");
	dim D = proc(p)*data(n/p); dim AA = array(D*D); dim DD = array(data(n/p)*data(n/p));
	dim Win = array(dim(size)*dim(size), AA);
	float *aa, *bb, *cc; 
	dimMalloc(AA, (void**)&aa, (void**)&bb, (void**)&cc, sizeof(float)); 
	setVal(DD, aa, bb, cc, (float)0.0);
	copyto(A, a, Win, aa, elem_size);
	copyto(A, b, Win, bb, elem_size);
	setPrint(1);
	sumGemm(rank, len, AA, aa, bb, cc);
	copyto(Win, cc, A, c, elem_size); 
}

//traditional block algorithm
//A: mlen * klen
//B: klen * llen
//C: mlen * llen
//on m*l processes
void blkGemm(int rank, int blk_len, _PA_Dim A, float* a, _PA_Dim B, float* b, _PA_Dim C, float* c){
	int len = blk_len; 
	int m = A.left->size/len; int l = B.right->size/len;
	int pnum = getProcs(A.psize(), B.psize()); if(pnum/(m*l) == 4){ len /=2; m *=2; l *=2; }
	double trans_time = 0.0; int elem_size = sizeof(float);
	//dims
	dim Mm = proc(m, l); dim Ml = proc(l, 1);
	dim Seg = array(data(len) * data(row(A).size));	
	dim A_Segs = multi(A, l);
	dim INA = Ml * Mm * Seg;

	dim cols = col(B) * get_low(row(B), len);
	dim nums = divide(row(B), len);
	dim B_Segs = multi(nums*cols, m);
	
	dim INB = Mm * Ml * Seg;
	dim Blk = data(len*len);
	dim AMB = Mm * Ml * Blk;
	dim C_Blk = block_major(C,len);
	//malloc
	float *ina, *inb, *amb; 
	dimMalloc(Seg, (void**)&ina, (void**)&inb, sizeof(float));  setVal(Seg, ina, inb, (float)0.0); 
	dimMalloc(Blk, (void**)&amb, sizeof(float)); setVal(Blk, amb, (float)0.0);
	//trans & compute
	trans_time += copyto(A_Segs, a, INA, ina, elem_size);
	trans_time += copyto(B_Segs, b, INB, inb, elem_size);
	MPI_Barrier(MPI_COMM_WORLD);	
	if(CHK_FLAG) gemm(Seg, ina, inb, amb);
	MPI_Barrier(MPI_COMM_WORLD);	
	trans_time += copyto(AMB, amb, C_Blk, c, elem_size);
	MPI_Barrier(MPI_COMM_WORLD);	
	free(ina);free(inb);free(amb);	
	if(rank==0){ 
		printf("thread num %d, all trans time is %.3f \n", _tnum, trans_time);
	}
}

