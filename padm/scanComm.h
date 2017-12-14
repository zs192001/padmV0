#ifndef scan_comm_h
#define scan_comm_h
//comm functions
int testReq(int launch, MPI_Request* req){
	int count = 0;
	for(int i = 0; i<launch; i++) if(msgTest(&req[i])){
		count++; req[i] = MPI_REQUEST_NULL;
	}
	return count;
}
int launchR(MPI_Request* req){
	int count = 0;
	int size = _gran * _elem_size;
	struct msg_info info; int dest, idx; int tag; void *buf = NULL;
	while(!recvMsgEmpty()){
		info = popQue(recv_msgs); getInfo(info, dest, tag, &buf); 
		idx = rbuf_idx.front(); rbuf_idx.pop();
		MPI_Irecv(buf, size, MPI_BYTE, dest, tag, MPI_COMM_WORLD, &req[idx]);
		//printf("testprt: launch recv, source %d, tag %d, size %d,\n", dest, tag, size);
		count++;
	}
	return count;
}
int launchS(MPI_Request* req){
	int count = 0;
	int size = _gran * _elem_size;
	struct msg_info info; int dest, idx; int tag; void *buf = NULL;
	while(!sendMsgEmpty()){
		info = popQue(send_msgs); getInfo(info, dest, tag, &buf); 
		idx = sbuf_idx.front(); sbuf_idx.pop();
		MPI_Isend(buf, size, MPI_BYTE, dest, tag, MPI_COMM_WORLD, &req[idx]);
		count++;
	}
	return count;
}
//_______________________________________________________________work functions
void comm(int rank){
	int num = _reqs; 
	int send_launch_count=0; int recv_launch_count=0; int send_count=0; int recv_count=0;
    bool send_launch_over = false; bool recv_launch_over = false; bool send_over = false; bool recv_over = false;
    MPI_Request* reqs_send = (MPI_Request*)malloc(sizeof(MPI_Request)*num); for(int i=0; i<num; i++) reqs_send[i]=MPI_REQUEST_NULL;
 	MPI_Request* reqs_recv = (MPI_Request*)malloc(sizeof(MPI_Request)*num); for(int i=0; i<num; i++) reqs_recv[i]=MPI_REQUEST_NULL;
	do{
        if(!recv_over){
            if(recv_launch_over && recv_count == recv_launch_count) recv_over=true;
			else recv_count += testReq(recv_launch_count, reqs_recv);
        }
        if(!send_over){
            if(send_launch_over && send_count == send_launch_count) send_over=true;
            else send_count += testReq(send_launch_count, reqs_send); 
        }
        if(!recv_launch_over){
            if(pack_over && recvMsgEmpty()) recv_launch_over = true;
			else recv_launch_count += launchR(reqs_recv);
        }
        if(!send_launch_over){
			if(pack_over && sendMsgEmpty()) send_launch_over = true;
			else send_launch_count += launchS(reqs_send);
        }
    }while(!recv_over || !send_over);
	//printf("comm over: rank %d, send %d, recv %d\n", rank, send_count, recv_count);
	free(reqs_send); free(reqs_recv);
	return;
}
void scan(int rank, int tid, _PA_Dim A, void *s, _PA_Dim B, void *r){ 
	long long int start = tid*_scope; 
	long long int end = start+_scope; if(end>A.size) end = A.size;
	int from, to; int sl, rl; 
    struct msg_info info; int tag; void *buf; void *tr,*ts; 
	int elem = _elem_size; int cpy_size = _gran*elem;
	int** tag_all = (int**)malloc(_procs*sizeof(int*));
	for(int i=0; i<_procs; i++){
		tag_all[i] = (int*)malloc(_procs*sizeof(int));
		for(int j=0; j<_procs; j++) tag_all[i][j] = 0;
	}
	for(long long int c=start; c<end; c+=_gran){
		from = A.nodes(c); to = B.nodes(c);
        sl = A.offset(c); rl = B.offset(c);
		tag_all[to][from]++;	
		//printf("rank %d, from %d, to %d, sl %d, rl %d\n", rank, from, to, sl, rl);
		if(from==rank && to==rank){
			tr = r+rl*elem; ts = s+sl*elem;
			if(tr != ts) memcpy(tr, ts, cpy_size);
		}
		else{
			if(from == rank){
            	tag = tag_all[to][from]; 
				buf = s+sl*elem;
				//printf("send: from %d to %d, tag %d, \n",rank, to, tag);
				regInfo(info, to, tag, buf); pushQue(send_msgs, info); 
        	}   
        	if(to == rank){
        	    tag = tag_all[to][from]; 
				buf = r+rl*elem;
				//printf("recv: from %d to %d, tag %d, \n",from, rank, tag);
        	    regInfo(info, from, tag, buf); pushQue(recv_msgs, info); 
        	}   
		}
    } 
    __sync_fetch_and_add(&packs_count, 1); 
    if(packs_count == _tnum) { 
		pack_over=true; 
		//printf("scan over: rank %d\n", rank);
	}
   	for(int i=0; i<_procs; i++) free(tag_all[i]); free(tag_all);	
	return;
}
void scanComm(int rank, _PA_Dim A, void *s, _PA_Dim B, void *r){
	launchWork();
    MPI_Barrier(MPI_COMM_WORLD);
    #pragma omp parallel num_threads(1+_tnum)
    {   int ttid = omp_get_thread_num();
        if(ttid==0) comm(rank);
        else{ int tid=ttid-1;
            scan(rank, tid, A, s, B, r);
        }
    }
	MPI_Barrier(MPI_COMM_WORLD);
	clearWork();
}
#endif
//end of scanComm.h
