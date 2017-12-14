#ifndef buf_comm_h
#define buf_comm_h
struct reqRInfo{ int tid; int blk; bool ready; void* ptr; MPI_Request req; };
struct reqSInfo{ bool ready; void *ptr; MPI_Request req; };
//_______________________________________________________________work functions
int testRecv(int rank, int launch, reqRInfo* info, bool* flag, void **adds){
	int count = 0; int tid, blk; void *ptr;
	for(int i = 0; i<launch && !info[i].ready; i++) if(msgTest(&info[i].req)){
		count++; info[i].req = MPI_REQUEST_NULL;
		tid = info[i].tid; blk = info[i].blk; ptr = info[i].ptr;
		markAndPost(tid, blk, ptr, flag, adds); //check_recv		
	}
	return count;
}
int testSend(int rank, int launch, reqSInfo* info){
	int count = 0;
	for(int i = 0; i<launch && !info[i].ready; i++) if(msgTest(&info[i].req)){
		count++; info[i].req = MPI_REQUEST_NULL;
		pushBufAdd(info[i].ptr);
	}
	return count;
}
int launchRecv(int rank, reqRInfo* info, int base){
	int idx=base; int dest, tag, len, tid, blk; void *ptr;
	while(!recvMsgEmpty()){
		getRecvMsg(dest, tag, len, &ptr, tid, blk);
		info[idx].tid = tid; info[idx].blk = blk; 
		info[idx].ready = false; info[idx].ptr = ptr;
		info[idx].req = MPI_REQUEST_NULL;
		MPI_Irecv(ptr, len*_elem_size, MPI_BYTE, dest, tag, MPI_COMM_WORLD, &info[idx].req);
		idx++;
	}
	return (idx-base);
}
int launchSend(int rank, reqSInfo* info, int base){
	int idx=base; int dest, tag, len; void *ptr; 
	while(!sendMsgEmpty()){
		getSendMsg(dest, tag, len, &ptr);
		info[idx].ready = false; info[idx].ptr = (void*)ptr; info[idx].req = MPI_REQUEST_NULL;
		MPI_Isend(ptr, len*_elem_size, MPI_BYTE, dest, tag, MPI_COMM_WORLD, &info[idx].req);
		idx++;
	}
	return (idx-base);
}
void comm_3(int rank, bool* flag, void **adds){
	int num = _reqs;
	int send_launch_count=0; int recv_launch_count=0; int send_count=0; int recv_count=0;
    bool send_launch_over = false; bool recv_launch_over = false; bool send_over = false; bool recv_over = false;
	
	struct reqSInfo* send_info; send_info = (reqSInfo*)malloc(num*sizeof(reqSInfo)); 
    struct reqRInfo* recv_info; recv_info = (reqRInfo*)malloc(num*sizeof(reqRInfo)); 
	do{
        if(!recv_over){
            if(recv_launch_over && recv_count == recv_launch_count) recv_over=true;
            else recv_count += testRecv(rank, recv_launch_count, recv_info, flag, adds); 
        }
        if(!send_over){
            if(send_launch_over && send_count == send_launch_count) send_over=true;
            else send_count += testSend(rank, send_launch_count, send_info); 
        }
        if(!recv_launch_over){
            if(pack_over && recv_msgs.empty()) recv_launch_over = true;
			else recv_launch_count += launchRecv(rank, recv_info, recv_launch_count);
        }
        if(!send_launch_over){
            if(pack_over && send_msgs.empty()) send_launch_over = true; 
			else send_launch_count += launchSend(rank, send_info, send_launch_count);
        }
    }while(!recv_over || !send_over);
	//printf("rank %d, comm over.recv %d, send %d,\n", rank, recv_count, send_count);
    free(send_info); free(recv_info);
	return;
}
void scanPack_3(int rank, int tid, _PA_Dim A, void *s, _PA_Dim B){
	long long int start = tid*_scope; long long int end = min(start+_scope, A.size); 
	int cnt = _len; int elem = _elem_size; int cnt_bytes = cnt*elem; int pz = _procs * sizeof(int);
	int cur_block = 0; int from, to, loc; 
    int* count_s = (int*)malloc(pz); memset(count_s, 0, pz); 
    int* count_r = (int*)malloc(pz); memset(count_r, 0, pz); 
	int* tag_s = (int*)malloc(pz); 
	int *tag_r = (int*)malloc(pz); int* tid_r = (int*)malloc(pz); int *blk_r = (int*)malloc(pz);
	void **ptr_s; ptr_s=(void**)malloc(_procs*sizeof(void*)); 
	void **ptr_r; ptr_r=(void**)malloc(_procs*sizeof(void*)); 
	int** count_all = (int**)malloc(_procs*sizeof(int*));
	int** tag_all = (int**)malloc(_procs*sizeof(int*));
	for(int i=0; i<_procs; i++){
		count_all[i] = (int*)malloc(_procs*sizeof(int));
		tag_all[i] = (int*)malloc(_procs*sizeof(int));
		for(int j=0; j<_procs; j++){ count_all[i][j]=0; tag_all[i][j]=0; }
	}
	for(long long int c=start; c<end; c+=cnt){
		from = A.nodes(c); to = B.nodes(c); loc = A.offset(c);
		if(count_all[to][from]==0) tag_all[to][from]++;
		count_all[to][from] += cnt;
		if(count_all[to][from] == _gran) count_all[to][from] = 0;
		if(from == rank){
            if(count_s[to]==0){ 
				ptr_s[to] = getBufAdd(_gran*elem); 
				tag_s[to] = tag_all[to][from]; //c; 
			}
            memcpy(ptr_s[to]+count_s[to]*elem, s+loc*elem, cnt_bytes); 
			count_s[to] += cnt;
            if(count_s[to] == _gran){
				regSendMsg(to, tag_s[to], _gran, ptr_s[to]); 
				count_s[to] = 0;
			}
		}
		if(to == rank){
			if(count_r[from] == 0){ 
				ptr_r[from] = getBufAdd(_gran*elem);
				tag_r[from] = tag_all[to][from]; //c; 
				tid_r[from] = tid; blk_r[from] = cur_block++; 
			} 
            count_r[from]+=cnt;
            if(count_r[from] == _gran){ 
				regRecvMsg(from, tag_r[from], _gran, ptr_r[from], tid_r[from], blk_r[from]); 
				count_r[from] = 0;
			}
		}
    }
	for(int i=0; i<_procs; i++){//smaller blocks
		if(rank+i==0 && count_s[0]>0  && _prt_flag) printf("pack_max %d <= gran size %d\n", count_s[0], _gran);
		if(count_s[i]>0){
			regSendMsg(i, tag_s[i], count_s[i], ptr_s[i]); 
			count_s[i]=0; 
		}
		if(count_r[i]>0){ 
			regRecvMsg(i, tag_r[i], count_r[i], ptr_r[i], tid_r[i], blk_r[i]); 
			count_r[i] = 0;
		}
	}
    __sync_fetch_and_add(&packs_count, 1);
    if(packs_count == _tnum){ 
		pack_over=true;
		//if(rank == 0) printf("pack over.\n");
	} 
	free(count_s); free(tag_s); free(count_r); free(tag_r); free(tid_r); free(blk_r);
	for(int i=0; i<_procs; i++){ ptr_s[i]=NULL; ptr_r[i]=NULL; free(count_all[i]); free(tag_all[i]); } 
	ptr_s = NULL; ptr_r = NULL; 
	free(ptr_s); free(ptr_r); free(count_all); free(tag_all);
    return;
}
void unpack_3(int rank, int tid, _PA_Dim A, _PA_Dim B, void *r, bool *flag, void **adds){
	long long int start = tid*_scope; long long int end = start+_scope; if(end>A.size) end = A.size;
	int cnt = _len; int from; int sp, rp, rl;
	int cur_block = 0; void *buf[_procs]; 
	int *count = (int*)malloc(_procs*sizeof(int)); memset(count, 0, _procs*sizeof(int)); 
	for(long long int c=start; c<end; c+=cnt){
        sp = A.nodes(c); rp = B.nodes(c); rl = B.offset(c); from = sp;
        if(rp == rank){
            if(count[from] == 0){
				buf[from] = getUnpackBuf(tid, cur_block, flag, adds);
				cur_block++;
            }
			memcpy(r+rl*_elem_size, buf[from]+count[from]*_elem_size, cnt*_elem_size); 
			count[from]+=cnt;
            if(count[from] == _gran) count[from] = 0;
        }
    }
    return;
}
void bufComm(int rank, _PA_Dim A, void *s, _PA_Dim B, void *r){
	launchWork();
	void *buf; buf=(void*)malloc(_buf_bytes); memset(buf, 0, _buf_bytes); 
	int buf_num = initBufAdd(buf); initList(_reqs*2-buf_num);
	int nums = _tnum*_blks;
	bool* recv_flag; recv_flag = (bool*)malloc(nums);
	void **unpack_adds; unpack_adds = (void**)malloc(sizeof(void*)*nums);
	for(int i=0; i<nums; i++){ recv_flag[i] = false; unpack_adds[i] = NULL; }

	MPI_Barrier(MPI_COMM_WORLD);
    #pragma omp parallel num_threads(1+_tnum*2)
    {   int ttid = omp_get_thread_num();
        if(ttid==0){
			comm_3(rank, recv_flag, unpack_adds);
		}
		else if(ttid<_tnum+1){ int tid=ttid-1;
            scanPack_3(rank, tid, A, s, B);
        }
		else{ int tid=ttid-1-_tnum;
			unpack_3(rank, tid, A, B, r, recv_flag, unpack_adds);
		}
    }
    MPI_Barrier(MPI_COMM_WORLD);
	free(buf); free(recv_flag); free(unpack_adds);
	clearWork();
}
#endif
//end of bufComm.h
