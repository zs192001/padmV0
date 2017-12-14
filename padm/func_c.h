#ifndef func_c_h
#define func_c_h
////----------------------------------------------
#include <omp.h>
#include <sched.h>
#include <queue>
#include <vector>
#include <algorithm>
#include <semaphore.h>
#include "args.h"
using namespace std;
//
struct msg_info{int dest; int tag; int len; int tid; int blk; void* buf;};
//
omp_lock_t qlock;
//
queue<void*> buf_add;
queue<int> sbuf_idx;
queue<int> rbuf_idx;
queue<msg_info> send_msgs;
queue<msg_info> recv_msgs;
//
sem_t * sem_recv;
int packs_count;
bool pack_over;
int** check_recv;
int* post_cur;
//
void** list;
int list_idx;
//_______________________________________________________________queue for requests & buffer index
msg_info popQue(queue<msg_info> &q){
    msg_info t; omp_set_lock(&qlock);
    t = q.front(); q.pop(); omp_unset_lock(&qlock);
    return t;
}
int popQue(queue<int> &q){
    int t; omp_set_lock(&qlock);
    t = q.front(); q.pop(); omp_unset_lock(&qlock); 
	return t;
}
void pushQue(queue<msg_info> &q, msg_info msg){
    omp_set_lock(&qlock);
    q.push(msg); omp_unset_lock(&qlock);
    return;
}
void pushQue(queue<int> &q, int idx){
    omp_set_lock(&qlock);
    q.push(idx); omp_unset_lock(&qlock);
    return;
}
void initQue(queue<int> &q, int num){
    while(!q.empty()){ q.pop(); }
    for(int i=0; i<num; i++){ q.push(i); }
    return;
}
void regInfo(msg_info &msg, int dest, int tag, int len, int tid, int blk){
    msg.dest = dest; msg.tag = tag; msg.len = len;
	msg.tid = tid; msg.blk = blk;
    return;
}
void regInfo(msg_info &msg, int dest, int tag, int len, int idx){
    msg.dest = dest; msg.tag = tag; msg.len = len; 
	msg.blk = idx;
    return;
}
void regInfo(msg_info &msg, int dest, int tag, void* buf){
    msg.dest=dest; msg.tag=tag; msg.buf=buf;
    return;
}
void getInfo(msg_info msg, int &dest, int &tag, int &len, int &tid, int &blk){
    dest = msg.dest; tag = msg.tag; len = msg.len;
	tid = msg.tid; blk = msg.blk;
    return;
}
void getInfo(msg_info msg, int &dest, int &tag, int &len, int &idx){
    dest = msg.dest; tag = msg.tag; len = msg.len; 
	idx = msg.blk;
    return;
}
void getInfo(msg_info msg, int &dest, int &tag, void **buf){
    dest = msg.dest; tag = msg.tag; *buf = msg.buf;
    return;
}

//---------------------------------to _msgs
bool sendMsgEmpty(){
	bool t;
    omp_set_lock(&qlock);
	t = send_msgs.empty(); 
	omp_unset_lock(&qlock);
	return t;
}
bool recvMsgEmpty(){
	bool t;
    omp_set_lock(&qlock);
	t = recv_msgs.empty(); 
	omp_unset_lock(&qlock);
	return t;
}
void regSendMsg(int dest, int tag, int len, void *buf){
    msg_info msg; 
	msg.dest = dest; msg.tag = tag; msg.len = len; msg.buf = buf;
	pushQue(send_msgs, msg);
	return;
}
void regRecvMsg(int dest, int tag, int len, void *buf, int tid, int blk){
    msg_info msg; 
	msg.dest = dest; msg.tag = tag; msg.len = len; msg.buf = buf;
	msg.tid = tid; msg.blk = blk;
    pushQue(recv_msgs, msg);
	return;
}
void getSendMsg(int &dest, int &tag, int &len, void **buf){
    msg_info msg = popQue(send_msgs);
	dest = msg.dest; tag = msg.tag; len = msg.len; *buf = msg.buf;
    return;
}
void getRecvMsg(int &dest, int &tag, int &len, void **buf, int &tid, int &blk){
    msg_info msg = popQue(recv_msgs);
	dest = msg.dest; tag = msg.tag; len = msg.len; *buf = msg.buf;
	tid = msg.tid; blk = msg.blk;
    return;
}
//_______________________________________________________________help functions
int initBufAdd(void *buf){
	int count=0;
	int size = _buf_bytes/_elem_size;
	while(!buf_add.empty()) buf_add.pop();//maybe usless
	for(int i=0; i+_gran<=size; i+=_gran){
		buf_add.push(buf+i*_elem_size);
		count++;
	}
	return count;	
}
void* getBufAdd(int bytes){
    void* t = NULL; 
	omp_set_lock(&qlock);
	if(!buf_add.empty()){
    	t = buf_add.front(); buf_add.pop(); 
	}
	else{
		t = list[list_idx++] = (void*)malloc(bytes);
	}
	omp_unset_lock(&qlock); 
	return t;
}
void pushBufAdd(void* ptr){
	omp_set_lock(&qlock);
	buf_add.push(ptr);
	omp_unset_lock(&qlock);
	return;
}
void initList(int num){
	list = (void**)malloc(sizeof(void*)*num);
	for(int i=0; i<num; i++) list[i]=NULL;
	return;
}
void clearList(){
	for(int i=0; i<list_idx; i++){
		if(list[i]!=NULL) printf("a new gran blk %p, \n", list[i]);
		free(list[i]); list[i] = NULL;
	}
}
void semInit(int num){
    if(sem_recv) sem_recv = (sem_t*)realloc(sem_recv, num*sizeof(sem_t)); 
    else sem_recv = (sem_t*)malloc(num*sizeof(sem_t));
    for(int i=0; i<num; i++) sem_init(&sem_recv[i], 0, 0);
    return;
}
void checkMalloc(int tnum, int blk){//malloc check_recv, post_cur
    if(check_recv){
		check_recv = (int**)realloc(check_recv, tnum*sizeof(int*));
		for(int i=0; i<tnum; i++) check_recv[i] = (int*)malloc(blk*sizeof(int));		
	} 
	else{
		check_recv = (int**)malloc(_tnum*sizeof(int*));
		for(int i=0; i<tnum; i++) check_recv[i] = (int*)malloc(blk*sizeof(int));
	}
	if(post_cur) post_cur = (int*)realloc(post_cur, tnum*sizeof(int));
	else post_cur = (int*)malloc(tnum*sizeof(int));
}

void initRecvPost(){
	for(int i=0; i<_tnum; i++){
		for(int j=0; j<_blks; j++) check_recv[i][j] = -1;
		post_cur[i] = 0; 
	}
}
void setInitVal(){
    initQue(sbuf_idx, _reqs);
    initQue(rbuf_idx, _reqs);
    initRecvPost();//check_recv and post_cur
	packs_count = 0;
    pack_over = false;
    return;
}
void launchWork(){
    omp_init_lock(&qlock);////
	semInit(_tnum);
    checkMalloc(_tnum, _blks);
	setInitVal();
}
void clearWork(){
    for(int i=0; i<_tnum; i++) sem_destroy(&sem_recv[i]);
	clearList();
}
//check_recv related
void* getUnpackBuf(int tid, int cur, bool* flag, void **adds){
	int idx = tid*_blks + cur;
	sem_wait(&sem_recv[tid]);
	while(flag[idx]){ 
		return adds[idx]; 
	}
	return NULL;	
}
void markAndPost(int tid, int cur, void *ptr, bool *flag, void **adds){
	int idx = tid*_blks + cur;
	flag[idx] = true; 
	adds[idx] = ptr;
	if(cur == post_cur[tid])
		for(int i=idx; i<(tid+1)*_blks; i++){
			if(adds[i]!=NULL){ 
				sem_post(&sem_recv[tid]); 
				post_cur[tid]++; 
			}
			else return;
		}
	return;			
}
void markAndPost(int tid, int cur, int idx){
	check_recv[tid][cur] = idx;
	if(cur == post_cur[tid]) 
		for(int i=cur; i<_blks; i++){
			if(check_recv[tid][i]!=-1){ 
				sem_post(&sem_recv[tid]); 
				post_cur[tid]++; 
			}
			else return;
		}
	return;			
}
//args functions
//set args by user
void setTnum(int tnum){ _tnum = tnum; }
void setGranBytes(int gran){ USR_GRAN_BYTES = gran; }
void setBufBytes(long long int bytes){ USR_BUF_BYTES = bytes; }
void setPrint(int flag){ _prt_flag = flag; }
//calculate args 
int getProcs(int ap, int bp){ return max(ap, bp); }
int getLen(int acnt, int bcnt){ return gcd(acnt, bcnt); }
int maxLocal(int s, int d, int p){
	int pnum = s/d;
	if(p<pnum) return ((pnum+p-1)/p)*d;
   	else return d;	
}
int maxPack(int ad, int bd, int p, int t){ return max(ad, bd)/(p*t);}
long long int getScope(long long int s, int g, int t){ return (((s+t-1)/t + g-1)/g)*g;}
int getBlks(long long int s, int g, int p){ return (s+g-1)/g + p; }
int getReq(int s, int g, int p, int t){ return (s+g-1)/g + p*t; }
int getGran(int cnt, int pack_max,int e_bytes){
	if(USR_GRAN_BYTES!=0) return (USR_GRAN_BYTES/e_bytes);
	int pack = (max(cnt, pack_max)/cnt)*cnt; int max=PER_MAX_GRAN, min=PER_MIN_GRAN;
	if(cnt<max) max = (max/cnt)*cnt; else max = cnt; //else if(cnt%max!=0) max = cnt;
	if(cnt<min) min = (min/cnt)*cnt; else min = cnt; //else if(cnt%min!=0) min = cnt;
	if(pack < min) return pack; 
	if(pack < max) max = pack;
	if(cnt>=max) return max;
	if(cnt>=min) return cnt;
	return min;
}
long long int getBytes(int reqs, int gran, int elem_bytes){
	if(USR_BUF_BYTES!=0) return USR_BUF_BYTES;
	long long int res = reqs*2;
	res *= gran; res *= elem_bytes;
	return res; 
}
void setArgs(_PA_Dim A, _PA_Dim B, int elem_size, int pnum=1){
	int a_dsize = A.dsize(); int b_dsize = B.dsize();
	int a_psize = A.psize(); int b_psize = B.psize();
	_len = getLen(A.cntLen(), B.cntLen());
	_procs = pnum;
	int max_pack = maxPack(a_dsize, b_dsize, getProcs(a_psize, b_psize), _tnum);
	_gran = getGran(_len, max_pack, elem_size);
	_scope = getScope(A.size, _gran, _tnum);
    _blks = getBlks(_scope, _gran, _procs);
	if(a_dsize == b_dsize && a_psize != b_psize);
   	int al = maxLocal(A.size, a_dsize, a_psize);
	int bl = maxLocal(B.size, b_dsize, b_psize);	
    _reqs = getReq(max(al, bl), _gran, _procs, _tnum); 
	_buf_bytes = getBytes(_reqs, _gran, elem_size);
	_elem_size = elem_size;
}
void prtArgs(){
	printf("procs %d, tnum %d; ", _procs, _tnum);
	printf("requests num  %d, ", _reqs);
	if(USR_GRAN_BYTES) printf("usr gran bytes 2^%d; ", log2(USR_GRAN_BYTES) );
	printf("cnt len "); if(_len&(_len-1)) printf("%d; ", _len); else printf("2^%d; ", log2(_len));
	printf("gran "); if(_gran&(_gran-1)) printf("%d; ", _gran); else printf("2^%d; ", log2(_gran));
	printf("scope %lld, ", _scope);
	printf("\n");
}
//end of func_c.h
#endif
