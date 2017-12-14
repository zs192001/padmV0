#ifndef args_h
#define args_h
#define PER_MAX_GRAN (1024*16)
#define PER_MIN_GRAN (1024*4)

long long int USR_GRAN_BYTES = 0; //irecv: int count;
long long int USR_BUF_BYTES = 0;

int _tnum = 4;
int _procs = 1;
int _reqs;
long long int _scope;
int _blks;

long long int _len; //cnt_len
int _gran;

long long int _buf_bytes;
int _elem_size;

int _prt_flag = 0;
int CHK_FLAG = 0;
#endif
