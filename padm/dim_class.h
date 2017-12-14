#ifndef class_dim_h
#define class_dim_h
#include <iostream>
#include "tool.h"
using namespace std;
#define dim _PA_Dim
#define data _PA_Data
#define proc _PA_Proc
#define node _PA_Node
#define ref _PA_Ref
#define array _PA_Array
const int TAG_DATA = 0;
const int TAG_PROC = 1;
const int TAG_DIMS = 2;
const int REF_DIMS = 3;
const int OUT_REFS = 4;
int _cnt_flag = true;
//-----class------
class _PA_Dim{
public:
    int tag;                //0-data,  1-proc, 2-two, 3-ref
    bool fixed;             //fixed:do not change the step
    bool is_ref;
    //int size;               //size of dim
    long long int size;               //size of dim
    int step;               //step   unused for TAG_DIMS
    int disp;				//displacement
	_PA_Dim *left, *right;  //left  & right dim
    _PA_Dim *base;   		//for ref
public:
    _PA_Dim(){ tag=TAG_DATA; size=step=0; disp=0; left=right=base=NULL; fixed=0; is_ref=0;}
    _PA_Dim(long long int s){ tag=TAG_DATA; size=s; step=1; disp=0; left=right=base=NULL; fixed=0; is_ref=0;}
    _PA_Dim(long long int s, int p){ tag=TAG_DATA; size=s; step=p; disp=0; left=right=base=NULL; fixed=1; is_ref=0;}
    _PA_Dim(long long int s, int p, int d){ tag=TAG_DATA; size=s; step=p; disp=d; left=right=base=NULL; fixed=1; is_ref=0;}
    void print(int k=0){
		//if(is_ref && tag!=TAG_DIMS){ printf("(%d,%d)", size, step); if(disp) printf("+%d", disp); printf(">>"); base->print(1);}
		if(is_ref && tag!=TAG_DIMS){ printf("(%lld,%d)_tag%d", size, step, tag); if(disp) printf("+%d", disp); printf(">>"); base->print(1);}
		else if(tag == TAG_DATA) { printf("(data %lld,%d)", size, step); if(disp) printf("+%d", disp); }   
        else if(tag == TAG_PROC){ printf("(proc %lld,%d)", size, step); if(disp) printf("+%d", disp); }   
        else if(tag == TAG_DIMS) { printf("["); left->print(1); printf(","); right->print(1); printf("]"); if(disp) printf("+%d", disp);} 
		else if(tag == OUT_REFS) { left->print(1); printf(" >> "); right->print(1); }
		if(k==0) printf("\n");
    }   
    void duplicate(_PA_Dim &d);
    void fix() {fixed=1; if(tag == TAG_DIMS){left->fix(); right->fix();} }
	void baseDim(_PA_Dim &d){
		base = new _PA_Dim; fixed =1; is_ref = 1;
		if(tag==TAG_DIMS){
            left->base = new _PA_Dim; right->base = new _PA_Dim;
			left->baseDim(d); right->baseDim(d);
            return;
        }
		if(d.tag!=TAG_DIMS) tag = d.tag; 
		else tag = REF_DIMS;
        base->duplicate(d); 
    }
    int calc_size(int chs){
		if(tag==TAG_DIMS) return left->calc_size(chs)*right->calc_size(chs);
		if(is_ref) { 
			if(size*base->size==1 || step%base->size==0) return 1;
			if(base->tag==chs){ int rg = modRange(step, base->calc_size(chs)); return min(size, rg); }//change
			if(base->tag==TAG_DIMS){
				_PA_Dim tmpr = _PA_Dim(size, step); tmpr.baseDim(*(base->right));	
				int rs = base->right->size; int rmod = modRange(step, rs); int ss = (size+rmod-1)/rmod;
				if(disp%rs!=0){
					int disp_range = size%rs; if(disp_range==0) disp_range=rs;
					disp_range = disp_range*step+disp-(step-1); if(disp_range>rs) ss++;	
				}
				int pp = (step+rs-1)/rs; 
				_PA_Dim tmpl = _PA_Dim(ss, pp); tmpl.baseDim(*(base->left)); 
				int rsize = tmpr.calc_size(chs); int lsize = tmpl.calc_size(chs); 
				return min(size, lsize*rsize);
			}
			return 1;
		}
		if(tag==chs) return size;
		return 1;
	}
	int calc_step(int chs){
		int p=0;
		if(tag==TAG_DIMS){
			p = right->calc_step(chs); if(p) return p;
			p = left->calc_step(chs); if(p) return p;
		} 
		if(is_ref) return step*base->calc_step(chs);
		if(tag==chs) return step;
		return 0;
	} 
	int calc_idx(long long int idx, int chs){
		if(tag==TAG_DIMS){
			long long int len = right->size; long long int x = idx/len; long long int y = idx%len;
			return left->calc_idx(x, chs) + right->calc_idx(y, chs) + disp;
		}
		if(tag==OUT_REFS){
			long long int t = left->calc_idx(idx, TAG_DATA);
			return right->calc_idx(t, chs);
		}
		int res = (idx+disp)*step;
		if(is_ref) return base->calc_idx(res%base->size, chs);
		if(tag==chs) return res;
		return 0;
	}
	void calc_mul(int n, int chs){
		if(n==0 || fixed) return;
		if(tag==TAG_DIMS){ left->calc_mul(n, chs); right->calc_mul(n, chs); return; }
		if(tag==chs) step *= n;
	}
	int dsize(){ return calc_size(TAG_DATA); }   
    int psize(){ return calc_size(TAG_PROC); } 
    int dstep(){ return calc_step(TAG_DATA); }
    int pstep(){ return calc_step(TAG_PROC); }
	int offset(long long int idx){ return calc_idx(idx, TAG_DATA); }
	int nodes(long long int idx){ return calc_idx(idx, TAG_PROC); }
	void print2(){
		for(int i=0; i<size; i++)
			printf("node %d, offset %d\n", nodes(i), offset(i));
	}
	void muld(int n) { return calc_mul(n, TAG_DATA); }
    void mulp(int n) { return calc_mul(n, TAG_PROC); }
    void mul(){
		if(tag!=TAG_DIMS) return;
		if(left->tag==TAG_DIMS) left->mul();
		if(right->tag==TAG_DIMS) right->mul();
		left->muld( right->dsize()*right->dstep() );
		left->mulp( right->psize()*right->pstep() );
	}
	int cntLen(int len=0){
		if(len==0){ _cnt_flag = true; len = 1;} //init flag
		if(!_cnt_flag) return len; //over
		if(tag == TAG_DIMS){ return left->cntLen( right->cntLen(len) ); }
		if(is_ref){
			if(disp){ _cnt_flag = false; return len; }//disp
			if(base->tag == TAG_DATA && len == step*base->step ){
				int bsize = base->size; 
				if(size<=bsize) return size*len;
				if(size>bsize){_cnt_flag = false; return gcd(size, bsize)*len;}
			}
			if(base->tag == TAG_DIMS){
				int brs = base->right->size; int rsize = min(size, brs); int lsize = size/rsize;
				_PA_Dim R = _PA_Dim(rsize, step); R.baseDim(*(base->right));
				_PA_Dim L = _PA_Dim(lsize, step); L.baseDim(*(base->left));
				if(size<=brs) return R.cntLen(len);//only right dim matters. 1108.
				int res = L.cntLen(R.cntLen(len));
				return gcd(size, res);
			}
			_cnt_flag = false; return len;
		}
		if(tag == TAG_DATA && step == len) return size*len;
		_cnt_flag = false; return len;
	}
	int width(){ if(tag == TAG_DIMS) return right->size;  else return 1;}
	int height(){ if(tag == TAG_DIMS) return left->size; else return 1;}
	_PA_Dim ld(){ return *(left); }
	_PA_Dim rd(){ return *(right); }
	int root(){
		if(this->psize()>1) return -1;
		if(tag==TAG_DIMS){
			int t; int res=0;
			t = left->root(); if(t>=0) res+=t;
			t = right->root(); if(t>=0) res+=t;
			return res;
		}
		if(tag==TAG_PROC) return this->nodes(0);//change
		return -1;
	}
	bool null(){ if(size==0) return true; return false; }
	void append(_PA_Dim &D){
		if(this->null()){ this->duplicate(D); return; }
		this->left = new _PA_Dim; this->left->duplicate(*this);
		this->right = new _PA_Dim; this->right->duplicate(D);
		this->step = 0; this->disp = 0; this->is_ref=0; this->tag = TAG_DIMS;
		this->size = this->left->size*D.size; 
	}
	int dims(){
		if(tag==TAG_DIMS) return left->dims()+right->dims();
		return 1;
	}
};

class _PA_Data : public _PA_Dim{
public:
	_PA_Data():_PA_Dim(1){ tag=TAG_DATA;}
	_PA_Data(long long int s):_PA_Dim(s){ tag = TAG_DATA;}
	_PA_Data(long long int s, int p):_PA_Dim(s, p){ tag = TAG_DATA;}
	_PA_Data(long long int s, int p, int d):_PA_Dim(s, p, d){ tag = TAG_DATA;}
};

class _PA_Proc : public _PA_Dim{
public:
    _PA_Proc():_PA_Dim(1){ tag=TAG_PROC;}
    _PA_Proc(int s):_PA_Dim(s){ tag=TAG_PROC;}
	_PA_Proc(int s, int p):_PA_Dim(s, p){ tag=TAG_PROC;}
	_PA_Proc(int s, int p, int d):_PA_Dim(s, p, d){ tag=TAG_PROC;}
};

class _PA_Node : public _PA_Proc{
public:
	_PA_Node():_PA_Proc(1){}
	_PA_Node(int id):_PA_Proc(1, 1, id){ fixed = 1; }
	_PA_Node(int id, int num):_PA_Proc(num){ is_ref = 1; fixed = 1; this->base = new _PA_Proc; this->base->disp = id; }
};

//for most outside reference
class _PA_Ref : public _PA_Dim{
public:
    _PA_Ref():_PA_Dim(){ is_ref=0;}
    _PA_Ref(_PA_Dim a, _PA_Dim b){
		is_ref=0; tag = OUT_REFS; size = a.size;
		left = new _PA_Dim; right = new _PA_Dim; 
		left->duplicate(a); 
		right->duplicate(b);
	}
};

class _PA_Array : public _PA_Dim{
public:
	_PA_Array():_PA_Dim(){ is_ref=1; fixed=1;}
	_PA_Array(_PA_Dim x){
		this->duplicate(x);
		if( !fixed && tag==TAG_DIMS) this->mul();
		this->fix();
	}
	_PA_Array(_PA_Dim x, _PA_Dim y){
		_PA_Dim a = _PA_Dim(); a.duplicate(x);
		if(!a.fixed && a.tag==TAG_DIMS){
			if(y.tag!=TAG_DIMS) a.mul();
			else{ a.left->mul(); a.right->mul(); }
		}
		this->duplicate(a); is_ref=1;
		_PA_Dim b = _PA_Array(y);
		if(a.tag==TAG_DIMS && b.tag==TAG_DIMS){ is_ref=0; left->baseDim(*b.left); right->baseDim(*b.right);}//debug //is_ref=0
		else this->baseDim(b);
	}
};

//class function
void _PA_Dim::duplicate(_PA_Dim &d){
    tag = d.tag; size = d.size; step = d.step; disp = d.disp;
    fixed = d.fixed; is_ref = d.is_ref;
	if(tag == TAG_DIMS){
    	right = new _PA_Dim; left = new _PA_Dim;
    	right->duplicate(*(d.right)); left->duplicate(*(d.left));
		return;
	}
    if(is_ref){ base = new _PA_Dim; base->duplicate(*d.base); }
	return;
}

//overload operator---------------------------------------------------------------------------------
_PA_Dim operator*(_PA_Dim a, _PA_Dim b){
    _PA_Dim c = _PA_Dim(a.size*b.size); c.tag = TAG_DIMS;
	c.left = new _PA_Dim; c.left->duplicate(a);
    c.right = new _PA_Dim; c.right->duplicate(b);
    return c;
}

_PA_Dim operator+(_PA_Dim a, int d){
	a.disp = d; return a;
}
#endif
//end of dim.h
