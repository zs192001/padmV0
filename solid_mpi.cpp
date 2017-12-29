#include<stdlib.h> 
#include<stdio.h>
#include<string.h> 
#include<time.h> 
#include<math.h>  //-lm
#include <mpi.h>
#include "padm/padm.h"
#define M_PI 3.14159265358979323846
#define W 9.24
#define H 18.48
#define L 18.48
#define D 1
#define sogalat (D/2.0)
#define soW ((int)(W / sogalat)+1)
#define soH ((int)(H / sogalat)+1)
#define soL ((int)(L / sogalat)+1)
#define extnum 3
#define ncells ((2 * extnum + 1)*(2 * extnum + 1)*(2 * extnum + 1) - 1) //7*7*7-1
#define DT 0.01
#define NASCPU 32
#define NP 4
#define NPROCS (NP*NP*NP)
#define NAS (NPROCS*NASCPU)
#define PID(px,py,pz) ((px)+(py)*NP+(pz)*NP*NP)
typedef double Vec[3];
const int fnx[26] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
const int fny[26] = { -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1 };
const int fnz[26] = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1 };
int nbCellx[ncells], nbCelly[ncells], nbCellz[ncells];
int soa[NAS]; int soaR[NAS]; int  somap[soH*soW*soL] = { 0 }; int latNumx, latNumy, latNumz; double vs, soma, resoma, lttcL;
double AA = 1.8308e3, BB = 4.7118e2, Mu = 1.7322, Lambda = 2.4799, cc = 1.0039e5, dd = 1.6218e1, hh = 0,//-5.9826e-1, 
Beta = 1.0999e-6, nexp = 7.8734e-1, IntCut = 2.7, ExtCut = 3.0, scalefactor = 2.35, nexpisqr = -1.0 / (2.0*7.8734e-1);
typedef struct { Vec s, v, a; } State;
State state[NAS]; 
State atoms[NAS]; 
State steps[1000][27*NASCPU]; 
static unsigned long next = 1;
double rnd() { next = next * 1103515245 + 12345; return (((unsigned)(next / 65536) % 32768)) / (double)(32768); }
inline double loop(double f, double m) { return fabs(f)<m/2.0 ? f : (f>0 ? f - m : f + m); }
inline void vzero(Vec a) { a[0] = a[1] = a[2] = 0; }
inline void vxv(Vec a, Vec b) { a[0] *= b[0];  a[1] *= b[1]; a[2] *= b[2]; }
inline void vxd(Vec a, double d, Vec b) { b[0] = d*a[0];  b[1] = d*a[1]; b[2] = d*a[2]; }
inline void vav(double d, Vec b, Vec a) { a[0] += d*b[0];  a[1] += d*b[1]; a[2] += d*b[2]; }
inline void vmav(double d, Vec c, Vec b, Vec a) { c[0] = a[0]+d*b[0];  c[1] = a[1]+d*b[1]; c[2] = a[2]+d*b[2]; }
inline void vsv(Vec b, Vec a) { a[0] -= b[0];  a[1] -= b[1]; a[2] -= b[2]; }
inline double dot_prod(Vec a, Vec b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline void vec_add(double q, Vec a, double p, Vec b, Vec c) { c[0] = q*a[0] + p*b[0];  c[1] = q*a[1] + p*b[1]; c[2] = q*a[2] + p*b[2]; }
void re_order(int px, int py, int pz, State* stt, State* ord) {
	int i, j=0; int p=PID(px,py,pz); double _X = 2.31, _Y = 4.62;
    for(i=0;i<NAS;i++) 
    if(stt[i].s[0]>px*_X+0.1 && stt[i].s[0]<px*_X+_X+0.1) 
    if(stt[i].s[1]>py*_Y-0.1 && stt[i].s[1]<py*_Y+_Y-0.1) 
    if(stt[i].s[2]>pz*_Y-0.1 && stt[i].s[2]<pz*_Y+_Y-0.1) {
        int _id = p*NASCPU+j;
		memcpy(ord[_id].s, stt[i].s, sizeof(Vec)); //here!
		memcpy(ord[_id].v, stt[i].v, sizeof(Vec)); //here!
		j++;
    } 
	if(j!=NASCPU) exit(-1);
}
void reOrder(State* stt, State* ord){
	for(int i=0; i<NP; i++) for(int j=0; j<NP; j++) for(int k=0; k<NP; k++)
		re_order(i, j, k, stt, ord);
}
void updateS(double dt, Vec ss, Vec nsv, Vec nss) {
	vmav(dt, nss, nsv, ss);
	if (nss[0] >= W) nss[0] -= W; else if (nss[0] < 0) nss[0] += W;
	if (nss[1] >= H) nss[1] -= H; else if (nss[1] < 0) nss[1] += H;
	if (nss[2] >= L) nss[2] -= L; else if (nss[2] < 0) nss[2] += L;
}
void init(double dt, State stt[NAS]) {
	vs = 0.279;  soma = 1; resoma = 1.0 / soma; lttcL = 2.31*D;
	latNumx = (int)(W / lttcL + 0.001); latNumy = (int)(H / lttcL + 0.001); latNumz = (int)(L / lttcL + 0.001);
	AA *= 1.0; BB *= 1.0; Mu *= (double)(scalefactor) / D; Lambda *= (double)(scalefactor) / D; 
	IntCut *= D / (double)(scalefactor); ExtCut *= D / (double)(scalefactor);
	int i = 0, j, k, m; for (j = -extnum; j <= extnum; j++) for (k = -extnum; k <= extnum; k++) for (m = -extnum; m <= extnum; m++)
	if (j != 0 || k != 0 || m != 0) { nbCellx[i] = j; nbCelly[i] = k; nbCellz[i] = m; i++; }
	int nn, p = 0; Vec ave = { 0 }; srand(872877);
	double StatePos[8][3] = { 	{ lttcL/4.0, lttcL/4.0, 3.*lttcL/ 4.0 }, { lttcL/4.0, 3.*lttcL/4.0, lttcL/4.0 }, 
								{ lttcL/2.0, 0, lttcL/2.0 },{ lttcL/2.0, lttcL / 2.0, 0 }, 
								{ 3.*lttcL / 4.0, lttcL / 4.0, lttcL / 4.0 }, { 3.*lttcL/4.0, 3.*lttcL / 4.0, 3.*lttcL / 4.0 }, 
								{ lttcL, 0, 0 }, { lttcL, lttcL / 2.0, lttcL / 2.0 } };
	Vec sp[NAS], sv[NAS];  
	for (j = 0; j<latNumx; j++) for (nn = 0; nn<4; nn++) for (i = 0; i<latNumy; i++) for (k = 0; k<latNumz; k++) for (m = 0; m<2; m++) {
		sp[p][0] = StatePos[2*nn+m][0] + j*lttcL; sp[p][1] = StatePos[2*nn+m][1] + i*lttcL; sp[p][2] = StatePos[2*nn+m][2] + k*lttcL;
		double ang1 = rnd()*M_PI * 2, ang2 = rnd()*M_PI;
		sv[p][0] = vs*sin(ang1)*cos(ang2); sv[p][1] = vs*sin(ang1)*sin(ang2); sv[p][2] = vs*cos(ang1); vav(1.0, sv[p], ave);  p++;
	}
	for (i=0; i<NAS; i++) { 
		memcpy(stt[i].s, sp[i], sizeof(stt[i].s)); 
		vec_add(1.0, sv[i], -1 / (double)NAS, ave, stt[i].v); 
	}
	for(i=0; i<NAS; i++) vav(0.5*dt, stt[i].a, stt[i].v); //stt.a is zero at this situation, useless?
}
inline double fc(double r) { if (r<IntCut) return 1.0; else return 0.5 *(1.0 + cos(M_PI * (r - IntCut) / (ExtCut - IntCut))); }
inline double fc1(double r) { if (r<IntCut) return 0.0; else return -0.5*M_PI / (ExtCut - IntCut)*sin(M_PI *(r - IntCut) / (ExtCut - IntCut)); }
inline void cos_arc(double c, Vec a, Vec b, double inv, Vec dcosj) { dcosj[0] = -(c*a[0] - b[0])*inv, dcosj[1] = -(c*a[1] - b[1])*inv; dcosj[2] = -(c*a[2] - b[2])*inv; }
double _distV(Vec v, Vec si, Vec sj){ 
	v[0] = loop(sj[0] - si[0], W), v[1] = loop(sj[1] - si[1], H), v[2] = loop(sj[2] - si[2], L);    
	return sqrt(dot_prod(v, v)); 
}
double distV(Vec si, Vec sj){
	Vec v; return _distV(v, si, sj);
}
int nblist[NAS][40], nbnum[NAS] = { 0 };
void make_nblist(State* stt, int rnas){
	int i, j; 
	for(i = 0; i < rnas; i++){
		somap[soaR[i]] = 0; vzero(stt[i].a);
		soaR[i] = (int)(fmod(stt[i].s[2]+L,L)/sogalat)*soW*soH + (int)(stt[i].s[0]/sogalat)*soH + (int)(stt[i].s[1]/sogalat);
		somap[soaR[i]] = i + 1;
	}
	for(i=0; i<rnas; i++) for(j = 0, nbnum[i] = 0; j<ncells; j++){
		int z = soaR[i]/soH/soW, x = (soaR[i]-z*soH*soW)/soH, y = soaR[i]-z*soH*soW-x*soH;
		int pj = somap[(z+nbCellz[j]+soL)%soL*soH*soW + (x+nbCellx[j]+soW)%soW*soH + (y+nbCelly[j]+soH)%soH]; if(pj==0) continue;
		pj--; 
		if (distV(stt[i].s, stt[pj].s) <= ExtCut) 
			nblist[i][nbnum[i]++] = pj;
	}
}
double computeZ(int i, int k, double invdij, Vec reij, Vec Zj, Vec Zkk, State nstt[NAS]) {
	Vec reik, dcosj, dcosk; 
	double dik = _distV(reik, nstt[i].s, nstt[nblist[i][k]].s);	
	double invdik = 1 / dik; vxd(reik, invdik, reik);
	double ccos = dot_prod(reij, reik); double GtermDenom = dd*dd + (hh - ccos)*(hh - ccos);
	double Gterm = 1 + cc*cc / dd / dd - cc*cc / GtermDenom; double GtermD = -2.0 * cc*cc * (hh - ccos) / (GtermDenom*GtermDenom);
	cos_arc(ccos, reij, reik, invdij, dcosj); cos_arc(ccos, reik, reij, invdik, dcosk);
	double fc1dik = fc1(dik)*Gterm, fcdik0 = fc(dik), fcdik = fcdik0*GtermD;
	vec_add(fc1dik, reik, fcdik, dcosk, Zkk); vav(fcdik, dcosj, Zj);
	return fcdik0*Gterm;// Z0 no way to be 0 for pow and / 
}
void computeTSF(int nbnum, double Z0, double dij, double &TF, double &SF){
	double BetaPowN = pow(Beta*Z0, nexp); double Bondfa = -BB*exp(-Mu*dij), Bondfr = AA*exp(-Lambda*dij); if (nbnum>1) Bondfa *= pow(1 + BetaPowN, nexpisqr);
	double fcdij = fc(dij); TF = -0.5*fcdij*(Lambda*Bondfr + Bondfa*Mu); if (dij >= IntCut) TF += 0.5*fc1(dij)*(Bondfr + Bondfa);
	SF = (nbnum <= 1 ? 0 : -0.25*fcdij*Bondfa / (1 + BetaPowN)*BetaPowN / Z0);
}

void compute_step(int i, double dt, State stt[NAS], State nstt[NAS]) {
	vzero(nstt[i].a); 
	for(int j = 0; j < nbnum[i]; j++) {
		int pj = nblist[i][j]; Vec reij; 
		double dij = _distV(reij, stt[i].s, stt[pj].s); 
		double invdij = 1 / dij; vxd(reij, invdij, reij);
		double Z0 = 0; Vec Zj = { 0 }, f = { 0 }, Zk[100] = { 0 };
		for(int k = 0; k < nbnum[i]; k++) if (j != k) Z0 += computeZ(i, k, invdij, reij, Zj, Zk[k], stt);
		double TF = 0, SF = 0; computeTSF(nbnum[i], Z0, dij, TF, SF); vec_add(TF, reij, SF, Zj, f); vav(1.0, f, nstt[i].a);
		for(int k = 0; k < nbnum[i]; k++) if (j != k) { vxd(Zk[k], SF, f); vav(1.0, f, nstt[i].a); } //compute_host
		for(int nbj = 0; nbj < nbnum[pj]; nbj++){
			int nbofj = nblist[pj][nbj]; Vec rejnbj; 
			double djnbj = _distV(rejnbj, stt[pj].s, stt[nbofj].s);	
			double invdjnbj = 1 / djnbj; vxd(rejnbj, invdjnbj, rejnbj);
			double Zjk0 = 0; Vec Znbj = { 0 }, fji = { 0 }, Zjk[100] = { 0 };
			for(int k = 0; k < nbnum[pj]; k++)	if (nbj != k) Zjk0 += computeZ(pj, k, invdjnbj, rejnbj, Znbj, Zjk[k], stt); //only use nstt.s
			double TFj = 0, SFj = 0; computeTSF(nbnum[pj], Zjk0, djnbj, TFj, SFj);
			if (nbofj == i){ vec_add(TFj, rejnbj, SFj, Znbj, fji); vsv(fji, nstt[i].a); } //compute_guest
			else{ for(int k = 0; k < nbnum[pj]; k++) if (i == nblist[pj][k]) { vxd(Zjk[k], SFj, fji); vsv(fji, nstt[i].a); } } //compute_observer
		}
	}
	vmav(dt, nstt[i].v, nstt[i].a, stt[i].v);
}
void compute_steps(int i, int t, double dt, State stt[NAS], State nstt[NAS]) {
	compute_step(i, dt, stt, nstt);
	updateS(dt, stt[i].s, nstt[i].v, nstt[i].s); //stt.s, nstt.v -> nstt.s 
}
//dimension
dim _B = data(27, 1);	
dim _Px = proc(NP,1);
dim _Py = proc(NP,NP);
dim _Pz = proc(NP,NP*NP);
dim _P = _Pz * _Py * _Px;
dim A1 = _P * get_low(cyclic(_B,13),1);
dim Blt = _Pz * _Py * cyclic(_Px, 1) * get_low(cyclic(_B,12),1);
dim Brt = _Pz * _Py * cyclic(_Px,-1) * get_low(cyclic(_B,14),1);
dim A3 = _P * get_low(cyclic(_B,12),3);
dim Bup = _Pz * cyclic(_Py, 1) * _Px * get_low(cyclic(_B,9),3);
dim Bdn = _Pz * cyclic(_Py,-1) * _Px * get_low(cyclic(_B,15),3);
dim A9 = _P * get_low(cyclic(_B,9),9);
dim Bfr = cyclic(_Pz, 1) * _Py * _Px * get_low(_B,9);
dim Bbk = cyclic(_Pz,-1) * _Py * _Px * get_low(cyclic(_B,18),9);
void exchange(int t){
	State* a = steps[t];
	//x dimension 	
	copyto(A1, a, Blt, a, sizeof(State)*NASCPU);
	copyto(A1, a, Brt, a, sizeof(State)*NASCPU);
	//y dimension
	copyto(A3, a, Bup, a, sizeof(State)*NASCPU);
	copyto(A3, a, Bdn, a, sizeof(State)*NASCPU);
	//z dimension
	copyto(A9, a, Bfr, a, sizeof(State)*NASCPU);
	copyto(A9, a, Bbk, a, sizeof(State)*NASCPU);
}
void copy27(State* stt, State* cubes, int x, int y, int z){
	dim _rX = get_low(cyclic(data(NP, 1), x-1), 3);
	dim _rY = get_low(cyclic(data(NP, NP), y-1), 3);	
	dim _rZ = get_low(cyclic(data(NP, NP*NP), z-1), 3);
	dim _A = _rZ*_rY*_rX;
	copyto(_A, stt, _B, cubes, NASCPU*sizeof(State));
}
int main() {
	MPI_Init(NULL, NULL); int rank; MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int pidx = rank%NP; int pidy = rank/NP%NP; int pidz = rank/NP/NP%NP;
	State *stt = state, *nstt; int i0 = 13*NASCPU;
	init(DT, stt); reOrder(stt, atoms); 
	copy27(atoms, steps[0], pidx, pidy, pidz); 
	make_nblist(steps[0], 27*NASCPU);
	
	//compute one more, update steps[0].s
	stt = steps[0];
	for(int i=i0; i<i0+NASCPU; i++) updateS(DT, stt[i].s, stt[i].v, stt[i].s);
	if(rank==0) printf("t = 0, test0_xx=%.16e\n", steps[0][i0].s[0]);
	
	for(int t=0; t<63; t++){
		stt = steps[t]; nstt = steps[t+1]; 
		MPI_Barrier(MPI_COMM_WORLD); exchange(t);
		for (int i = i0; i < i0+NASCPU; i++) compute_steps(i, t, DT, stt, nstt); 
		if(rank==0) printf("t = %d, test0_xx=%.16e\n", t, nstt[i0].s[0]);
	}
	if(rank==0) printf("sequence version: 5.9511132573079328\n");	
	MPI_Finalize(); return 0;
}
/*process note:
 *1, init: stt.s; stt.s->stt.v; stt.s+stt.v-> stt'.s
 *2, exchange(t): stt( or only stt'.s)
 *3, compute_steps(stt, nstt):	compute: stt'.s -> nstt.a; stt.v(own)+nstt.a -> nstt.v 
 *								updateS: stt'.s+nstt.v -> nstt.s
 */



