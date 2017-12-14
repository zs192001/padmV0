#include <math.h>
//greatest common divisor
int gcd(int a, int b){
	int x, y, t;
	if(a>b){ x=a; y=b;}
	else{ x=b; y=a;}
	while(y!=0){
		t = x%y; x=y; y=t;
	}
	return x;	
}

int log2(int x){
	return log(x)/log(2);
}
long long int log2(long long int x){
	return log(x)/log(2);
}
int max(int x, int y){
	return x>y ? x : y;
}
int min(int x, int y){
	return x>y ? y : x;
}
int modRange(int val, int mod){
	if(val%2 || mod%2) return mod;
	int i = val; int count =1;
	while(i%mod){ count++; i+=val; }
	return count;
}
//end of tool.h ..
