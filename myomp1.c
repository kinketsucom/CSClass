#include<stdio.h>
#include<omp.h>

#define N 20000
#define K 20000
#define P 1001001011

int S[N][K];

int main(){
  int n, k;
  double t1, t2;

  t1 = omp_get_wtime();
  S[0][0] = 1;
 for(n=1;n<N;n++){
#pragma omp parallel for  
  for(k=1;k<K;k++){
      S[n][k] = (S[n-1][k-1] + (long long)(n-1)*S[n-1][k]) % P;	
	}
  }
  t2 = omp_get_wtime();

  printf("%d %d %d\n", S[N-1][K/4], S[N-1][K/3], S[N-1][K/2]);
  printf("time %f\n", t2-t1);
  
 return 0;
}
