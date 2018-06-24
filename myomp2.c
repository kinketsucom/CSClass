#include<stdio.h>
#include<omp.h>

#define L 1000LL
#define R 2000000LL

int main(){
  long long n, i, j, t;
  long long answer;
  double t1, t2;

  t1 = omp_get_wtime();

  answer = 0;
#pragma omp parallel
{ 
 #pragma omp for private(j,t) reduction(+:answer) schedule(dynamic,1)

  for(i=L; i<=R; i++){
  t = i;
  for(j=2; j*j<=i; j++){
     if(i%j==0){
        t = j;
        break;
      }
    }
    answer += t;
  }
}

  t2 = omp_get_wtime();

  printf("answer %lld\n", answer);
  printf("time %f\n", t2-t1);
  
  return 0;
}
