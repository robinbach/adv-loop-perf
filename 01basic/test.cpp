/* add.cpp
* a simple C program
*/
  
#include <stdio.h>
#define LAST 10

// int addone(int a, int b)
// {
// 	return a + b + 1;
// }
int main()
{
    // int i, sum = 0;
   
    // for ( i = 1; i <= LAST; i++ ) {
    //   sum += i;
    // } /*-for-*/
    // printf("sum = %d\n", sum);

    // printf("add one = %d\n", addone(1, 2));

    // return 0;
  int in[1000]; 
  int i,j;
  FILE* myfile;

  for (i = 0; i < 1000; i++)
  {
    in[i] = 0;
  }   

  for (j = 0; j < 10; j++)
  {
   in[j]+= 10;
    for (int k = 0; k < 10; k++)
    {
      in[k] = 0;
    }   

    for (i = 0; i < 10; i++)
    {
      in[i] = 0;
    }
  }

  
  for (i = 0; i< 1000; i++)
    fprintf(stdout,"%d\n", in[i]);
  
  return 1;
}
