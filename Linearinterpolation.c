Enter file contents here
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;
int main(void)
{
	/*signal 생성*/
	int L = 8;     // 신호길이
	int xqcol_vals = 32;     // 신호길이
	double pi = 3.14;

	double *x = (double*)malloc(sizeof(double) * L); // 메모리 동적 할당
	for (int i = 0; i < L; i++)
	{
		x[i] = (double)(i + 1) * pi / 4; 
		//printf("%lf\n", x[i]);
	}
	
	double *v = (double*)malloc(sizeof(double) * L);
	int count1 = 0;
	for (int i = 0; i < L; i++)
	{
		v[i] = (double)sin(x[i]);
		//printf("%lf\n", v[i]);
		count1++;
	}
	//printf("%d\n", count1);

	double *xq = (double*)malloc(sizeof(double) * xqcol_vals); // 메모리 동적 할당
	for (int i = 0; i < xqcol_vals; i++)
	{
		xq[i] = (double)(i + 1) * pi / 16; 
		//printf("%lf\n", xq[i]);
	}

	int k;
	k = -1;
	double *interpout = (double*)malloc(sizeof(double) * xqcol_vals);
	for (int i = 0; i < L-1; i++) 
	{
		for (int j = 1+k; j < xqcol_vals; j++)
		{
			if (xq[j] <= x[1 + i]) 
			{
				interpout[j] = (v[1 + i] - v[i]) / (x[1 + i] - x[i]) * (xq[j] - x[i]) +	v[i];
				k = j;
			}
		}
	}

	return 0;
}
