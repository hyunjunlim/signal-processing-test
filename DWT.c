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
	int L = 512;     // 신호길이
	int	Fs = 512; // sample rate
	double pi = 3.14;

	double *time = (double*)malloc(sizeof(double) * L); // 메모리 동적 할당
	for (int i = 0; i < L; i++)
	{
		time[i] = (double) (i+1)/Fs; // 시간
		//printf("%lf\n", time[i]);
	}
	double *signal = (double*)malloc(sizeof(double) * L);
	int count1 = 0;
	for (int i = 0; i < L; i++)
	{
		signal[i] = (double) cos(2 * pi * 200 * time[i]) + cos(2 * pi * 100 * time[i]);
		//printf("%lf\n", signal[i]);
		count1++;
	}
	//printf("%d\n", count1);

	static const double lowd[8] = { -0.0105974017849973, 0.0328830116669829, 0.030841381835987, -0.187034811718881, -0.0279837694169839, 0.63088076792959, 0.714846570552542, 0.230377813308855 };
	static const double highd[8] = { -0.230377813308855, 0.714846570552542, -0.63088076792959, -0.0279837694169839, 0.187034811718881, 0.030841381835987, -0.0328830116669829, -0.0105974017849973 };

	int shift = 0;

	// Compute sizes and shape.
	int lf = 8;
	int lx = 512; 

	/*Signal Extensions: Dealing with Border Distortion */
	// Extend, Decompose &  Extract coefficients.
	int first1 = 2 - shift;
	int lenEXT = lf - 1;
	int last1 = lx + lf - 1;
	int *lenEXTarray = (int*)malloc(sizeof(int) * 7);
	int *lxarray = (int*)malloc(sizeof(int) * 512);
	int *lxlenEXTarray = (int*)malloc(sizeof(int) * 7);

	for (int i = 0; i < lenEXT; i++)
	{
		lenEXTarray[i] = lenEXT-i-1;
		//printf("%d\n", lenEXTarray[i]);
	}
	
	for (int i = 0; i < L; i++)
	{
		lxarray[i] = i;
		//printf("%d\n", lxarray[i]);
	}
	
	for (int i = 0; i < lenEXT; i++)
	{
		lxlenEXTarray[i] = lx-i-1;
		//printf("%d\n", lxlenEXTarray[i]);
	}
	
	int *l = (int*)malloc(sizeof(int) * 526);
	int count2 = 0;
	for (int i = 0; i < 526; i++)
	{
		if(i< lenEXT)
			l[i] = lenEXTarray[i];
		else if(lenEXT<=i && i<lenEXT +L)
			l[i] = lxarray[i- lenEXT];
		else if(lenEXT + L<=i)
			l[i] = lxlenEXTarray[i- lenEXT - L];
		count2++;
		//printf("%d\n", l[i]);
	}
	//printf("%d\n", count2);
	double *signal_extention = (double*)malloc(sizeof(double) * 526);
	int count3 = 0;
	for (int i = 0; i < 526; i++)
	{
		signal_extention[i] = signal[l[i]];

		//printf("%lf\n", signal_extention[i]);
		count3++;
	}
	//printf("%d\n", count3);

	/* Compute coefficients of approximation. 1-D Convolution. downsampling */
	double *signal_lowdconv = (double*)malloc(sizeof(double) * 519);
	double s;
	int k;
	int count4 = 0;
	for (int i = 0; i < 526-7; i++)
	{
		s = 0.0;
		for (int j = 1; j-1 < 8; j++)
		{
			k = (i + j) - 8;
			s += signal_extention[k + 7] * lowd[i - k];
		}
		signal_lowdconv[i] = s;
		//printf("%lf\n", signal_lowdconv[i]);
		count4++;
	}
	//printf("%d\n", count4);
	double *appro_coeff = (double*)malloc(sizeof(double) * 259);
	int count5 = 0;
	for (int i = 0; i < 259; i++)
	{
		appro_coeff[i] = signal_lowdconv[1 + (i << 1)];
		//printf("%lf\n", appro_coeff[i]);
		count5++;
	}
	//printf("%d\n", count5);

	/* Compute coefficients of detail. 1-D Convolution. downsampling */
	double *signal_highdconv = (double*)malloc(sizeof(double) * 519);
	int count6 = 0;
	for (int i = 0; i < 526 - 7; i++)
	{
		s = 0.0;
		for (int j = 1; j - 1 < 8; j++)
		{
			k = (i + j) - 8;
			s += signal_extention[k + 7] * highd[i - k];
		}
		signal_highdconv[i] = s;
		//printf("%lf\n", signal_highdconv[i]);
		count6++;
	}
	//printf("%d\n", count6);
	double *detail_coeff = (double*)malloc(sizeof(double) * 259);
	int count7 = 0;
	for (int i = 0; i < 259; i++)
	{
		detail_coeff[i] = signal_highdconv[1 + (i << 1)];
		//printf("%lf\n", detail_coeff[i]);
		count7++;
	}
	//printf("%d\n", count7);

	/* Reconstruction from coefficients of approximation. 1-D Convolution. upsampling */
	static const double lowr[8] = { 0.230377813308855, 0.714846570552542, 0.63088076792959, -0.0279837694169839, -0.187034811718881, 0.030841381835987, 0.0328830116669829, -0.0105974017849973 };
	double *appro_dyadup = (double*)malloc(sizeof(double) * 517);
	int count8 = 0;
	for (int i = 0; i < 517; i++)
	{
		appro_dyadup[i] = 0;
		//printf("%lf\n", appro_dyadup[i]);
	}
	for (int i = 0; i < 259; i++) 
	{
		appro_dyadup[i << 1] = appro_coeff[i];
	}
	for (int i = 0; i < 517; i++)
	{
		//printf("%lf\n", appro_dyadup[i]);
		count8++;
	}
	//printf("%d\n", count8);
	double *appro_lowrconv = (double*)malloc(sizeof(double) * 524);
	int i1, i2, i3;
	for (int i = 0; i < 524; i++) 
	{
		if (8 < 2 + i) 
			i1 = i - 6;
		else 
			i1 = 1;
		s = 0.0;
		      
		if (517 < 1 + i) 
			i2 = 516;
		else 
			i2 = i;
		i3 = (i2 - i1) + 1;      
		for (int j = -1; j + 1 <= i3; j++) 
		{
			k = i1 + j;
			s += appro_dyadup[k] * lowr[i - k];
		}
		appro_lowrconv[i] = s;
		//printf("%lf\n", appro_lowrconv[i]);
	}

	int dl = (524 - lx) / 2;
	int first2 = (int) 1+floor(dl);
	int last2 = (int) 524-ceil(dl);
	double *appro_lowrsignal = (double*)malloc(sizeof(double) * 512);
	int count9 = 0;
	for (int i = 0; i < last2-first2+1; i++)
	{
		appro_lowrsignal[i] = appro_lowrconv[first2-1 + i];
		//printf("%lf\n", appro_lowrsignal[i]);
		count9++;
	}
	//printf("%d\n", count9);
	
	/* Reconstruction from coefficients of detail. 1-D Convolution. upsampling */
	static const double highr[8] = { -0.0105974017849973,-0.0328830116669829,0.0308413818359870,0.187034811718881,-0.0279837694169839,-0.630880767929590,0.714846570552542,-0.230377813308855 };
	double *detail_dyadup = (double*)malloc(sizeof(double) * 517);
	int count10 = 0;
	for (int i = 0; i < 517; i++)
	{
		detail_dyadup[i] = 0;
		//printf("%lf\n", detail_dyadup[i]);
	}
	for (int i = 0; i < 259; i++)
	{
		detail_dyadup[i << 1] = detail_coeff[i];
	}
	for (int i = 0; i < 517; i++)
	{
		//printf("%lf\n", detail_dyadup[i]);
		count10++;
	}
	//printf("%d\n", count10);
	double *detail_highrconv = (double*)malloc(sizeof(double) * 524);
	for (int i = 0; i < 524; i++)
	{
		if (8 < 2 + i)
			i1 = i - 6;
		else
			i1 = 1;
		s = 0.0;

		if (517 < 1 + i)
			i2 = 516;
		else
			i2 = i;
		i3 = (i2 - i1) + 1;
		for (int j = -1; j + 1 <= i3; j++)
		{
			k = i1 + j;
			s += detail_dyadup[k] * highr[i - k];
		}
		detail_highrconv[i] = s;
		//printf("%lf\n", appro_lowrconv[i]);
	}

	double *detail_highrsignal = (double*)malloc(sizeof(double) * 512);
	int count11 = 0;
	for (int i = 0; i < last2 - first2 + 1; i++)
	{
		detail_highrsignal[i] = detail_highrconv[first2 - 1 + i];
		//printf("%lf\n", detail_highrsignal[i]);
		count11++;
	}
	printf("%d\n", count11);

	ofstream DataOut;
	DataOut.open("C:\\Users\\hhyunjjun\\Desktop\\appro_lowrsignal.txt");
	for (int j = 0; j<512; j++)
	{
		DataOut << appro_lowrsignal[j];
		DataOut << "\n";
	}
	DataOut.close();

	DataOut.open("C:\\Users\\hhyunjjun\\Desktop\\detail_highrsignal.txt");
	for (int j = 0; j<512; j++)
	{
		DataOut << detail_highrsignal[j];
		DataOut << "\n";
	}
	DataOut.close();
	return 0;
}
