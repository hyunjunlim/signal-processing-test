Enter file contents here
/*RR Highfrequency extraction */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#define MAX_COLS 10000
using namespace std;

struct emxArray_real_T
{
	double *data;
	int *size;
	int allocatedsize;
	int numdimensions;
};
typedef struct emxArray_real_T emxArray_real_T;

struct emxArray__common
{
	void *data;
	int *size;
	int allocatedsize;
	int numdimensions;
};
typedef struct emxArray__common emxArray__common;
void emxInit_real_T(emxArray_real_T **pEmxArray, int numdimensions);
void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int elementSize);
void mydwt(emxArray_real_T *x, double lowd[], double highd[], emxArray_real_T *appro_coeff, emxArray_real_T *detail_coeff);
void myhighwrcoef(emxArray_real_T *C, emxArray_real_T *CL, int level, emxArray_real_T *upsy);
int main(int argc, char* argv[])
{

	/* RRseries load */
	FILE *in;
	char s[MAX_COLS]; // 행이 1줄씩 임시로 저장될 버퍼
	char filea;

	in = fopen("C:/Users/hhyunjjun/Desktop/matlabC/SPR006_PACU_RRseries.csv", "r"); // data 불러오기

	int line = 0;
	while ((filea = fgetc(in)) != EOF)
	{  // 읽은 data 전체 개수 세기
		if (filea == '\n') line++;
	}
	//printf("%d\n", line);

	int count1 = 0;

	double *rrdata = (double*)malloc(sizeof(double) * line); // 메모리 동적 할당
	fseek(in, 0, SEEK_SET);
	while (fgets(s, MAX_COLS, in) != NULL)
	{
		rrdata[count1] = atof(s);
		count1++;
	}
	//printf("%lf\n", data[479]);

	/*  Initialization. */
	double lowd[8] = { -0.0105974017849973, 0.0328830116669829, 0.030841381835987, -0.187034811718881, -0.0279837694169839, 0.63088076792959, 0.714846570552542, 0.230377813308855 };
	double highd[8] = { -0.230377813308855, 0.714846570552542, -0.63088076792959, -0.0279837694169839, 0.187034811718881, 0.030841381835987, -0.0328830116669829, -0.0105974017849973 };
	int n = 8;
	int k = 0;
	int i, j;
	int loop_ub;

	emxArray_real_T *y;
	emxInit_real_T(&y, 2);
	i = y->size[0] * y->size[1];
	y->size[0] = 1;
	y->size[1] = 480;
	emxEnsureCapacity((emxArray__common *)y, i, (int)sizeof(double));
	loop_ub = 1 * 480;
	for (i = 0; i < loop_ub; i++)
	{
		y->data[i] = rrdata[i];
	}

	emxArray_real_T *c;
	emxInit_real_T(&c, 2);
	i = c->size[0] * c->size[1];
	c->size[0] = 1;
	c->size[1] = 0;
	emxEnsureCapacity((emxArray__common *)c, i, (int)sizeof(double));

	emxArray_real_T *cl;
	emxInit_real_T(&cl, 2);
	i = cl->size[0] * cl->size[1];
	cl->size[0] = 1;
	cl->size[1] = 10;
	emxEnsureCapacity((emxArray__common *)cl, i, (int)sizeof(double));
	loop_ub = 8 + 2;
	for (i = 0; i < loop_ub; i++)
	{
		cl->data[i] = 0.0;
	}
	cl->data[9] = y->size[1];

	/*  decomposition / store detail, length */
	emxArray_real_T *a;
	emxArray_real_T *d;
	emxArray_real_T *b_y;
	emxArray_real_T *b_d;


	emxInit_real_T(&a, 2);
	emxInit_real_T(&d, 2);
	emxInit_real_T(&b_y, 2);
	emxInit_real_T(&b_d, 2);

	i = b_y->size[0] * b_y->size[1];
	b_y->size[0] = 1;
	b_y->size[1] = y->size[1];
	emxEnsureCapacity((emxArray__common *)b_y, i, (int)sizeof(double));
	loop_ub = y->size[0] * y->size[1];
	for (i = 0; i < loop_ub; i++)
	{
		b_y->data[i] = y->data[i];
	}
	while (k <= 7)
	{
		mydwt(b_y, lowd, highd, a, d); // dwt

		i = b_y->size[0] * b_y->size[1];
		b_y->size[0] = 1;
		b_y->size[1] = a->size[1];
		emxEnsureCapacity((emxArray__common *)b_y, i, (int)sizeof(double));
		loop_ub = a->size[0] * a->size[1];
		for (i = 0; i < loop_ub; i++)
		{
			b_y->data[i] = a->data[i];
		}

		i = b_d->size[0] * b_d->size[1];
		b_d->size[0] = 1;
		b_d->size[1] = d->size[1] + c->size[1];
		emxEnsureCapacity((emxArray__common *)b_d, i, (int)sizeof(double));
		loop_ub = d->size[1];
		for (i = 0; i < loop_ub; i++)
		{
			b_d->data[i] = d->data[i];
		}
		loop_ub = c->size[1];
		for (i = 0; i < loop_ub; i++)
		{
			b_d->data[i + d->size[1]] = c->data[i];
		}

		i = c->size[0] * c->size[1];
		c->size[0] = 1;
		c->size[1] = b_d->size[1];
		emxEnsureCapacity((emxArray__common *)c, i, (int)sizeof(double));
		loop_ub = b_d->size[1];
		for (i = 0; i < loop_ub; i++)
		{
			c->data[i] = b_d->data[i];
		}
		cl->data[10 - (1 + k) - 1] = d->size[1];
		k++;
	}

	emxArray_real_T *c_y;
	emxInit_real_T(&c_y, 2);

	/*  Last approximation. */
	i = c_y->size[0] * c_y->size[1];
	c_y->size[0] = 1;
	c_y->size[1] = b_y->size[1] + c->size[1];
	emxEnsureCapacity((emxArray__common *)c_y, i, (int)sizeof(double));
	loop_ub = c_y->size[1];
	for (i = 0; i < loop_ub; i++)
	{
		c_y->data[i] = b_y->data[i];
	}
	loop_ub = c->size[1];
	for (i = 0; i < loop_ub; i++)
	{
		c_y->data[i + b_y->size[1]] = c->data[i];
	}

	i = c->size[0] * c->size[1];
	c->size[0] = 1;
	c->size[1] = c_y->size[1];
	emxEnsureCapacity((emxArray__common *)c, i, (int)sizeof(double));
	loop_ub = c_y->size[1];
	for (i = 0; i < loop_ub; i++)
	{
		c->data[i] = c_y->data[i];
	}
	cl->data[0] = b_y->size[1];

	/* Wrcoef */
	emxArray_real_T *details;
	emxInit_real_T(&details, 2);

	i = details->size[0] * details->size[1];
	details->size[0] = 1;
	details->size[1] = 480;
	emxEnsureCapacity((emxArray__common *)details, i, (int)sizeof(double));
	loop_ub = details->size[1];

	int height = 8, width = 480;
	double **details1;
	details1 = (double**)malloc(sizeof(double*) * height);
	for (i = 0; i < height; i++)
	{
		details1[i] = (double*)malloc(sizeof(double) * width);
	}
	
	
	for (i = 0; i < 8; i++)
	{
		myhighwrcoef(c, cl, i+1, details);
		for (j = 0; j < 480; j++)
		{
			details1[i][j] = details->data[j];
		}
	}
	/*
	myhighwrcoef(c, cl, 8, details);
	for (j = 0; j < 480; j++)
	{
		printf("%lf\n", details->data[j]);
	}
	*/
	double *RRHF = (double*)malloc(sizeof(double) * width); // 메모리 동적 할당
	for (j = 0; j < width; j++)
	{
		RRHF[j] = details1[3][j]+ details1[4][j];
		printf("%lf\n", RRHF[j]);
	}

	ofstream DataOut;
	DataOut.open("C:\\Users\\hhyunjjun\\Desktop\\RRHF.txt");
	for (j = 0; j < 480; j++)
	{
		DataOut << RRHF[j];
		DataOut << "\n";
	}
	DataOut.close();
	
	return 0;
}

void emxInit_real_T(emxArray_real_T **pEmxArray, int numdimensions)
{
	emxArray_real_T *emxArray;
	int i;
	*pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
	emxArray = *pEmxArray;
	emxArray->data = (double *)NULL;
	emxArray->numdimensions = numdimensions;
	emxArray->size = (int *)malloc((unsigned int)(sizeof(int) * numdimensions));
	emxArray->allocatedsize = 0;
	for (i = 0; i < numdimensions; i++) 
	{
		emxArray->size[i] = 0;
	}
}

void emxEnsureCapacity(emxArray__common *emxArray, int oldNumel, int elementSize)
{
     int newNumel;
     int i;
     void *newData;
     newNumel = 1;
     for (i = 0; i < emxArray->numdimensions; i++) 
	 {
       newNumel *= emxArray->size[i];
	 }

     if (newNumel > emxArray->allocatedsize)
	 {
       i = emxArray->allocatedsize;
       if (i < 16) 
	   {
         i = 16;
	   }

       while (i < newNumel) 
	   {
         i <<= 1;
	   }

       newData = calloc((unsigned int)i, (unsigned int)elementSize);
       if (emxArray->data != NULL) 
	   {
         memcpy(newData, emxArray->data, (unsigned int)(elementSize * oldNumel));
}
       emxArray->data = newData;
       emxArray->allocatedsize = i;
}
   }

void mydwt(emxArray_real_T *x, double low[], double high[], emxArray_real_T *appro_coeff, emxArray_real_T *detail_coeff)
{
	int shift = 0;

	// Compute sizes and shape.
	int lf = 8;
	int lx = x->size[1];

	/*Signal Extensions: Dealing with Border Distortion */
	// Extend, Decompose &  Extract coefficients.
	int first1 = 2 - shift;
	int lenEXT = lf - 1;
	int last1 = lx + lf - 1;
	int *lenEXTarray = (int*)malloc(sizeof(int) * 7);
	int *lxarray = (int*)malloc(sizeof(int) * lx);
	int *lxlenEXTarray = (int*)malloc(sizeof(int) * 7);

	for (int i = 0; i < lenEXT; i++)
	{
		lenEXTarray[i] = lenEXT - i - 1;
	}
	for (int i = 0; i < lx; i++)
	{
		lxarray[i] = i;
	}
	for (int i = 0; i < lenEXT; i++)
	{
		lxlenEXTarray[i] = lx - i - 1;
	}
	int *l = (int*)malloc(sizeof(int) * (last1 + lenEXT));
	for (int i = 0; i < (last1 + lenEXT); i++)
	{
		if (i < lenEXT)
			l[i] = lenEXTarray[i];
		else if (lenEXT <= i && i < lenEXT + lx)
			l[i] = lxarray[i - lenEXT];
		else if (lenEXT + lx <= i)
			l[i] = lxlenEXTarray[i - lenEXT - lx];
	}

	emxArray_real_T *signal_extention;
	emxInit_real_T(&signal_extention, 2);
	int emxsize1 = signal_extention->size[0] * signal_extention->size[1];
	signal_extention->size[0] = 1;
	signal_extention->size[1] = (last1 + lenEXT);
	emxEnsureCapacity((emxArray__common *)signal_extention, emxsize1, (int)sizeof(double));
	for (int i = 0; i < (last1 + lenEXT); i++)
	{
		signal_extention->data[i] = x->data[l[i]];
	}

	/* Compute coefficients of approximation. 1-D Convolution. downsampling */
	emxArray_real_T *signal_lowdconv;
	emxInit_real_T(&signal_lowdconv, 2);
	int emxsize2 = signal_lowdconv->size[0] * signal_lowdconv->size[1];
	signal_lowdconv->size[0] = 1;
	signal_lowdconv->size[1] = last1;
	emxEnsureCapacity((emxArray__common *)signal_lowdconv, emxsize2, (int)sizeof(double));

	double s;
	int k;
	for (int i = 0; i < last1; i++)
	{
		s = 0.0;
		for (int j = 1; j - 1 < 8; j++)
		{
			k = (i + j) - 8;
			s += signal_extention->data[k + 7] * low[i - k];
		}
		signal_lowdconv->data[i] = s;
	}

	int emxsize3 = appro_coeff->size[0] * appro_coeff->size[1];
	appro_coeff->size[0] = 1;
	appro_coeff->size[1] = (((int)last1 - 2) >> 1) + 1;
	emxEnsureCapacity((emxArray__common *)appro_coeff, emxsize3, (int)sizeof(double));

	for (int i = 0; i < (((int)last1 - 2) >> 1) + 1; i++)
	{
		appro_coeff->data[i] = signal_lowdconv->data[1 + (i << 1)];
	}

	/* Compute coefficients of detail. 1-D Convolution. downsampling */
	emxArray_real_T *signal_highdconv;
	emxInit_real_T(&signal_highdconv, 2);
	int emxsize4 = signal_highdconv->size[0] * signal_highdconv->size[1];
	signal_highdconv->size[0] = 1;
	signal_highdconv->size[1] = last1;
	emxEnsureCapacity((emxArray__common *)signal_highdconv, emxsize4, (int)sizeof(double));

	for (int i = 0; i < last1; i++)
	{
		s = 0.0;
		for (int j = 1; j - 1 < 8; j++)
		{
			k = (i + j) - 8;
			s += signal_extention->data[k + 7] * high[i - k];
		}
		signal_highdconv->data[i] = s;
	}

	int emxsize5 = detail_coeff->size[0] * detail_coeff->size[1];
	detail_coeff->size[0] = 1;
	detail_coeff->size[1] = (((int)last1 - 2) >> 1) + 1;
	emxEnsureCapacity((emxArray__common *)detail_coeff, emxsize5, (int)sizeof(double));

	for (int i = 0; i < (((int)last1 - 2) >> 1) + 1; i++)
	{
		detail_coeff->data[i] = signal_highdconv->data[1 + (i << 1)];
	}
}

void myhighwrcoef(emxArray_real_T *C, emxArray_real_T *CL, int level, emxArray_real_T *upsy)
{
	double x_first[10];
	double xfirst;
	int k;
	int ix = 0;
	double b_x[10];
	double lastcl[8];
	int i;
	double first[8];
	double last[8];
	int j;
	int j1, j2, j3;
	int loop_ub;
	
	double s;
	static const double highr[8] = { -0.0105974017849973, -0.0328830116669829,
		0.030841381835987, 0.187034811718881, -0.0279837694169839, -0.63088076792959,
		0.714846570552542, -0.230377813308855 };

	static const double lowr[8] = { 0.230377813308855, 0.714846570552542,
		0.63088076792959, -0.0279837694169839, -0.187034811718881, 0.030841381835987,
		0.0328830116669829, -0.0105974017849973 };

	x_first[0] = CL->data[0];
	xfirst = CL->data[0];
	for (k = 0; k < 9; k++) 
	{
		ix++;
		xfirst += CL->data[ix];
		x_first[ix] = xfirst;
	}

	for (i = 0; i < 10; i++)
	{
		b_x[i] = x_first[i] + 1.0;
	}

	for (i = 0; i < 8; i++) 
	{
		first[i] = b_x[7 - i];
	}

	for (i = 0; i < 8; i++)
	{
		lastcl[i] = CL->data[8 - i];
	}

	for (i = 0; i < 8; i++) 
	 {
		last[i] = (first[i] + lastcl[i]) - 1.0;
	}

	int i1 = (int)first[level - 1] - 1;
	int i2 = (int)last[level - 1] ;

	/*  Compute Upsampling and Convolution. */
	emxArray_real_T *dyadup;
	emxInit_real_T(&dyadup, 2);

	int ll = 2 * (i2 - i1) - 1;
	j = dyadup->size[0] * dyadup->size[1];
	dyadup->size[0] = 1;
	dyadup->size[1] = ll;
	emxEnsureCapacity((emxArray__common *)dyadup, j, (int)sizeof(double));
	loop_ub = ll;
	for (j = 0; j < loop_ub; j++) {
		dyadup->data[j] = 0.0;
	}
	loop_ub = i2 - i1;
	for (i = 0; i < loop_ub; i++) 
	{
		dyadup->data[2 * i] = C->data[i1 + i];
	}

	emxArray_real_T *arg1;
	emxInit_real_T(&arg1, 2);

	i = arg1->size[0] * arg1->size[1];
	arg1->size[0] = 1;
	arg1->size[1] = dyadup->size[1];
	emxEnsureCapacity((emxArray__common *)arg1, i, (int)sizeof(double));
	
	int xls;
	for (i = 0; i < arg1->size[1]; i++)
	{
		arg1->data[i] = dyadup->data[i];
	}
	j = dyadup->size[0] * dyadup->size[1];
	dyadup->size[0] = 1;
	dyadup->size[1] = arg1->size[1]+ 7U;
	emxEnsureCapacity((emxArray__common *)dyadup, j, (int)sizeof(double));
	
	for (i = 0; i < dyadup->size[1]; i++)
	{
		if (8 < 2 + i)
			j1 = i - 6;
		else
			j1 = 1;

		xls = arg1->size[1];
		if (xls < 1+i)
		{
			j2 = xls-1;
		}
		else
			j2 = i;

		s = 0.0;
		j3 = j2 - j1 + 1;
		for (j = -1; j+1 <= j3; j++) 
		{
			k = j1 + j;
			s += arg1->data[k] * highr[i-k];
		}
		dyadup->data[i] = s;
		
	}
	
	double dl = (dyadup->size[1] - CL->data[10 - level])/ 2;
	
	int first2 = (int)1 + floor(dl);
	int last2 = (int)dyadup->size[1] - ceil(dl);
	
	i = upsy->size[0] * upsy->size[1];
	upsy->size[0] = 1;
	upsy->size[1] = last2 - first2 + 1;
	emxEnsureCapacity((emxArray__common *)upsy, i, (int)sizeof(double));
	loop_ub = upsy->size[1];
	for (i = 0; i < loop_ub; i++) 
	{
		upsy->data[i] = dyadup->data[first2 - 1 + i];
	}
	
	for (k = 1; k < level; k++)
	{
		ll = 2 * upsy->size[1] - 1;
		
		i = dyadup->size[0] * dyadup->size[1];
		dyadup->size[0] = 1;
		dyadup->size[1] = ll;
		emxEnsureCapacity((emxArray__common *)dyadup, i, (int)sizeof(double));
		loop_ub = ll;
		for (i = 0; i < loop_ub; i++) 
		{
			dyadup->data[i] = 0.0;
		}
		loop_ub = upsy->size[1];
		for (i = 0; i < loop_ub; i++) 
		{
			dyadup->data[2 * i] = upsy->data[i];
		}

		emxArray_real_T *arg1;
		emxInit_real_T(&arg1, 2);

		i = arg1->size[0] * arg1->size[1];
		arg1->size[0] = 1;
		arg1->size[1] = dyadup->size[1];
		emxEnsureCapacity((emxArray__common *)arg1, i, (int)sizeof(double));
		int xls;
		for (i = 0; i <  arg1->size[1]; i++)
		{
			arg1->data[i] = dyadup->data[i];
		}

		j = dyadup->size[0] * dyadup->size[1];
		dyadup->size[0] = 1;
		dyadup->size[1] = arg1->size[1] + 7U;
		emxEnsureCapacity((emxArray__common *)dyadup, j, (int)sizeof(double));
		for (i = 0; i < dyadup->size[1]; i++)
		{
			if (8 < 2 + i)
				j1 = i - 6;
			else
				j1 = 1;

			xls = arg1->size[1];
			if (xls < 1 + i)
			{
				j2 = xls-1;
			}
			else
				j2 = i;

			s = 0.0;
			j3 = j2 - j1 + 1;
			for (j = -1; j+1 <= j3; j++)
			{
				int z = j1 + j;
				s += arg1->data[z] * lowr[i - z];
			}
			dyadup->data[i] = s;
		}
		
		dl = (dyadup->size[1] - CL->data[10 - level + k]) / 2;
		
		first2 = (int)1 + floor(dl);
		last2 = (int)dyadup->size[1] - ceil(dl);
		i = upsy->size[0] * upsy->size[1];
		upsy->size[0] = 1;
		upsy->size[1] = last2 - first2 + 1;
		emxEnsureCapacity((emxArray__common *)upsy, i, (int)sizeof(double));
		loop_ub = upsy->size[1];
		for (i = 0; i < loop_ub; i++)
		{
			upsy->data[i] = dyadup->data[first2 - 1 + i];
		}
	}
}
