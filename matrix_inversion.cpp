/* matrix inversion */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matrix_inversion.h"

void swap(double* a, double* b)
{
	double c;
	c = *a;
	*a = *b;
	*b = c;
}

int inv(double* p, int n)
{
	int *is, *js;
	int i, j, k;

#if 0
	// print matrix;
	for (i = 0; i < n; ++i)
	{
		putchar('\n');
		for (j = 0; j < n; ++j)
			printf("%f   ", *(p + i * n + j));
	}
	puts("\n\n\n\n");
#endif

	double temp, fmax;
	is = (int*)malloc(n * sizeof(int));
	js = (int*)malloc(n * sizeof(int));

	for (k = 0; k < n; ++k)
	{
		fmax = 0.0;
		for (i = k; i < n; ++i)
		{
			for (j = k; j < n; ++j)
			{
				temp = fabs(*(p + i*n + j));
				if (temp > fmax)
				{
					fmax = temp;
					is[k] = i;
					js[k] = j;
				}
			}
		}
		if ((fmax + 1.0) == 1.0)
		{
			free(is);
			free(js);
			return 0;
		}

		if ((i = is[k]) != k)
			for (j = 0; j < n; ++j)
				swap(p+k*n+j, p+i*n+j);
		if ((j = js[k]) != k)
			for (i = 0; i < n; ++i)
				swap(p+i*n+k, p+i*n+j);
		p[k*n+k] = 1.0/p[k*n+k];

		for (j = 0; j < n; ++j)
			if (j != k)
				p[k*n+j] *= p[k*n+k];
		for (i = 0; i < n; ++i)
			if (i != k)
				for (j = 0; j < n; ++j)
					if (j != k)
						p[i*n+j] = p[i*n+j] - p[i*n+k]*p[k*n+j];
		for (i = 0; i < n; ++i)
			if (i != k)
				p[i*n+k] *= -p[k*n+k];
	}
	for (k = n-1; k >= 0; --k)
	{
		if ((j = js[k]) != k)
			for (i = 0; i < n; ++i)
				swap(p+j*n+i, p+k*n+i);
		if ((i=is[k]) != k)
			for (j = 0; j < n; ++j)
				swap(p+j*n+i, p+j*n+k);
	}
	free(is);
	free(js);
	return 1;
}
