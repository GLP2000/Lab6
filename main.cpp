#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include<stdlib.h>
#include <math.h>
#include "winbgi2.h"
#include "gauss.h"

void HilbertMatrix(int N, double** H);
void displayMatrix(int N, double** H);
void computeVec(int N, double** H, double* B);
void plotVec(int N, double* V);

void computeMatrix2(int N, double** K);
void displayMatrix2(int N, double** K);
void computeVector2(int N, double* F);
void displayVector2(int N, double* F);
void grafik(int N, double* T);
double determinant(int N, double** K);

double lambda = 58.0;




int main()
{
	int w, N;
	printf("Wybierz 1. Macierz Huberta lub 2.Rozklad temperatury w precie\n");
	scanf("%d", &w);
	

	if (w == 1)
	{
		printf("podaj wartosc N-liczba rownan\n");
		scanf("%d", &N);

		double** H;
		H = (double**)malloc(N * sizeof(double*));
		for (int i = 0; i < N; i++)
		{
			H[i] = (double*)malloc(N * sizeof(double));
		}
		double* X, * B;
		X = (double*)malloc(N * sizeof(double));
		B = (double*)malloc(N * sizeof(double));

		for (int i = 0; i < N; i++)
		{
			B[i] = 0;
			X[i] = 0;
		}

		HilbertMatrix(N, H);
		//displayMatrix(N, H);
		computeVec(N, H, B);

		gauss(N, H, X, B);
		plotVec(N, X);
		for (int i = 0; i < N; i++)
		{
			free(H[i]);
		}
		free(H);
		free(X);
		free(B);
	}

	else if (w == 2)
	{
		printf("podaj wartosc N-liczba rownan\n");
		scanf("%d", &N);

		N = N + 1;

		double** K;
		K = (double**)malloc(N * sizeof(double*));
		for (int i = 0; i < N; i++)
		{
			K[i] = (double*)malloc(N * sizeof(double));
		}
		double* T, * F;
		T = (double*)malloc(N * sizeof(double));
		F = (double*)malloc(N * sizeof(double));

		for (int i = 0; i < N; i++)
		{
			T[i] = 0;
			F[i] = 0;
		}

		computeMatrix2(N, K);
		displayMatrix2(N, K);
		computeVector2(N, F);
		displayVector2(N, F);
		gauss(N, K, T, F);
		displayVector2(N, T);
		grafik(N, T);

		determinant(N, K);
		printf("det=%lf", determinant(N, K));




		for (int i = 0; i < N; i++)
		{
			free(K[i]);
		}
		free(K);
		free(T);
		free(F);
	}

	else
	{
		printf("Wybierz znowu 1. Macierz Huberta lub 2.Rozklad temperatury w precie\n");
		scanf("%d", &w);
		printf("podaj wartosc N-liczba rownan\n");
		scanf("%d", &N);
	}



	wait();
}

void HilbertMatrix(int N, double** H)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			H[i][j] = 1.0 / ((double)(1 + i + j));
		}
	}
}
void displayMatrix(int N, double** H)
{
	printf("drukujemy dla macierzy H\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%lf\t", H[i][j]);
		}
		printf("\n");
	}
}
void computeVec(int N, double** H, double* B)
{
	for (int i = 0; i < N; i++)
	{
		B[i] = 0;
		for (int j = 0; j < N; j++)
		{
			B[i] += H[i][j];
		}
	}
}
void plotVec(int N, double* V)
{
	printf("Macierz rowna sie:\n");
	for (int i = 0; i < N; i++)
	{
		printf("%lf\n", V[i]);

	}
	printf("\n");
}
void computeMatrix2(int N, double** K)
{
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			K[i][j] = 0.0;
		}
	}
	for (int u = 0; u < (N - 2); u++)
	{
		K[u + 1][u] = 1.0;
		K[u + 1][u + 1] = -2.0;
		K[u + 1][u + 2] = 1.0;
	}
	K[0][0] = 1.0;
	K[N - 1][N - 1] = 1.0;
}
void displayMatrix2(int N, double** K)
{
	printf("drukujemy dla macierzy K\n");
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			printf("%lf\t", K[i][j]);
		}
		printf("\n");
	}
}
void computeVector2(int N, double* F)
{
	double h = 1.0 / N; 
	double X = 0;
	for (int i = 0; i < N; i++)
	{
		X += h;
		F[i] = (-(double)(pow(10.0, 4.0) * (double)sin(X * 3.1415)) / (double)lambda) * (double)pow(h, 2.0);
	}
	F[0] = 273;// tp=273 K
	F[N - 1] = 300;//tk=300 K
}
void displayVector2(int N, double* F)
{
	printf("Macierz rowna sie:\n");
	for (int i = 0; i < N; i++)
	{
		printf("%lf\n", F[i]);

	}
	printf("\n");
}
void grafik(int N, double* T)
{
	graphics(800, 600);
	scale(0, 275, 1, 315);
	double p = 0;
	for (int i = 0; i < N; i++)
	{
		setcolor(0.5);
		point(p, T[i]);
		lineto(p, T[i]);
		p += 1.0 / N;
	}
}
double determinant(int N, double** K)
{
	double wsp, det=1;
	for (int i = 0; i < (N - 1); i++)
	{
		for (int j =0 ; j < N; j++)
		{
			if (j > i)
			{
				wsp = K[j][i] / K[i][i];
				for (int k = 0; k < N; k++)
				{
					K[j][k] -= wsp* K[i][k];
				}
			}
		}
	}
	
	for (int i = 1; i < N; i++)
	{
			det *= K[i][i];
	}
	return det;
}


