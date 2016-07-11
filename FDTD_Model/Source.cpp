// cpp_compiler_options_openmp.cpp
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <windows.h>


void calcMatrix(double*** Ez, int nx, int ny, int nz)
{
	srand((unsigned)time(NULL));
//	printf_s("ok\n");
#pragma omp parallel for
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			for (int k = 0; k < nz; ++k) {
				double sigma = rand() / double(RAND_MAX)+1.0;
				double epsilon = rand() / double(RAND_MAX)+1.0;
				Ez[i][j][k] = (2 * epsilon - sigma) / (2 * epsilon + sigma) * Ez[i][j][k] + (sigma/(2*epsilon+sigma))*(Ez[i][j][k]-epsilon * Ez[i][j][k]);
				//Ez[i][j][k] = Ez[i][j][k] + 0.1;
				//printf("%.2lf ", Ez[i][j][k]);
			}
		}
	}
}
int main()
{
	double*** Ez = NULL;
	int nx = 700, ny = 700, nz = 700;
	Ez = new double**[nx];
	for (int i = 0; i < nx; ++i) {
		Ez[i] = new double*[ny];
		for (int j = 0; j < ny; ++j) {
			Ez[i][j] = new double[nz];
			memset(Ez[i][j], 0, sizeof(double)*nz);
		}
	}
//	printf_s("calc done\n");
	DWORD dwStart = GetTickCount();
//	printf_s("start\n");
	for (int t = 0; t < 1; ++t) {
		calcMatrix(Ez, nx, ny, nz);
	}
//	printf_s("end\n");
	printf_s("%d ms\n", GetTickCount() - dwStart);
	system("pause");
	return 0;
}