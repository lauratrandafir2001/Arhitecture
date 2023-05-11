/*
 * Tema 2 ASC
 * 2023 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double *my_solver(int N, double *A, double *B)
{
	double *AxB = calloc(N * N, sizeof(double));
	double *At = calloc(N * N, sizeof(double));
	double *Bt = calloc(N * N, sizeof(double));
	double *AxBxAt_product = calloc(N * N, sizeof(double));
	double *BtxBt_product = calloc(N * N, sizeof(double));
	double *C_final_product = calloc(N * N, sizeof(double));

	if (AxBxAt_product == NULL || AxB == NULL || BtxBt_product == NULL || At == NULL || Bt == NULL || C_final_product == NULL)
		return NULL;

	// Transpose A and B
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{

			Bt[j * N + i] = B[i * N + j];
		}
	}
	// Compute AxB where A is upper triangular

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				if (i <= k)
					AxB[i * N + j] += A[i * N + k] * B[k * N + j];
			}
		}
	}

	// compute AxBxAt_product where A is upper triangular
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				if (j <= k)
					AxBxAt_product[i * N + j] += AxB[i * N + k] * A[j * N + k];
			}
		}
	}

	// compute BtxBt_product

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			for (int k = 0; k < N; k++)
			{
				BtxBt_product[i * N + j] += Bt[i * N + k] * Bt[k * N + j];
			}
		}
	}

	// compute C_final_product where C_final_product = AxBxAt_product + BtxBt_product
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			C_final_product[i * N + j] = AxBxAt_product[i * N + j] + BtxBt_product[i * N + j];
		}
	}

	free(AxB);
	free(At);
	free(Bt);
	free(AxBxAt_product);
	free(BtxBt_product);
	return C_final_product;
}
