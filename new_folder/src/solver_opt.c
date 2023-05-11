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
	double *Bt = calloc(N * N, sizeof(double));
	double *AxBxAt_product = calloc(N * N, sizeof(double));
	double *BtxBt_product = calloc(N * N, sizeof(double));
	double *C_final_product = calloc(N * N, sizeof(double));

	if (AxBxAt_product == NULL || AxB == NULL || BtxBt_product == NULL || Bt == NULL || C_final_product == NULL)
		return NULL;

	// Transpose A and B
	for (register int i = 0; i < N; i++)
	{

		for (int j = 0; j < N; j++)
		{

			Bt[j * N + i] = B[i * N + j];
		}
	}
	// Compute AxB where A is upper triangular
	register double *pAxB = AxB;
	for (register int i = 0; i < N; i++)
	{
		register double *pA_initial = &A[i * N];
		for (register int j = 0; j < N; j++)
		{
			register double *pA = pA_initial + i;
			register double *pB = &B[i * N + j];
			register double sum = 0.0;

			for (register int k = i; k < N; k++)
			{

				// AxB[i * N + j] += A[i * N + k] * B[k * N + j];
				sum += *pA * *pB;
				pA++;
				pB += N;
			}
			*pAxB = sum;
			pAxB++;
		}
	}

	// compute AxBxAt_product where A is upper triangular
	register double *pAxBxAt_product = AxBxAt_product;
	for (register int i = 0; i < N; i++)
	{
		register double *orig_pAxB = &AxB[i * N];

		for (register int j = 0; j < N; j++)
		{
			register double *pAxB = orig_pAxB + j;
			register double *pAt = &A[j * N + j];
			register double sum = 0.0;
			for (register int k = j; k < N; k++)
			{

				sum += *pAxB * *pAt;
				pAxB++;
				pAt++;
			}
			*pAxBxAt_product = sum;
			pAxBxAt_product++;
		}
	}

	// compute BtxBt_product
	register double *pBtxBt_product = BtxBt_product;
	for (register int i = 0; i < N; i++)
	{
		register double *orig_pB_T = &Bt[i * N];
		for (register int j = 0; j < N; j++)
		{
			register double *pB_T = orig_pB_T;
			register double *pB_T_second = &Bt[j];
			register double sum = 0.0;
			for (register int k = 0; k < N; k++)
			{
				// BtxBt_product[i * N + j] += Bt[i * N + k] * Bt[k * N + j];
				sum += *pB_T * *pB_T_second;
				pB_T++;
				pB_T_second += N;
			}
			*pBtxBt_product = sum;
			pBtxBt_product++;
		}
	}

	// compute C_final_product where C_final_product = AxBxAt_product + BtxBt_product
	for (register int i = 0; i < N; i++)
	{
		for (register int j = 0; j < N; j++)
		{
			C_final_product[i * N + j] = AxBxAt_product[i * N + j] + BtxBt_product[i * N + j];
		}
	}

	free(AxB);
	free(Bt);
	free(AxBxAt_product);
	free(BtxBt_product);
	return C_final_product;
}
