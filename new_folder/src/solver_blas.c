
#include <stdio.h>
#include <string.h>
#include "utils.h"
#include <cblas.h>

double *my_solver(int N, double *A, double *B)
{
	double *AxBxAt_product = calloc(N * N, sizeof(double));
	double *AxB_product = calloc(N * N, sizeof(double));
	double *BtxBt_product = calloc(N * N, sizeof(double));
	if (AxB_product == NULL || AxBxAt_product == NULL || BtxBt_product == NULL)
		return NULL;

	memcpy(AxB_product, B, N * N * sizeof(double));
	// compute AxB and store it in AxB_product
	cblas_dtrmm(
		CblasRowMajor,
		CblasLeft,
		// A is upper triangular
		CblasUpper,
		CblasNoTrans,
		CblasNonUnit,
		N,
		N,
		1,
		A, N,
		AxB_product, N);
	// compute AxB_productxAt and store it in AxBxAt_product
	memcpy(AxBxAt_product, AxB_product, N * N * sizeof(double));
	cblas_dtrmm(
		CblasRowMajor,
		CblasRight,
		CblasUpper,
		// AT
		CblasTrans,
		CblasNonUnit,
		N,
		N,
		1,
		A, N,
		AxBxAt_product, N);

	// compute BtxBt and store it in BtxBt_product
	cblas_dgemm(
		CblasRowMajor,
		CblasTrans,
		CblasTrans,
		N, N,
		N, 1.0,
		B, N,
		B, N,
		1.0, BtxBt_product, N);
	// compute the sum of AxBxAt_product and BtxBt_product and store it in BtxBt_product
	cblas_daxpy(N * N, 1.0, AxBxAt_product, 1, BtxBt_product, 1);
	free(AxB_product);
	free(AxBxAt_product);

	return BtxBt_product;
}
