# Code Optimisation for Matrix Multiplication

# About
Computing Systems Architecture course
https://ocw.cs.pub.ro/courses/asc/teme/tema2

May 2023

Student:  Trandafir Laura;

# General Description:

We needed to implement C=A×B×At+Bt×Bt where A is upper triangular and At and Bt 
are A and B transpose
We had to implement this operation in 3 different ways:
blas :which uses BLAS Atlas functions in order to make the operations
neopt- in the simple way, the classical one
opt_m - it's the optimaized implementation for neaop, using pointers
and other tricks in order for the time of running to be shorter

# Implementation:
    -solver_blas: First, I used the cblas_dtrmm function, to calculate
            the product between A and B and stores it in AxB_product
            After that I also used cblas_dtrmm, to do AxB_product x At
            but this time I specified using CblasTrans that A is transpose.
            I chose to use this function because it considers that A is upper triangular
            I store the result in AxBxAt_product
            For the Bt x Bt product I used cblas_dgemm, also specifing that B is transpose
            using the CblasTrans argument, and in the end I store it in BtxBt_product
            I use cblas_daxpy to compute the final operations, the sum between
            AxBxAt_product and BtxBt_product
            I free the memory, the above matrixes used, and then return the final result
    -solver_neopt: For the simple, classical approach, we multiply the matrixes element by element
                    First I use to loops to transpose the B matrix, after that I compute AxB
                    considering A is upper triangular we will only consider the elements 
                    which have the row number smaller than k
                    After this, we will compute AxBxAt where we will notice that At is lower triangular
                    and will consider only the elements which have the colomn number smaller than k
                    In the end we will compute BtxBt using the transpose we calculated at the begining
                    and make the final sum, storing it into C_final_product

    -solver_opt: This approach, by far the fastest one(considering the running time) but also the most difficult
                I also transpose the B matrix
                I computes AxB product, but it memories a pointer to  &A[i * N] right after the loop for rows,
                to help the computer not to calculate it every time the colomn number changes
                inside the j loop we add i to the reference, the matrix being superior triangular
                and also take another reference for B matrix &B[i * N + j]
                Inside the 3rd(k loop) we calculate the product and increase the pointers
                I also save a pointer to the matrix where I want to save the result and increase it step
                by step
                I do the same thing for AxBxAt, and also used "register as sugessted in the laboratories
                to make it faster
                In the end I compute the final sum and store it into C_final_product
# CacheGrind
We have 3 files neopt.cache, blas.cache and opt_m.cache that contains the output of the following valgrind command
"–tool=cachegrind –branch-sim=yes"
The output contains:

I refs: This refers to the number of instruction references, i.e., the number of times an instruction was executed.
D refs: This refers to the number of data references, i.e., the number of times data was read from or written to memory.
Branches: This refers to the number of branch instructions executed, i.e., the number of times the program had to 
decide which instruction to execute next based on a condition.
Mispredicts: This refers to the number of times the branch predictor made an incorrect prediction about which instruction to execute next.
Mispred rate: This is the percentage of branch mispredictions out of the total number of branches executed.
# Graphics
I also attached 4 graphics, one for each implementation and one with all the implementation.

# Conclusions:
The optimized implementations is the best one, but also the blas one is really close.
I consider that the blas was also easier to make, given those functions, and the optimized implementation
is harder to debug, you can use it for multiplying matrixes but if you want to make something more twisted,
the implementation using pointers will become pretty hard to read and understand by other(if more people work
on the same project). On the other hand, the "register" key-word made really significant changes on the times, and
it's easy to use.
References:
- ["Matrix Multiplication improved"](https://ocw.cs.pub.ro/courses/asc/laboratoare/05)
- [Blas Function](https://netlib.org/blas/)
- [Matrix multiplication using pointers](https://www.tutorialspoint.com/how-to-multiply-two-matrices-using-pointers-in-c)


