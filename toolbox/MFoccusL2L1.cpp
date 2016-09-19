/* MFoccusL2L1.CPP Group sparse coding
 * 
 * This code solves the optimization problem for each group of signals
 *
 * min_P 0.5*||X - D*P||^2 + gamma*SUM_1:K ||P(k,:)||_2
 *
 * Syntax:      P = MFoccusL2L1(X, D, G, P0, param)
 *              P = MFoccusL2L1(X, D, G, P0)
 *
 * Inputs:
 *          D:          Dictionary
 *          X:          Signals
 *          P0:         Initial Guess
 *          G:          Array of zero-based indexes into start of each group
 *          param:      Struct array of following variables
 *              gamma:      Regularization parameter (default: 1)
 *              tol:        Tolerance (default: 10^-2)
 *              max_iter:   Maximum number of iterations (default: 100)
 *              nthreads:   Number of threads to use (default: No. of CPUs)
 *
 * Output:
 *          P:          Group Sparse Coefficients
 *
 *
 * Written by Ali Soltani-Farani <a_soltani@ce.sharif.edu>
 * Copyright 2012 by Ali Soltani-Farani                                 */

 
/* MATLAB Code for a single group
 * function P = MFOCCUS_L2L1_N(X, D, P, gamma, tol, max_iter )
 *
 *   gI = gamma*eye(size(D,1),size(D,1));
 *   
 *   delta = 1;
 *   iter = 1;
 *   while delta > tol && iter <= max_iter
 *       L = sum(P.^2,2).^0.25;
 *       Dhat = bsxfun(@times,D,L');
 *       newP = bsxfun( @times, Dhat'*((gI+Dhat*Dhat')\X), L );
 *       
 *       delta = sum(sum((P-newP).^2))/sum(sum(P.^2));
 *       P = newP;
 *       iter = iter + 1;
 *   end
 *   
 * end                                                                  */

#include <mex.h>
#include <blas.h>
#include <lapack.h>
#include <math.h>
#include <new>

#ifdef _OPENMP
#include <omp.h>
#endif

int my_lapack_dgesv( int N, int NRHS, double *A, int LDA, mwSignedIndex *ipiv, double *B, int LDB);
void my_blas_dgemm ( const char TransA, const char TransB, const int M, const int N, const int K, const double ALPHA, const double * A, const int LDA, const double * B, const int LDB, const double BETA, double * C, const int LDC);
void my_blas_daxpy( int N, double alpha, double *X, int INCX, double *Y, int INCY);
double my_blas_ddot( int N, double *X, int INCX, double *Y, int INCY);
double my_blas_dnrm2(int N,double *X,int INCX);
void my_blas_dcopy(int N, double *X, int INCX, double* Y, int INCY );
void my_blas_dscal( int N, double ALPHA, double * X, int INCX );
void core(double *X, int M, double *D, int K, double *P0, double *P, double *G, int N_G, double gamma, double tol, double max_iter, int nthreads );
void mfoccus_core(double *X, int M, double *D, int K, double *P0, double *P, double *G, int N_G, double gamma, double tol, double max_iter, double *L, double *DL, double *DLDt, double *DLDt_inv_X, mwSignedIndex *ipiv, int tid, int nthreads);


using namespace std;

#define IS_REAL_2D_FULL_DOUBLE(P) (!mxIsComplex(P) && \
mxGetNumberOfDimensions(P) == 2 && !mxIsSparse(P) && mxIsDouble(P))
#define IS_REAL_SCALAR(P) (IS_REAL_2D_FULL_DOUBLE(P) && mxGetNumberOfElements(P) == 1)

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    /* Macros for the input and output arguments */
    #define P_OUT         plhs[0]       
    #define X_IN          prhs[0]
    #define D_IN          prhs[1]
    #define G_IN          prhs[2]
    #define P0_IN         prhs[3]
    #define param_IN      prhs[4]           
            
    double *D, *X, *P0, *P, *G, gamma, tol;
    int M, N, K, N_G, nthreads, max_iter;
    
    if( nrhs < 4 ) /* Check number of arguments */
        mexErrMsgTxt("Wrong number of input arguments.");
    else if( nlhs > 1)
        mexErrMsgTxt("Wrong number of output arguments.");
    
    if( !IS_REAL_2D_FULL_DOUBLE(D_IN) ) /* Check D */
        mexErrMsgTxt("D must be a real 2D full double array.");
    else if( !IS_REAL_2D_FULL_DOUBLE(X_IN) ) /* Check X */
        mexErrMsgTxt("X must be a real 2D full double array.");
    else if( !IS_REAL_2D_FULL_DOUBLE(P0_IN) ) /* Check P0 */
        mexErrMsgTxt("P0 must be a real 2D full double array.");
    else if( !IS_REAL_2D_FULL_DOUBLE(G_IN) ) /* Check G */
        mexErrMsgTxt("G must be a real full double array.");
    
    /* Get the input variables */
    if( nrhs == 4)
    {
        gamma = 1;
        tol = 0.01;
        max_iter = 100;
#ifdef _OPENMP
        nthreads = omp_get_num_procs();
#else
        nthreads = 1;
#endif
    }
    else if( mxIsStruct(param_IN) )
    {
        if( mxGetFieldNumber(param_IN, "gamma") == -1 )
            gamma = 1;
        else
        {
            mxArray *gamma_IN = mxGetField(param_IN, 0, "gamma");
            if( !IS_REAL_SCALAR(gamma_IN) ) /* Check gamma */
                mexErrMsgTxt("gamma must be a real double scalar.");
            gamma = mxGetScalar(gamma_IN);
            if( gamma <= 0 )
                mexErrMsgTxt("gamma must be positive.");
        }
        if( mxGetFieldNumber(param_IN, "tol") == -1 )
            tol = 0.01;
        else
        {
            mxArray *tol_IN = mxGetField(param_IN, 0, "tol");
            if( !IS_REAL_SCALAR(tol_IN) ) /* Check tol */
                mexErrMsgTxt("tol must be a real double scalar.");
            tol = mxGetScalar(tol_IN);
            if( tol <= 0 )
                mexErrMsgTxt("tol must be positive.");
        }
        if( mxGetFieldNumber(param_IN, "max_iter") == -1 )
            max_iter = 100;
        else
        {
            mxArray *max_iter_IN = mxGetField(param_IN, 0, "max_iter");
            if( !IS_REAL_SCALAR(max_iter_IN) ) /* Check max_iter */
                mexErrMsgTxt("max_iter must be a real double scalar.");
            max_iter = mxGetScalar(max_iter_IN);
            if( max_iter <= 0 )
                mexErrMsgTxt("max_iter must positive.");
        }
        if( mxGetFieldNumber(param_IN, "nthreads") == -1 )
#ifdef _OPENMP
        nthreads = omp_get_num_procs();
#else
        nthreads = 1;
#endif
        else
        {
            mxArray *nthreads_IN = mxGetField(param_IN, 0, "nthreads");
            if( !IS_REAL_SCALAR(nthreads_IN) ) /* Check nthreads */
                mexErrMsgTxt("nthreads must be a real double scalar.");
            nthreads = mxGetScalar(nthreads_IN);
            if( nthreads <= 0 )
                mexErrMsgTxt("nthreads must positive.");
        }
    }
    else
        mexErrMsgTxt("param must be a structure array.");
    
    M = mxGetM(D_IN);  /* Get D */
    K = mxGetN(D_IN);
    D = mxGetPr(D_IN);
    
    N = mxGetN(X_IN); /* Get X */
    X = mxGetPr(X_IN);
    
    if( M != mxGetM(X_IN) ) 
        mexErrMsgTxt("Error regarding D and X: Inner matrix dimensions must agree.");
    
    
    N_G = mxGetM(G_IN) - 1; /* Get G */
    G = mxGetPr(G_IN);
    
    P0 = mxGetPr(P0_IN); /* Get P0 */
    if( K != mxGetM(P0_IN) || N != mxGetN(P0_IN) )
        mexErrMsgTxt("Error regarding P0 and X or D: Inner matrix dimensions must agree.");
    
    P_OUT = mxCreateDoubleMatrix( K, N, mxREAL ); /* Create Output Matrix */
    P = mxGetPr(P_OUT);
    
    if( (K==0) || (N==0) )
    {        
        my_blas_dcopy( K*N, P0, 1, P, 1 );
        return;
    }
    
    core(X, M, D, K, P0, P, G, N_G, gamma, tol, max_iter, nthreads );
}

void core(double *X, int M, double *D, int K, double *P0, double *P, double *G, int N_G, double gamma, double tol, double max_iter, int nthreads )
{
    
    int tid, sz = 0;
    for( int g = 0; g < N_G; g++ )
        if( sz < G[g+1] - G[g] )
            sz = G[g+1] - G[g];
#ifdef _OPENMP
    omp_set_num_threads(nthreads);
#else
    nthreads = 1;
#endif

    double * L = (double*)malloc(nthreads*K*sizeof(double));
    if( L==NULL )
    {
        printf( "You requested %d MB of memory for L.\n", nthreads*K*sizeof(double)/1024.0/1024.0 );
        mexErrMsgTxt("OUT OF MEMORY.");  
    }
    double * DL = (double*)malloc(nthreads*M*K*sizeof(double));
    if( DL==NULL )
    {
        printf( "You requested %d MB of memory for DL.\n", nthreads*M*K*sizeof(double)/1024.0/1024.0 );
        mexErrMsgTxt("OUT OF MEMORY.");  
    }
    double * DLDt = (double*)malloc(nthreads*M*M*sizeof(double));
    if( DLDt==NULL )
    {
        printf( "You requested %d MB of memory for DLDt.\n", nthreads*M*M*sizeof(double)/1024.0/1024.0 );
        mexErrMsgTxt("OUT OF MEMORY.");  
    }
    double * DLDt_inv_X = (double*)malloc(nthreads*M*sz*sizeof(double));
    if( DLDt_inv_X==NULL )
    {
        printf( "You requested %d MB of memory for DLDt_inv_X.\n", nthreads*M*sz*sizeof(double)/1024.0/1024.0 );
        mexErrMsgTxt("OUT OF MEMORY.");  
    }
    mwSignedIndex *ipiv = (mwSignedIndex*)malloc(nthreads*M*sizeof(mwSignedIndex));
    if( ipiv==NULL )
    {
        printf( "You requested %d MB of memory for ipiv.\n", nthreads*M*sizeof(double)/1024.0/1024.0 );
        mexErrMsgTxt("OUT OF MEMORY.");  
    }
    
//     if( L==NULL || DL==NULL || DLDt==NULL || DLDt_inv_X==NULL || ipiv==NULL )
//         mexErrMsgTxt("OUT OF MEMORY.");
    
    #pragma omp parallel shared(X, M, D, K, P0, P, G, N_G, gamma, tol, max_iter, L, DL, DLDt, DLDt_inv_X) private(nthreads,tid)
    {
#ifdef _OPENMP
        tid = omp_get_thread_num();
        nthreads = omp_get_num_threads();
#else
        tid = 0;
        nthreads = 1;
#endif        
        mfoccus_core(X, M, D, K, P0, P, G, N_G, gamma, tol, max_iter, L + tid*K, DL + tid*M*K, DLDt + tid*M*M, DLDt_inv_X + tid*M*sz, ipiv + tid*M, tid, nthreads);
    }

    free((void*)DLDt_inv_X);
    free((void*)DL);
    free((void*)DLDt);
    free((void*)L);
    free((void*)ipiv);
}

void mfoccus_core(double *X, int M, double *D, int K, double *P0, double *P, double *G, int N_G, double gamma, double tol, double max_iter, double *L, double *DL, double *DLDt, double *DLDt_inv_X, mwSignedIndex *ipiv, int tid, int nthreads)
{
    int N;
    double *mX, *mP, *mP0;
    for( int g = tid; g < N_G; g+=nthreads )
    {
        N = G[g];
        mX = X + N*M;
        mP = P + N*K;
        mP0 = P0 + N*K;
        
        N = G[g+1] - G[g];
        double new_gamma = gamma*(sqrt(G[g+1] - G[g])+1/sqrt(G[g+1] - G[g]));
        
        double delta = 1;
        int iter = 1;
        while( (delta > tol) && (iter <= max_iter) )
        {
            my_blas_dcopy(M*K, D, 1, DL, 1);   
            for( int i = 0; i < K; i++ )
            {
                L[i] = my_blas_dnrm2(N,mP0+i,K);
                my_blas_dscal( M, L[i], DL+i*M, 1);
            }

            my_blas_dgemm ( 'N', 'T', M, M, K, 1, DL, M, D, M, 0, DLDt, M);
            for( int i = 0; i < M*M; i+=M+1)
                DLDt[i] += new_gamma;
            my_blas_dcopy(M*N,mX,1,DLDt_inv_X,1);
            

            //put the new stuff into DLLDt_inv_X = DLLDt\X
            if( my_lapack_dgesv( M, N, DLDt, M, ipiv, DLDt_inv_X, M) > 0 )
                mexErrMsgTxt("Singular matrix encountered");

            my_blas_dgemm ( 'T', 'N', K, N, M, 1, DL, M, DLDt_inv_X, M, 0, mP, K);

            my_blas_daxpy(K*N, -1.0, mP, 1, mP0, 1);
            delta = my_blas_ddot(K*N, mP0, 1, mP0, 1)/my_blas_ddot(K, L, 1, L, 1);
            // copy P into P0
            my_blas_dcopy( K*N, mP, 1, mP0, 1 );
        }
    }
}

void my_blas_dgemm ( const char TransA, const char TransB, const int M, const int N, const int K, const double ALPHA, const double * A, const int LDA, const double * B, const int LDB, const double BETA, double * C, const int LDC)
{
    mwSignedIndex m = M;
    mwSignedIndex n = N;
    mwSignedIndex k = K;
    mwSignedIndex lda = LDA;
    mwSignedIndex ldb = LDB;
    mwSignedIndex ldc = LDC;

    double alpha = ALPHA;
    double beta = BETA;
    
    char * trA = (char*)&TransA;
    char * trB = (char*)&TransB;
    
    double *a = (double*)A;
    double *b = (double*)B;
    
    dgemm(trA, trB, &m, &n, &k, &alpha, a, &lda, b, &ldb, &beta, C, &ldc);
}

double my_blas_dnrm2(int N,double *X,int INCX)
{
    mwSignedIndex n = N;
    mwSignedIndex incX = INCX;
    return dnrm2( &n, X, &incX );
}

void my_blas_dcopy(int N, double *X, int INCX, double* Y, int INCY )
{
    mwSignedIndex n = N;
    mwSignedIndex incX = INCX;
    mwSignedIndex incY = INCY;
    dcopy( &n, X, &incX, Y, &incY );
}

void my_blas_dscal( int N, double ALPHA, double * X, int INCX )
{
    mwSignedIndex n = N;
    mwSignedIndex incX = INCX;
    dscal( &n, &ALPHA, X, &incX );
}

double my_blas_ddot( int N, double *X, int INCX, double *Y, int INCY)
{
    mwSignedIndex n = N;
    mwSignedIndex incX = INCX;
    mwSignedIndex incY = INCY;
    return ddot( &n, X, &incX, Y, &incY);
}

void my_blas_daxpy( int N, double alpha, double *X, int INCX, double *Y, int INCY)
{
    mwSignedIndex n = N;
    mwSignedIndex incX = INCX;
    mwSignedIndex incY = INCY;
    daxpy( &n, &alpha, X, &incX, Y, &incY);
}

int my_lapack_dgesv( int N, int NRHS, double *A, int LDA, mwSignedIndex *ipiv, double *B, int LDB)
{
    mwSignedIndex info;
    mwSignedIndex n = N;
    mwSignedIndex nrhs = NRHS;
    mwSignedIndex lda = LDA;
    mwSignedIndex ldb = LDB;
    dgesv( &n, &nrhs, A, &lda, ipiv, B, &ldb, &info );
    return info;
}
