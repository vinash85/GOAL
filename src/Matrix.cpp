#include "Matrix.h"

////////////////////////////////////////////////////////////////////////////////

void rsyevd(char jobz, char uplo, int n, double* a, int lda, double* w, double* work, int lwork, int* iwork, int liwork, int* info)
{ return dsyevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, info); }

void rgesvd(char jobu, char jobvt, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt, double*work, int lwork, int* info)
{ return dgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, info); }

void rgesdd(char jobz, int m, int n, double* a, int lda, double* s, double* u, int ldu, double* vt, int ldvt, double* work, int lwork, int* iwork, int* info)
{ return dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info); }

////////////////////////////////////////////////////////////////////////////////

#ifndef DISABLE_SINGLE

void rsyevd(char jobz, char uplo, int n, float* a, int lda, float* w, float* work, int lwork, int* iwork, int liwork, int* info)
{ return ssyevd_(&jobz, &uplo, &n, a, &lda, w, work, &lwork, iwork, &liwork, info); }

void rgesvd(char jobu, char jobvt, int m, int n, float* a, int lda, float* s, float* u, int ldu, float* vt, int ldvt, float*work, int lwork, int* info)
{ return sgesvd_(&jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, info); }

void rgesdd(char jobz, int m, int n, float* a, int lda, float* s, float* u, int ldu, float* vt, int ldvt, float* work, int lwork, int* iwork, int* info)
{ return sgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt, work, &lwork, iwork, info); }

#endif

////////////////////////////////////////////////////////////////////////////////
