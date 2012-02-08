
/******************************************************************************

  Copyright (C) 2011-2012 by Lorenzo Stella <lorenzo DOT stella77 AT gmail.com>

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in
  all copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
  THE SOFTWARE.

******************************************************************************/

#ifndef __PolyCEID_LAPACK__
#define __PolyCEID_LAPACK__

#include "PolyCEID_complex.h"
#include "PolyCEID_blas.h"


void F77_FUNC( zheev, ZHEEV )( char *jobz, char *uplo, int *n, double *a, 
			       int *lda, double *w, double *work, int *lwork, double *rwork, int *info);

void F77_FUNC( zgeev, ZGEEV )( char *jobvl, char *jobvr, int *n, double *a, int *lda, double *w, double *vl,  
			       int *ldvl, double *vr, int *ldvr, double *work, int *lwork, double *rwork, 
			       int *info );
/*
void F77_FUNC( zgesvd, ZGESVD)( char *jobu, char *jobvt, int *m, int *n, double *a, int *lda, double *s, 
				double *u, int *ldu, double *vt, int *ldvt, double *work, int *lwork, 
				double *rwork, int *info );
*/

void F77_FUNC( zgetrf, ZGETRF )( int *n, int *m, double *a, int *lda, int *ipiv, int *info);

void F77_FUNC( zgetri, ZGETRI )( int *n, double *a, int *lda, int *ipiv, double *work, int *lwork, int *info);

#endif /* __PolyCEID_LAPACK__ */
