
/******************************************************************************

  Copyright (C) 2011 by Lorenzo Stella <lorenzo.stella77@gmail.com>

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

#ifndef __PolyCEID_BLAS__
#define __PolyCEID_BLAS__

#include "PolyCEID_complex.h"


#define F77_FUNC( x, y ) x##_ //WARNING! it must be defined in a safer way !


complex F77_FUNC( zdotc, ZDOTC)( int *n, double *x, int *incx, double *y, int *incy );


void    F77_FUNC( zgemv, ZGEMV )( char *trans, int *m, int *n, double *alpha, 
				  double *a, int *lda, double *x, int *incx, double *beta, 
				  double *y, int *incy );

void    F77_FUNC( zgemm, ZGEMM )( char *transa, char *transb, int *m, int *n, int *k,
				  double *alpha, double *a, int *lda, double *b, int *ldb, double *beta, 
				  double *c, int *ldc );

void    F77_FUNC( daxpy, DAXPY )( int *n, double *alpha, double *x, int *incx, double *y, int *incy );

double  F77_FUNC( zlange, ZLANGE )( char *norm, int *m, int *n, double *a, int *lda, double* work );

double  F77_FUNC( dznrm2, DZNRM2 )( int *n, double *x, int *incx );

double  F77_FUNC( dnrm2,  DNRM2  )( int *n, double *x, int *incx );

#endif /* __PolyCEID_BLAS__ */

