
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

#ifndef __PolyCEID_COMPLEX__
#define __PolyCEID_COMPLEX__ 

#include "../algebra/PolyCEID_constants.h"
#include "../utils/my_error.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>


typedef struct{

  double z[2]; 

} __my_complex;

typedef __my_complex   complex;
typedef __my_complex*  complex_p;
typedef __my_complex** complex_pp;


complex __my_complex_cast( const double );

#define CMPLX( z ) __my_complex_cast( (z) )

complex __my_complex_initialize( const double, const double );

#define CMPLX_INIT( x, y ) __my_complex_initialize( (x), (y) )

double  __my_complex_real_part( const complex );

#define REAL( z ) __my_complex_real_part( (z) )

double  __my_complex_imaginary_part( const complex );

#define IMAG( z ) __my_complex_imaginary_part( (z) )

/*
double  __my_complex_norm( const complex );

  #define CMPLX_NORM( z )  __my_complex_norm( (z) )
*/

complex __my_complex_times_I( const complex );

#define CMPLX_TIMES_I( z ) __my_complex_times_I( (z) )

complex __my_complex_conjugate( const complex );

#define CONJ( z ) __my_complex_conjugate( (z) )

complex __my_complex_inverse( const complex );

#define CMPLX_INVERSE( z )  __my_complex_inverse( (z) )

complex __my_complex_sum( const complex, const complex );

#define CMPLX_SUM( z1, z2 ) __my_complex_sum( (z1), (z2) )

complex __my_complex_dif( const complex, const complex );

#define CMPLX_DIF( z1, z2 ) __my_complex_dif( (z1), (z2) )

complex __my_complex_product( const complex, const complex );

#define CMPLX_PRODUCT( z1, z2 ) __my_complex_product( (z1), (z2) )

complex __my_complex_division( const complex, const complex );

#define CMPLX_DIVISION(z1, z2 ) __my_complex_division( (z1), (z2) )

complex __my_complex_phase_fix( const complex );

#define CMPLX_PHASE_FIX( z ) __my_complex_phase_fix( (z) )

int     __my_complex_copy( complex_p, const complex );

#define CMPLX_COPY( z1, z2 )    FUNCTION_CHECK( __my_complex_copy( &(z1), (z2) ), __my_complex_copy )

int     __my_complex_compare( const complex, const complex );

#define CMPLX_COMPARE( z1, z2 )    FUNCTION_CHECK( __my_complex_compare( (z1), (z2) ), __my_complex_compare )

int     __my_complex_read( FILE*, complex_p );

#define CMPLX_READ( x, y ) FUNCTION_CHECK( __my_complex_read( (x), &(y) ),  __my_complex_read )

int     __my_complex_print( FILE*, const complex );

#define CMPLX_PRINT( x, y ) FUNCTION_CHECK( __my_complex_print( (x), (y) ),  __my_complex_print )

int     __my_complex_print_plus( FILE*, const complex );

#define CMPLX_PRINT_PLUS( x, y ) FUNCTION_CHECK( __my_complex_print_plus( (x), (y) ),  __my_complex_print_plus )

int     __my_complex_print_plus_var( FILE*, const complex );

#define CMPLX_PRINT_PLUS_VAR( x, y ) FUNCTION_CHECK( __my_complex_print_plus_var( (x), (y) ),  __my_complex_print_plus_var )


/* constant */

#define CMPLX_ZERO              CMPLX_INIT( 0.0e0, 0.0e0 )

#define CMPLX_ZERO_VAR( c )     (c).z[0]=0.0e0; (c).z[1]=0.0e0    //WARNING: faster

#define CMPLX_UNIT              CMPLX_INIT( 1.0e0, 0.0e0 )

#define CMPLX_UNIT_VAR( c )     (c).z[0]=1.0e0; (c).z[1]=0.0e0    //WARNING: faster

#define CMPLX_I                 CMPLX_INIT( 0.0e0, 1.0e0 )

#define CMPLX_I_VAR( c )        (c).z[0]=0.0e0; (c).z[1]=1.0e0    //WARNING: faster

#define CMPLX_HUGE_VAL          CMPLX( HUGE_VAL )

#define CMPLX_HUGE_VAL_VAR( c ) (c).z[0]=HUGE_VAL; (c).z[1]=0.0e0 //WARNING: faster

#endif /*  __PolyCEID_COMPLEX__ */

