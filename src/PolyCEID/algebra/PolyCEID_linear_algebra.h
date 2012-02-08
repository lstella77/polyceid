
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

#ifndef __PolyCEID_LINEAR_ALGEBRA__
#define __PolyCEID_LINEAR_ALGEBRA__

#include "PolyCEID_lapack.h"
#include "PolyCEID_blas.h"
#include "PolyCEID_vector.h"
#include "PolyCEID_matrix.h"


#define MAX( a, b )  ( (a) > (b) ) ? (a) : (b)
#define MIN( a, b )  ( (a) < (b) ) ? (a) : (b)


/*----------------------------------------- */

/* private working space utilities */

/*----------------------------------------- */

int       __my_working_space_check( const int );

#define WORKING_SPACE_CHECK( n ) FUNCTION_CHECK(  __my_working_space_check( (n) ),  __my_working_space_check  )

/*----------------------------------------- */

int       __my_vector_transcription( const vector, double* );

#define VECTOR_TRANSCRIPTION( a, ap ) FUNCTION_CHECK( __my_vector_transcription( (a), (ap) ), __my_vector_transcription )

/*----------------------------------------- */

int       __my_matrix_transcription( const matrix, double* );

#define MATRIX_TRANSCRIPTION( a, ap ) FUNCTION_CHECK( __my_matrix_transcription( (a), (ap) ), __my_matrix_transcription )

/*----------------------------------------- */

int       __my_vector_pointer_check( const vector, const vector);

#define VECTOR_POINTER_CHECK( v, w ) FUNCTION_CHECK(  __my_vector_pointer_check( (v), (w) ), __my_vector_pointer_check )

/*----------------------------------------- */

int       __my_matrix_pointer_check( const matrix, const matrix);

#define MATRIX_POINTER_CHECK( a, b ) FUNCTION_CHECK(  __my_matrix_pointer_check( (a), (b) ), __my_matrix_pointer_check )

/*----------------------------------------- */

/* utilities */

/*----------------------------------------- */

int     __my_matrix_vector_product_consinstency( const matrix, const vector, const vector );

#define MATRIX_VECTOR_PRODUCT_CONSINSTENCY( a, v, w )  FUNCTION_CHECK( __my_matrix_vector_product_consinstency( (a), (v), (w) ), \
									       __my_matrix_vector_product_consinstency )
/*----------------------------------------- */

int     __my_matrix_matrix_product_consinstency( const matrix, const matrix, const matrix );

#define MATRIX_MATRIX_PRODUCT_CONSINSTENCY( a, b, c )  FUNCTION_CHECK( __my_matrix_matrix_product_consinstency( (a), (b), (b) ), \
									       __my_matrix_matrix_product_consinstency )
/*----------------------------------------- */

/* blas & lapack */

/*----------------------------------------- */

double  __my_complex_norm( const complex );

#define CMPLX_NORM( z )  __my_complex_norm( (z) )

/*----------------------------------------- */

complex __my_scalar_product( const vector, const vector );

#define SCALAR_PRODUCT( a, b ) __my_scalar_product( (a), (b) )

/*----------------------------------------- */

int     __my_matrix_vector_product( const matrix, const vector, vector_p );

#define MATRIX_VECTOR_PRODUCT( a, v, w ) FUNCTION_CHECK( __my_matrix_vector_product( (a), (v), &(w) ),  __my_matrix_vector_product )

/*----------------------------------------- */

int     __my_matrix_matrix_product( const matrix, const matrix, matrix_p );

#define MATRIX_MATRIX_PRODUCT( a, b, c ) FUNCTION_CHECK( __my_matrix_matrix_product( (a), (b), &(c) ),  __my_matrix_matrix_product )
 
/*----------------------------------------- */

double  __my_vector_norm( const vector );

#define VECTOR_NORM( v ) __my_vector_norm( (v) )

/*----------------------------------------- */

double  __my_rvector_norm( const rvector );

#define RVECTOR_NORM( v ) __my_rvector_norm( (v) )

/*----------------------------------------- */

double  __my_matrix_norm( const matrix );

#define MATRIX_NORM( a ) __my_matrix_norm( (a) )

/*----------------------------------------- */

int     __my_diagonalisation( const matrix, matrix_p, rvector_p );

#define DIAGONALISATION( a, vr, w ) FUNCTION_CHECK( __my_diagonalisation( (a), &(vr), &(w) ), \
									    __my_diagonalisation ) 
/*----------------------------------------- */

int     __my_non_symmetric_diagonalisation( const matrix, matrix_p, matrix_p, vector_p );

#define NON_SYMMETRIC_DIAGONALISATION( a, vl, vr, w ) FUNCTION_CHECK( __my_non_symmetric_diagonalisation( (a), &(vl), &(vr), &(w) ), \
									    __my_non_symmetric_diagonalisation ) 
/*----------------------------------------- */

int     __my_matrix_inverse( const matrix, matrix_p );

#define MATRIX_INVERSE( mat1, mat2 ) FUNCTION_CHECK( __my_matrix_inverse( (mat1), &(mat2) ), __my_matrix_inverse )

/*----------------------------------------- */

int     __my_rvector_increase( rvector_p, const rvector, double );

#define RVECTOR_INCREASE( vec1, vec2, alpha ) FUNCTION_CHECK( __my_rvector_increase( &(vec1), (vec2), (alpha) ), __my_rvector_increase ) 

/*----------------------------------------- */

/* non-blas & non-lapack functions */

/*----------------------------------------- */

int     __my_matrix_adjoint( const matrix, matrix_p );

#define MATRIX_ADJOINT( mat1, mat2 ) FUNCTION_CHECK( __my_matrix_adjoint( (mat1), &(mat2) ), __my_matrix_adjoint )

/*----------------------------------------- */

int     __my_matrix_conjugate( const matrix, matrix_p );

#define MATRIX_CONJ( mat1, mat2 ) FUNCTION_CHECK( __my_matrix_conjugate( (mat1), &(mat2) ), __my_matrix_conjugate )

/*----------------------------------------- */

int     __my_matrix_transpose( const matrix, matrix_p );

#define MATRIX_TRANS( mat1, mat2 ) FUNCTION_CHECK( __my_matrix_transpose( (mat1), &(mat2) ), __my_matrix_transpose )

/*----------------------------------------- */

complex __my_matrix_trace( const matrix );

#define MATRIX_TRACE( mat ) __my_matrix_trace( (mat) )

/*----------------------------------------- */

int     __my_matrix_scalar_product( const matrix, const complex, matrix_p );

#define  MATRIX_SCALAR_PRODUCT( a, z, b ) FUNCTION_CHECK(__my_matrix_scalar_product( (a), (z), &(b) ), __my_matrix_scalar_product )

/*----------------------------------------- */

int     __my_matrix_scalar_division( const matrix, const complex, matrix_p );

#define  MATRIX_SCALAR_DIVISION( a, z, b ) FUNCTION_CHECK(__my_matrix_scalar_division( (a), (z), &(b) ), __my_matrix_scalar_division )

/*----------------------------------------- */

int     __my_matrix_addition( const matrix, const matrix, matrix_p );

#define MATRIX_SUM( a, b, c ) FUNCTION_CHECK(  __my_matrix_addition( (a), (b), &(c) ),  __my_matrix_addition )

/*----------------------------------------- */

int     __my_matrix_subtraction( const matrix, const matrix, matrix_p );

#define MATRIX_DIF( a, b, c ) FUNCTION_CHECK(  __my_matrix_subtraction( (a), (b), &(c) ),  __my_matrix_subtraction )

/*----------------------------------------- */

int     __my_matrix_commutator( const matrix, const matrix, matrix_p );

int     __my_matrix_commutator_alt( const matrix, const matrix, matrix_p );

#define MATRIX_COMMUTATOR( a, b, c ) FUNCTION_CHECK( __my_matrix_commutator_alt( (a), (b), &(c) ), __my_matrix_commutator_alt )

/*----------------------------------------- */

int     __my_matrix_anticommutator( const matrix, const matrix, matrix_p );

int     __my_matrix_anticommutator_alt( const matrix, const matrix, matrix_p );

#define MATRIX_ANTICOMMUTATOR( a, b, c ) FUNCTION_CHECK( __my_matrix_anticommutator_alt( (a), (b), &(c) ), __my_matrix_anticommutator_alt )

/*----------------------------------------- */

int     __my_matrix_commutator_hermitian( const matrix, const matrix, matrix_p );

#define MATRIX_COMMUTATOR_HERMITIAN( a, b, c ) FUNCTION_CHECK( __my_matrix_commutator_hermitian( (a), (b), &(c) ), __my_matrix_commutator_hermitian )

/*----------------------------------------- */

int     __my_matrix_anticommutator_hermitian( const matrix, const matrix, matrix_p );

#define MATRIX_ANTICOMMUTATOR_HERMITIAN( a, b, c ) FUNCTION_CHECK( __my_matrix_anticommutator_hermitian( (a), (b), &(c) ), __my_matrix_anticommutator_hermitian )

/*----------------------------------------- */

#endif /* __PolyCEID_LINEAR_ALGEBRA__ */

