
/******************************************************************************

  Copyright (C) 2011 by Lorenzo Stella <lorenzo DOT stella77 AT gmail DOT com>

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

#ifndef __PolyCEID_MATRIX__
#define __PolyCEID_MATRIX__


#include "../algebra/PolyCEID_complex.h"
#include "../utils/my_error.h"
#include <stdlib.h>


/* matrix */

typedef struct{

  int                matrix_dim_row;
  int                matrix_dim_column;
  int                matrix_dim;
  complex*           matrix;
  complex**          column_p;

} __my_matrix;

typedef __my_matrix   matrix;
typedef __my_matrix*  matrix_p;
typedef __my_matrix** matrix_pp;

/* matrix_array */

typedef struct{

  int      dim;
  matrix* array; 

} __my_matrix_array;

typedef __my_matrix_array   matrix_array;
typedef __my_matrix_array*  matrix_array_p;

/* rmatrix */

typedef struct{

  int                rmatrix_dim_row;
  int                rmatrix_dim_column;
  int                rmatrix_dim;
  double*            rmatrix;
  double**           rcolumn_p;

} __my_rmatrix;

typedef __my_rmatrix   rmatrix;
typedef __my_rmatrix*  rmatrix_p;
typedef __my_rmatrix** rmatrix_pp;

/* imatrix */

typedef struct{

  int                imatrix_dim_row;
  int                imatrix_dim_column;
  int**              imatrix;

} __my_imatrix;

typedef __my_imatrix    imatrix;
typedef __my_imatrix*   imatrix_p;


int __my_matrix_allocate( const int, const int, matrix_p );

#define MATRIX_ALLOCATE( n, m, mat ) FUNCTION_CHECK( __my_matrix_allocate( (n), (m), &(mat) ), __my_matrix_allocate )

int       __my_matrix_free( matrix_p );

#define MATRIX_FREE( mat )        FUNCTION_CHECK( __my_matrix_free( &(mat) ),  __my_matrix_free )

int       __my_matrix_copy( matrix_p, const matrix );

#define MATRIX_COPY( mat1, mat2 )      FUNCTION_CHECK( __my_matrix_copy( &(mat1), (mat2) ), __my_matrix_copy )

int       __my_matrix_compare( const matrix, const matrix );

#define MATRIX_COMPARE( mat1, mat2 )      FUNCTION_CHECK( __my_matrix_compare( (mat1), (mat2) ), __my_matrix_compare )

int       __my_matrix_read( FILE*, matrix_p );

#define MATRIX_READ( x, y )           FUNCTION_CHECK( __my_matrix_read( (x), &(y) ), __my_matrix_read )

int       __my_matrix_print( FILE*, const matrix );

#define MATRIX_PRINT( x, y )           FUNCTION_CHECK( __my_matrix_print( (x), (y) ), __my_matrix_print )

int       __my_matrix_print_plus( FILE*, const matrix );

#define MATRIX_PRINT_PLUS( x, y )           FUNCTION_CHECK( __my_matrix_print_plus( (x), (y) ), __my_matrix_print_plus )

/* conversion to a vectorised form */  

int __my_matrix_arg( const int, const int, const int, const int );

#ifdef __DEBUG__

#define MATRIX_ARG( i, j ,n ,m )  __my_matrix_arg( (i), (j), (n), (m) ) // WARNING! Do not add FUNCTION_CHECK, it must return -1 !

#else /* __DEBUG__ */

#define MATRIX_ARG( i, j ,n ,m ) ( (i)<0 || (i)>(n)-1 || (j)<0 || (j)>(m)-1 ) ? -1 : (n) *(j) +(i)

#endif /* __DEBUG__ */

/* conversion from a vectorised form */  

int __my_matrix_arg_inverse( const int, int*, int*, const int, const int);

#define MATRIX_ARG_INVERSE( index, i, j ,n ,m ) FUNCTION_CHECK( __my_matrix_arg_inverse( (index), &(i), &(j), (n), (m) ), __my_matrix_arg_inverse )

/* Errors handling */

int       __my_matrix_check( const matrix );

#define MATRIX_CHECK( a )            FUNCTION_CHECK( __my_matrix_check( (a) ), __my_matrix_check )

int       __my_matrix_consinstency( const matrix, const matrix );

#define MATRIX_CONSINSTENCY( a, b )  FUNCTION_CHECK( __my_matrix_consinstency( (a), (b) ), __my_matrix_consinstency )

int       __my_matrix_square_check( const matrix );

#define MATRIX_SQUARE_CHECK( a )     FUNCTION_CHECK( __my_matrix_square_check( (a) ), __my_matrix_square_check )

/* notable matrices */

int    __my_matrix_set_to_zero( matrix_p );

#define MATRIX_ZERO( mat )                   FUNCTION_CHECK( __my_matrix_set_to_zero( &(mat) ), __my_matrix_set_to_zero )

int    __my_matrix_set_to_unit( matrix_p ); 

#define MATRIX_UNIT( mat )                   FUNCTION_CHECK( __my_matrix_set_to_unit( &(mat) ), __my_matrix_set_to_unit )

int    __my_matrix_set_to_sigma_x( matrix_p ); 

#define MATRIX_SIGMA_X( mat )                FUNCTION_CHECK( __my_matrix_set_to_sigma_x( &(mat) ), __my_matrix_set_to_sigma_x )

int    __my_matrix_set_to_sigma_y( matrix_p ); 

#define MATRIX_SIGMA_Y( mat )                FUNCTION_CHECK( __my_matrix_set_to_sigma_y( &(mat) ), __my_matrix_set_to_sigma_y )

int    __my_matrix_set_to_sigma_z( matrix_p ); 

#define MATRIX_SIGMA_Z( mat )                FUNCTION_CHECK( __my_matrix_set_to_sigma_z( &(mat) ), __my_matrix_set_to_sigma_z )

/* matrix_array */

int     __my_matrix_array_allocate( int, int, int, matrix_array_p );

#define MATRIX_ARRAY_ALLOCATE( n, dim_row, dim_column, az )   FUNCTION_CHECK(  __my_matrix_array_allocate( (n), (dim_row), (dim_column), &(az) ), __my_matrix_array_allocate )

int     __my_matrix_array_free( matrix_array_p );

#define MATRIX_ARRAY_FREE( az )   FUNCTION_CHECK(  __my_matrix_array_free( &(az) ),  __my_matrix_array_free )

int     __my_matrix_array_copy( matrix_array_p, const matrix_array );

#define MATRIX_ARRAY_COPY( az1, az2 )   FUNCTION_CHECK(  __my_matrix_array_copy( &(az1), (az2) ),  __my_matrix_array_copy )

int     __my_matrix_array_compare( const matrix_array, const matrix_array );

#define MATRIX_ARRAY_COMPARE( az1, az2 )   FUNCTION_CHECK(  __my_matrix_array_compare( (az1), (az2) ),  __my_matrix_array_compare )

int     __my_matrix_array_read( FILE*, matrix_array_p );

#define MATRIX_ARRAY_READ( fp, az )   FUNCTION_CHECK(  __my_matrix_array_read( fp, &(az) ),  __my_matrix_array_read )

int     __my_matrix_array_print( FILE*, const matrix_array );

#define MATRIX_ARRAY_PRINT( fp, az )   FUNCTION_CHECK(  __my_matrix_array_print( fp, (az) ),  __my_matrix_array_print )

int     __my_matrix_array_print_plus( FILE*, const matrix_array );

#define MATRIX_ARRAY_PRINT_PLUS( fp, az )   FUNCTION_CHECK(  __my_matrix_array_print_plus( fp, (az) ),  __my_matrix_array_print_plus )

int     __my_matrix_array_zero( matrix_array_p );

#define MATRIX_ARRAY_ZERO( az )   FUNCTION_CHECK(  __my_matrix_array_zero( &(az) ),  __my_matrix_array_zero )

/* imatrix */

int     __my_imatrix_check( const imatrix );

#define IMATRIX_CHECK( a )            FUNCTION_CHECK( __my_imatrix_check( (a) ), __my_imatrix_check )

int       __my_imatrix_consinstency( const imatrix, const imatrix );

#define IMATRIX_CONSINSTENCY( a, b )  FUNCTION_CHECK( __my_imatrix_consinstency( (a), (b) ), __my_imatrix_consinstency )

int     __my_imatrix_allocate( const int, const int, imatrix_p );

#define IMATRIX_ALLOCATE( n, m, mat ) FUNCTION_CHECK( __my_imatrix_allocate( (n), (m), &(mat) ), __my_imatrix_allocate )

int     __my_imatrix_free( imatrix_p );

#define IMATRIX_FREE( mat )        FUNCTION_CHECK( __my_imatrix_free( &(mat) ),  __my_imatrix_free )

int       __my_imatrix_copy( imatrix_p, const imatrix );

#define IMATRIX_COPY( mat1, mat2 )      FUNCTION_CHECK( __my_imatrix_copy( &(mat1), (mat2) ), __my_imatrix_copy )

int       __my_imatrix_compare( const imatrix, const imatrix );

#define IMATRIX_COMPARE( mat1, mat2 )      FUNCTION_CHECK( __my_imatrix_compare( (mat1), (mat2) ), __my_imatrix_compare )

int       __my_imatrix_read( FILE*, imatrix_p );

#define IMATRIX_READ( x, y )           FUNCTION_CHECK( __my_imatrix_read( (x), &(y) ), __my_imatrix_read )

int       __my_imatrix_print( FILE*, const imatrix );

#define IMATRIX_PRINT( x, y )           FUNCTION_CHECK( __my_imatrix_print( (x), (y) ), __my_imatrix_print )

int       __my_imatrix_print_plus( FILE*, const imatrix );

#define IMATRIX_PRINT_PLUS( x, y )           FUNCTION_CHECK( __my_imatrix_print_plus( (x), (y) ), __my_imatrix_print_plus )

int    __my_imatrix_set_to_zero( imatrix_p );

#define IMATRIX_ZERO( mat )                   FUNCTION_CHECK( __my_imatrix_set_to_zero( &(mat) ), __my_imatrix_set_to_zero )

/* rmatrix */

int     __my_rmatrix_check( const rmatrix );

#define RMATRIX_CHECK( a )            FUNCTION_CHECK( __my_rmatrix_check( (a) ), __my_rmatrix_check )

int       __my_rmatrix_consinstency( const rmatrix, const rmatrix );

#define RMATRIX_CONSINSTENCY( a, b )  FUNCTION_CHECK( __my_rmatrix_consinstency( (a), (b) ), __my_rmatrix_consinstency )

int     __my_rmatrix_allocate( const int, const int, rmatrix_p );

#define RMATRIX_ALLOCATE( n, m, mat ) FUNCTION_CHECK( __my_rmatrix_allocate( (n), (m), &(mat) ), __my_rmatrix_allocate )

int     __my_rmatrix_free( rmatrix_p );

#define RMATRIX_FREE( mat )        FUNCTION_CHECK( __my_rmatrix_free( &(mat) ),  __my_rmatrix_free )

int       __my_rmatrix_copy( rmatrix_p, const rmatrix );

#define RMATRIX_COPY( mat1, mat2 )      FUNCTION_CHECK( __my_rmatrix_copy( &(mat1), (mat2) ), __my_rmatrix_copy )

int       __my_rmatrix_compare( const rmatrix, const rmatrix );

#define RMATRIX_COMPARE( mat1, mat2 )      FUNCTION_CHECK( __my_rmatrix_compare( (mat1), (mat2) ), __my_rmatrix_compare )

int       __my_rmatrix_read( FILE*, rmatrix_p );

#define RMATRIX_READ( x, y )           FUNCTION_CHECK( __my_rmatrix_read( (x), &(y) ), __my_rmatrix_read )

int       __my_rmatrix_print( FILE*, const rmatrix );

#define RMATRIX_PRINT( x, y )           FUNCTION_CHECK( __my_rmatrix_print( (x), (y) ), __my_rmatrix_print )

int       __my_rmatrix_print_plus( FILE*, const rmatrix );

#define RMATRIX_PRINT_PLUS( x, y )           FUNCTION_CHECK( __my_rmatrix_print_plus( (x), (y) ), __my_rmatrix_print_plus )

int    __my_rmatrix_set_to_zero( rmatrix_p );

#define RMATRIX_ZERO( mat )                   FUNCTION_CHECK( __my_rmatrix_set_to_zero( &(mat) ), __my_rmatrix_set_to_zero )

#endif /* __PolyCEID_MATRIX__ */
