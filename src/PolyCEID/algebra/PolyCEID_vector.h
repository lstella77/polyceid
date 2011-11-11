
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

#ifndef __PolyCEID_VECTOR__
#define __PolyCEID_VECTOR__

#include "../utils/my_error.h"
#include "../algebra/PolyCEID_complex.h"
#include <stdlib.h>

/* complex vectors */

typedef struct{

  int                vector_dim;
  complex*           vector;

} __my_vector;

typedef __my_vector  vector;
typedef __my_vector* vector_p;

typedef struct{

  int      dim;
  vector*  array; 

} __my_vector_array;

typedef __my_vector_array   vector_array;
typedef __my_vector_array*  vector_array_p;


// vector

int __my_vector_allocate( const int, vector_p );

#define VECTOR_ALLOCATE( n, vec ) FUNCTION_CHECK( __my_vector_allocate( (n), &(vec) ), __my_vector_allocate )

int       __my_vector_free( vector_p );

#define VECTOR_FREE( vec )     FUNCTION_CHECK( __my_vector_free( &(vec) ),  __my_vector_free )

int       __my_vector_copy( vector_p, const vector );

#define VECTOR_COPY( vec1, vec2 )   FUNCTION_CHECK(  __my_vector_copy( &(vec1), (vec2) ), __my_vector_copy )

int       __my_vector_compare( const vector, const vector );

#define VECTOR_COMPARE( vec1, vec2 )   FUNCTION_CHECK(  __my_vector_compare( (vec1), (vec2) ), __my_vector_compare )

int       __my_vector_read( FILE*, const vector );

#define VECTOR_READ( x, y )        FUNCTION_CHECK( __my_vector_read( (x), (y) ), __my_vector_read )

int       __my_vector_print( FILE*, const vector );

#define VECTOR_PRINT( x, y )        FUNCTION_CHECK( __my_vector_print( (x), (y) ), __my_vector_print )

int       __my_vector_print_plus( FILE*, const vector );

#define VECTOR_PRINT_PLUS( x, y )        FUNCTION_CHECK( __my_vector_print_plus( (x), (y) ), __my_vector_print_plus )

/* Errors handling */

int       __my_vector_check( const vector );

#define VECTOR_CHECK( a )   FUNCTION_CHECK( __my_vector_check( (a) ), __my_vector_check )

int       __my_vector_consinstency( const vector, const vector );

#define VECTOR_CONSINSTENCY( a, b )  FUNCTION_CHECK( __my_vector_consinstency( (a), (b) ), __my_vector_consinstency )

int    __my_vector_set_to_zero( vector_p );

#define VECTOR_ZERO( vec )                   FUNCTION_CHECK( __my_vector_set_to_zero( &(vec) ), __my_vector_set_to_zero )

int    __my_vector_set_to_unit( vector_p ); 

#define VECTOR_UNIT( vec )                   FUNCTION_CHECK( __my_vector_set_to_unit( &(vec) ), __my_vector_set_to_unit )

/* ------------------------------ */

int     __my_vector_array_allocate( int, int, vector_array_p );

#define VECTOR_ARRAY_ALLOCATE( n, dim_vec, az )   FUNCTION_CHECK(  __my_vector_array_allocate( (n), (dim_vec), &(az) ), __my_vector_array_allocate )

int     __my_vector_array_free( vector_array_p );

#define VECTOR_ARRAY_FREE( az )   FUNCTION_CHECK(  __my_vector_array_free( &(az) ),  __my_vector_array_free )

int     __my_vector_array_copy( vector_array_p, const vector_array );

#define VECTOR_ARRAY_COPY( az1, az2 )   FUNCTION_CHECK(  __my_vector_array_copy( &(az1), (az2) ),  __my_vector_array_copy )

int     __my_vector_array_compare( const vector_array, const vector_array );

#define VECTOR_ARRAY_COMPARE( az1, az2 )   FUNCTION_CHECK(  __my_vector_array_compare( (az1), (az2) ),  __my_vector_array_compare )

int     __my_vector_array_read( FILE*, vector_array_p );

#define VECTOR_ARRAY_READ( fp, az )   FUNCTION_CHECK(  __my_vector_array_read( fp, &(az) ),  __my_vector_array_read )

int     __my_vector_array_print( FILE*, const vector_array );

#define VECTOR_ARRAY_PRINT( fp, az )   FUNCTION_CHECK(  __my_vector_array_print( fp, (az) ),  __my_vector_array_print )

int     __my_vector_array_print_plus( FILE*, const vector_array );

#define VECTOR_ARRAY_PRINT_PLUS( fp, az )   FUNCTION_CHECK(  __my_vector_array_print_plus( fp, (az) ),  __my_vector_array_print_plus )

int     __my_vector_array_zero( vector_array_p );

#define VECTOR_ARRAY_ZERO( az )   FUNCTION_CHECK(  __my_vector_array_zero( &(az) ),  __my_vector_array_zero )

/* ------------------------------ */

/* int vectors */

typedef struct{

  int                ivector_dim;
  int*               ivector;

} __my_ivector;

typedef __my_ivector  ivector;
typedef __my_ivector* ivector_p;

int __my_ivector_allocate( const int, ivector_p );

#define IVECTOR_ALLOCATE( n, vec ) FUNCTION_CHECK( __my_ivector_allocate( (n), &(vec) ), __my_ivector_allocate )

int       __my_ivector_free( ivector_p );

#define IVECTOR_FREE( vec )     FUNCTION_CHECK( __my_ivector_free( &(vec) ),  __my_ivector_free )

int       __my_ivector_copy( ivector_p, const ivector );

#define IVECTOR_COPY( vec1, vec2 )   FUNCTION_CHECK(  __my_ivector_copy( &(vec1), (vec2) ), __my_ivector_copy )

int       __my_ivector_compare( const ivector, const ivector );

#define IVECTOR_COMPARE( vec1, vec2 )   FUNCTION_CHECK(  __my_ivector_compare( (vec1), (vec2) ), __my_ivector_compare )

int       __my_ivector_compare_mute( const ivector, const ivector );

#define IVECTOR_COMPARE_MUTE( vec1, vec2 )   __my_ivector_compare_mute( (vec1), (vec2) )

int       __my_ivector_read( FILE*, const ivector );

#define IVECTOR_READ( x, y )        FUNCTION_CHECK( __my_ivector_read( (x), (y) ), __my_ivector_read )

int       __my_ivector_print( FILE*, const ivector );

#define IVECTOR_PRINT( x, y )        FUNCTION_CHECK( __my_ivector_print( (x), (y) ), __my_ivector_print )

int       __my_ivector_print_plus( FILE*, const ivector );

#define IVECTOR_PRINT_PLUS( x, y )        FUNCTION_CHECK( __my_ivector_print_plus( (x), (y) ), __my_ivector_print_plus )

/* Errors handling */

int       __my_ivector_check( const ivector );

#define IVECTOR_CHECK( a )   FUNCTION_CHECK( __my_ivector_check( (a) ), __my_ivector_check )

int       __my_ivector_consinstency( const ivector, const ivector );

#define IVECTOR_CONSINSTENCY( a, b )  FUNCTION_CHECK( __my_ivector_consinstency( (a), (b) ), __my_ivector_consinstency )

int    __my_ivector_set_to_zero( ivector_p );

#define IVECTOR_ZERO( vec )                   FUNCTION_CHECK( __my_ivector_set_to_zero( &(vec) ), __my_ivector_set_to_zero )

int    __my_ivector_set_to_unit( ivector_p ); 

#define IVECTOR_UNIT( vec )                   FUNCTION_CHECK( __my_ivector_set_to_unit( &(vec) ), __my_ivector_set_to_unit )

int __my_ivector_are_equal( const ivector, const ivector );

#define IVECTOR_ARE_EQUAL( vec1, vec2 ) FUNCTION_CHECK( __my_ivector_are_equal( (vec1), (vec2) ), __my_ivector_are_equal )

/* ------------------------------ */

/* real vectors */

typedef struct{

  int                rvector_dim;
  double*            rvector;

} __my_rvector;

typedef __my_rvector  rvector;
typedef __my_rvector* rvector_p;

int __my_rvector_allocate( const int, rvector_p );

#define RVECTOR_ALLOCATE( n, vec ) FUNCTION_CHECK( __my_rvector_allocate( (n), &(vec) ), __my_rvector_allocate )

int       __my_rvector_free( rvector_p );

#define RVECTOR_FREE( vec )     FUNCTION_CHECK( __my_rvector_free( &(vec) ),  __my_rvector_free )

int       __my_rvector_copy( rvector_p, const rvector );

#define RVECTOR_COPY( vec1, vec2 )   FUNCTION_CHECK(  __my_rvector_copy( &(vec1), (vec2) ), __my_rvector_copy )

int       __my_rvector_compare( const rvector, const rvector );

#define RVECTOR_COMPARE( vec1, vec2 )   FUNCTION_CHECK(  __my_rvector_compare( (vec1), (vec2) ), __my_rvector_compare )

int       __my_rvector_read( FILE*, const rvector );

#define RVECTOR_READ( x, y )        FUNCTION_CHECK( __my_rvector_read( (x), (y) ), __my_rvector_read )

int       __my_rvector_print( FILE*, const rvector );

#define RVECTOR_PRINT( x, y )        FUNCTION_CHECK( __my_rvector_print( (x), (y) ), __my_rvector_print )

int       __my_rvector_print_plus( FILE*, const rvector );

#define RVECTOR_PRINT_PLUS( x, y )        FUNCTION_CHECK( __my_rvector_print_plus( (x), (y) ), __my_rvector_print_plus )

/* Errors handling */

int       __my_rvector_check( const rvector );

#define RVECTOR_CHECK( a )   FUNCTION_CHECK( __my_rvector_check( (a) ), __my_rvector_check )

int       __my_rvector_consinstency( const rvector, const rvector );

#define RVECTOR_CONSINSTENCY( a, b )  FUNCTION_CHECK( __my_rvector_consinstency( (a), (b) ), __my_rvector_consinstency )

int    __my_rvector_set_to_zero( rvector_p );

#define RVECTOR_ZERO( vec )                   FUNCTION_CHECK( __my_rvector_set_to_zero( &(vec) ), __my_rvector_set_to_zero )

int    __my_rvector_set_to_unit( rvector_p ); 

#define RVECTOR_UNIT( vec )                   FUNCTION_CHECK( __my_rvector_set_to_unit( &(vec) ), __my_rvector_set_to_unit )

#endif /* __PolyCEID_VECTOR__ */
