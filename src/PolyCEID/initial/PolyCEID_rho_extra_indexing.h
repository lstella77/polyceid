
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

#ifndef __PolyCEID_RHO_EXTRA_INDEXING__
#define __PolyCEID_RHO_EXTRA_INDEXING__


#include "../structures/PolyCEID_structures_constants.h"
#include "../initial/PolyCEID_indexing.h"
#include <math.h>


/*********************
  PUBLIC VARIABLES
*********************/

imatrix  rho_extra_index_aux_table; // [sqrt_max_rho_index,sqrt_max_rho_index]

ivector  dummy_rho_extra_index1;    // [N_coor_red]
ivector  dummy_rho_extra_index2;    // [N_coor_red]

ivector  rho_index_aux_translate;   // [sqrt_max_rho_index]


/*********************
  FUNCTIONS & MACROS
*********************/

//------------------------------------------
//                  RHO_EXTRA
//------------------------------------------

/* rho_extra_index */

int rho_extra_index( ivector, ivector );

#define RHO_EXTRA_INDEX( i, j ) rho_extra_index( (i), (j) )

//------------------------------------------

/* rho_extra_index_aux */

int rho_extra_index_aux( int, int );

#define RHO_EXTRA_INDEX_AUX( i, j ) rho_extra_index_aux( (i), (j) )

//------------------------------------------
//------------------------------------------

/* rho_extra_index_change1_plus1_aux */

int rho_extra_index_change1_plus1_aux( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( index, i ) rho_extra_index_change1_plus1_aux( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_minus1_aux */

int rho_extra_index_change1_minus1_aux( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_MINUS1_AUX( index, i ) rho_extra_index_change1_minus1_aux( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_plus2_aux */

int rho_extra_index_change1_plus2_aux( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( index, i ) rho_extra_index_change1_plus2_aux( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_minus2_aux */

int rho_extra_index_change1_minus2_aux( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_MINUS2_AUX( index, i ) rho_extra_index_change1_minus2_aux( (index), (i) )

//------------------------------------------

/* rho_extra_index_change2_plus1_plus1_aux */

int rho_extra_index_change2_plus1_plus1_aux( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( index, i, j ) rho_extra_index_change2_plus1_plus1_aux( (index), (i), (j) )

//------------------------------------------

/* rho_extra_index_change2_plus1_minus1_aux */

int rho_extra_index_change2_plus1_minus1_aux( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_AUX( index, i, j ) rho_extra_index_change2_plus1_minus1_aux( (index), (i), (j) )

//------------------------------------------

/* rho_extra_index_change2_minus1_minus1_aux */

int rho_extra_index_change2_minus1_minus1_aux( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_AUX( index, i, j ) rho_extra_index_change2_minus1_minus1_aux( (index), (i), (j) )

//------------------------------------------
//------------------------------------------

/* rho_extra_index_change1_plus1_first */

int rho_extra_index_change1_plus1_first( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_PLUS1_FIRST( index, i ) rho_extra_index_change1_plus1_first( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_minus1_first */

int rho_extra_index_change1_minus1_first( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_MINUS1_FIRST( index, i ) rho_extra_index_change1_minus1_first( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_plus2_first */

int rho_extra_index_change1_plus2_first( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_PLUS2_FIRST( index, i ) rho_extra_index_change1_plus2_first( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_minus2_first */

int rho_extra_index_change1_minus2_first( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_MINUS2_FIRST( index, i ) rho_extra_index_change1_minus2_first( (index), (i) )

//------------------------------------------

/* rho_extra_index_change2_plus1_plus1_first */

int rho_extra_index_change2_plus1_plus1_first( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_FIRST( index, i, j ) rho_extra_index_change2_plus1_plus1_first( (index), (i), (j) )

//------------------------------------------

/* rho_extra_index_change2_plus1_minus1_first */

int rho_extra_index_change2_plus1_minus1_first( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( index, i, j ) rho_extra_index_change2_plus1_minus1_first( (index), (i), (j) )

//------------------------------------------

/* rho_extra_index_change2_minus1_minus1_first */

int rho_extra_index_change2_minus1_minus1_first( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index, i, j ) rho_extra_index_change2_minus1_minus1_first( (index), (i), (j) )

//------------------------------------------
//------------------------------------------

/* rho_extra_index_change1_plus1_second */

int rho_extra_index_change1_plus1_second( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_PLUS1_SECOND( index, i ) rho_extra_index_change1_plus1_second( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_minus1_second */

int rho_extra_index_change1_minus1_second( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_MINUS1_SECOND( index, i ) rho_extra_index_change1_minus1_second( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_plus2_second */

int rho_extra_index_change1_plus2_second( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_PLUS2_SECOND( index, i ) rho_extra_index_change1_plus2_second( (index), (i) )

//------------------------------------------

/* rho_extra_index_change1_minus2_second */

int rho_extra_index_change1_minus2_second( int index, int i );

#define RHO_EXTRA_INDEX_CHANGE1_MINUS2_SECOND( index, i ) rho_extra_index_change1_minus2_second( (index), (i) )

//------------------------------------------

/* rho_extra_index_change2_plus1_plus1_second */

int rho_extra_index_change2_plus1_plus1_second( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_SECOND( index, i, j ) rho_extra_index_change2_plus1_plus1_second( (index), (i), (j) )

//------------------------------------------

/* rho_extra_index_change2_plus1_minus1_second */

int rho_extra_index_change2_plus1_minus1_second( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_SECOND( index, i, j ) rho_extra_index_change2_plus1_minus1_second( (index), (i), (j) )

//------------------------------------------

/* rho_extra_index_change2_minus1_minus1_second */

int rho_extra_index_change2_minus1_minus1_second( int index, int i, int j );

#define RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_SECOND( index, i, j ) rho_extra_index_change2_minus1_minus1_second( (index), (i), (j) )

//------------------------------------------

/* utilities */

//------------------------------------------
//                  GENERAL
//------------------------------------------

/* initialise rho_extra_indexing  */

int rho_extra_indexing_allocate( constants );

#define RHO_EXTRA_INDEXING_ALLOCATE( c ) FUNCTION_CHECK( rho_extra_indexing_allocate( (c) ), rho_extra_indexing_allocate )

//------------------------------------------

/* rho_extra_indexing free  */

int rho_extra_indexing_free( void );

#define RHO_EXTRA_INDEXING_FREE() FUNCTION_CHECK( rho_extra_indexing_free(), rho_extra_indexing_free )

//------------------------------------------
//                  RHO_EXTRA
//------------------------------------------

/* rho_extra_index_inverse */

int rho_extra_index_inverse( int, ivector_p, ivector_p );

#define RHO_EXTRA_INDEX_INVERSE( index, i, j ) FUNCTION_CHECK( rho_extra_index_inverse( (index), &(i), &(j) ), rho_extra_index_inverse )

//------------------------------------------

/* rho_extra_index_inverse_aux */

int rho_extra_index_inverse_aux( int, ivector_p );

#define RHO_EXTRA_INDEX_INVERSE_AUX( index, i ) FUNCTION_CHECK( rho_extra_index_inverse_aux( (index), &(i)  ), rho_extra_index_inverse_aux )

//------------------------------------------

/* compute_rho_extra_index_table */

int compute_rho_extra_index_table( unsigned short int  );

#define COMPUTE_RHO_EXTRA_INDEX_TABLE( flag ) FUNCTION_CHECK( compute_rho_extra_index_table( (flag) ), compute_rho_extra_index_table ) 

//------------------------------------------

/* search_rho_extra_index_table */

int search_rho_extra_index_table( ivector_p );

#define SEARCH_RHO_EXTRA_INDEX_TABLE( ind ) search_rho_extra_index_table( &(ind) )

//------------------------------------------

/* compute_rho_extra_index_other_tables */

int compute_rho_extra_index_other_tables( void );

#define COMPUTE_RHO_EXTRA_INDEX_OTHER_TABLES() FUNCTION_CHECK( compute_rho_extra_index_other_tables(), compute_rho_extra_index_other_tables ) 

//------------------------------------------
//------------------------------------------

/* rho_index_aux_to_rho_extra_index_aux */

int rho_index_aux_to_rho_extra_index_aux( int );

#define RHO_INDEX_AUX_TO_RHO_EXTRA_INDEX_AUX( i ) rho_index_aux_to_rho_extra_index_aux( (i) )

//------------------------------------------

/* compute_rho_index_aux_translate*/

int compute_rho_index_aux_translate( const constants );

#define COMPUTE_RHO_INDEX_AUX_TRANSLATE( c ) FUNCTION_CHECK( compute_rho_index_aux_translate( (c) ), compute_rho_index_aux_translate )

//------------------------------------------
//------------------------------------------

/* check rho_extra_indexing */

int check_rho_extra_indexing( void );

#define CHECK_RHO_EXTRA_INDEXING() FUNCTION_CHECK( check_rho_extra_indexing(), check_rho_extra_indexing )

//------------------------------------------

#endif /* __PolyCEID_RHO_EXTRA_INDEXING__ */
