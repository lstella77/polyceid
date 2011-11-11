
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

#ifndef __PolyCEID_INDEXING__
#define __PolyCEID_INDEXING__


#include "../structures/PolyCEID_structures_constants.h"
#include <math.h>


/*********************
  PUBLIC VARIABLES
*********************/

imatrix  rho_index_aux_table;            // [sqrt_max_rho_index,sqrt_max_rho_index]
imatrix  coordinate_index_table;         // [N_coor,N_coor]
imatrix  reduced_coordinate_index_table; // [N_coor_red,N_coor_red]
imatrix  atom_index_table;               // [N_atoms,N_atoms]
imatrix  electron_many_index_table;      // [N_levels_many,N_levels_many]
imatrix  electron_single_index_table;    // [N_levelss_single,N_levels_single]

ivector  diagonal_rho_indices;     // [sqrt_max_rho_index]
ivector  rho_index_conjugate;      // [max_rho_index]
ivector  rho_index_border;         // [rho_index_border_length]
ivector  rho_index_next_to_border; // [rho_index_next_to_border_length]


ivector  dummy_rho_index1; // [N_coor_red]
ivector  dummy_rho_index2; // [N_coor_red]

/*********************
  FUNCTIONS & MACROS
*********************/

//------------------------------------------
//                  COORDINATES
//------------------------------------------

/* coordinate_index */

int coordinate_index( int i, int j );

#define COORDINATE_INDEX( i, j ) coordinate_index( (i), (j) )

//------------------------------------------

/* reduced_coordinate_index */

int reduced_coordinate_index( int i, int j );

#define REDUCED_COORDINATE_INDEX( i, j ) reduced_coordinate_index( (i), (j) )

//------------------------------------------
//                  ATOMS
//------------------------------------------

/* atom_index */

int atom_index( int i, int j );

#define ATOM_INDEX( i, j ) atom_index( (i), (j) )

//------------------------------------------
//                  ELECTRONS
//------------------------------------------

/* electron_many_index */

int electron_many_index( int i, int j );

#define ELECTRON_MANY_INDEX( i, j ) electron_many_index( (i), (j) )

//------------------------------------------

/* electron_single_index */

int electron_single_index( int i, int j );

#define ELECTRON_SINGLE_INDEX( i, j ) electron_single_index( (i), (j) )


//------------------------------------------
//                  RHO
//------------------------------------------

/* rho_index */

int rho_index( ivector, ivector );

#define RHO_INDEX( i, j ) rho_index( (i), (j) )

//------------------------------------------

/* rho_index_aux */

int rho_index_aux( int, int );

#define RHO_INDEX_AUX( i, j ) rho_index_aux( (i), (j) )

//------------------------------------------
//------------------------------------------

/* rho_index_change1_plus1_aux */

int rho_index_change1_plus1_aux( int index, int i );

#define RHO_INDEX_CHANGE1_PLUS1_AUX( index, i ) rho_index_change1_plus1_aux( (index), (i) )

//------------------------------------------

/* rho_index_change1_minus1_aux */

int rho_index_change1_minus1_aux( int index, int i );

#define RHO_INDEX_CHANGE1_MINUS1_AUX( index, i ) rho_index_change1_minus1_aux( (index), (i) )

//------------------------------------------

/* rho_index_change1_plus2_aux */

int rho_index_change1_plus2_aux( int index, int i );

#define RHO_INDEX_CHANGE1_PLUS2_AUX( index, i ) rho_index_change1_plus2_aux( (index), (i) )

//------------------------------------------

/* rho_index_change1_minus2_aux */

int rho_index_change1_minus2_aux( int index, int i );

#define RHO_INDEX_CHANGE1_MINUS2_AUX( index, i ) rho_index_change1_minus2_aux( (index), (i) )

//------------------------------------------

/* rho_index_change2_plus1_plus1_aux */

int rho_index_change2_plus1_plus1_aux( int index, int i, int j );

#define RHO_INDEX_CHANGE2_PLUS1_PLUS1_AUX( index, i, j ) rho_index_change2_plus1_plus1_aux( (index), (i), (j) )

//------------------------------------------

/* rho_index_change2_plus1_minus1_aux */

int rho_index_change2_plus1_minus1_aux( int index, int i, int j );

#define RHO_INDEX_CHANGE2_PLUS1_MINUS1_AUX( index, i, j ) rho_index_change2_plus1_minus1_aux( (index), (i), (j) )

//------------------------------------------

/* rho_index_change2_minus1_minus1_aux */

int rho_index_change2_minus1_minus1_aux( int index, int i, int j );

#define RHO_INDEX_CHANGE2_MINUS1_MINUS1_AUX( index, i, j ) rho_index_change2_minus1_minus1_aux( (index), (i), (j) )

//------------------------------------------
//------------------------------------------

/* rho_index_change1_plus1_first */

int rho_index_change1_plus1_first( int index, int i );

#define RHO_INDEX_CHANGE1_PLUS1_FIRST( index, i ) rho_index_change1_plus1_first( (index), (i) )

//------------------------------------------

/* rho_index_change1_minus1_first */

int rho_index_change1_minus1_first( int index, int i );

#define RHO_INDEX_CHANGE1_MINUS1_FIRST( index, i ) rho_index_change1_minus1_first( (index), (i) )

//------------------------------------------

/* rho_index_change1_plus2_first */

int rho_index_change1_plus2_first( int index, int i );

#define RHO_INDEX_CHANGE1_PLUS2_FIRST( index, i ) rho_index_change1_plus2_first( (index), (i) )

//------------------------------------------

/* rho_index_change1_minus2_first */

int rho_index_change1_minus2_first( int index, int i );

#define RHO_INDEX_CHANGE1_MINUS2_FIRST( index, i ) rho_index_change1_minus2_first( (index), (i) )

//------------------------------------------

/* rho_index_change2_plus1_plus1_first */

int rho_index_change2_plus1_plus1_first( int index, int i, int j );

#define RHO_INDEX_CHANGE2_PLUS1_PLUS1_FIRST( index, i, j ) rho_index_change2_plus1_plus1_first( (index), (i), (j) )

//------------------------------------------

/* rho_index_change2_plus1_minus1_first */

int rho_index_change2_plus1_minus1_first( int index, int i, int j );

#define RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( index, i, j ) rho_index_change2_plus1_minus1_first( (index), (i), (j) )

//------------------------------------------

/* rho_index_change2_minus1_minus1_first */

int rho_index_change2_minus1_minus1_first( int index, int i, int j );

#define RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index, i, j ) rho_index_change2_minus1_minus1_first( (index), (i), (j) )

//------------------------------------------
//------------------------------------------

/* rho_index_change1_plus1_second */

int rho_index_change1_plus1_second( int index, int i );

#define RHO_INDEX_CHANGE1_PLUS1_SECOND( index, i ) rho_index_change1_plus1_second( (index), (i) )

//------------------------------------------

/* rho_index_change1_minus1_second */

int rho_index_change1_minus1_second( int index, int i );

#define RHO_INDEX_CHANGE1_MINUS1_SECOND( index, i ) rho_index_change1_minus1_second( (index), (i) )

//------------------------------------------

/* rho_index_change1_plus2_second */

int rho_index_change1_plus2_second( int index, int i );

#define RHO_INDEX_CHANGE1_PLUS2_SECOND( index, i ) rho_index_change1_plus2_second( (index), (i) )

//------------------------------------------

/* rho_index_change1_minus2_second */

int rho_index_change1_minus2_second( int index, int i );

#define RHO_INDEX_CHANGE1_MINUS2_SECOND( index, i ) rho_index_change1_minus2_second( (index), (i) )

//------------------------------------------

/* rho_index_change2_plus1_plus1_second */

int rho_index_change2_plus1_plus1_second( int index, int i, int j );

#define RHO_INDEX_CHANGE2_PLUS1_PLUS1_SECOND( index, i, j ) rho_index_change2_plus1_plus1_second( (index), (i), (j) )

//------------------------------------------

/* rho_index_change2_plus1_minus1_second */

int rho_index_change2_plus1_minus1_second( int index, int i, int j );

#define RHO_INDEX_CHANGE2_PLUS1_MINUS1_SECOND( index, i, j ) rho_index_change2_plus1_minus1_second( (index), (i), (j) )

//------------------------------------------

/* rho_index_change2_minus1_minus1_second */

int rho_index_change2_minus1_minus1_second( int index, int i, int j );

#define RHO_INDEX_CHANGE2_MINUS1_MINUS1_SECOND( index, i, j ) rho_index_change2_minus1_minus1_second( (index), (i), (j) )

//------------------------------------------

/* utilities */

//------------------------------------------
//                  GENERAL
//------------------------------------------

/* initialise indexing  */

int indexing_allocate( constants, unsigned short int );

#define INDEXING_ALLOCATE( c, flag ) FUNCTION_CHECK( indexing_allocate( (c), (flag) ), indexing_allocate )

//------------------------------------------

/* indexing free  */

int indexing_free( unsigned short int  );

#define INDEXING_FREE( flag ) FUNCTION_CHECK( indexing_free( (flag) ), indexing_free )

//------------------------------------------
//                  COORDINATES
//------------------------------------------

/* coordinate_index_inverse */

int coordinate_index_inverse( int, int*, int* );

#define COORDINATE_INDEX_INVERSE( index, i, j ) FUNCTION_CHECK( coordinate_index_inverse( (index), &(i), &(j) ), coordinate_index_inverse )

//------------------------------------------

/* reduced_coordinate_index_inverse */

int reduced_coordinate_index_inverse( int, int*, int* );

#define REDUCED_COORDINATE_INDEX_INVERSE( index, i, j ) FUNCTION_CHECK( reduced_coordinate_index_inverse( (index), &(i), &(j) ), reduced_coordinate_index_inverse )

//------------------------------------------
//                  ATOMS
//------------------------------------------

/* atom_index_inverse */

int atom_index_inverse( int, int*, int* );

#define ATOM_INDEX_INVERSE( index, i, j ) FUNCTION_CHECK( atom_index_inverse( (index), &(i), &(j) ), atom_index_inverse )

//------------------------------------------
//                  ELECTRONS
//------------------------------------------

/* electron_many_index_inverse */

int electron_many_index_inverse( int, int*, int* );

#define ELECTRON_MANY_INDEX_INVERSE( index, i, j ) FUNCTION_CHECK( electron_many_index_inverse( (index), &(i), &(j) ), electron_many_index_inverse )

//------------------------------------------

/* electron_single_index_inverse */

int electron_single_index_inverse( int, int*, int* );

#define ELECTRON_SINGLE_INDEX_INVERSE( index, i, j ) FUNCTION_CHECK( electron_single_index_inverse( (index), &(i), &(j) ), electron_single_index_inverse )

//------------------------------------------
//                  RHO
//------------------------------------------

int compute_rho_dimensions( constants_p, int, int );

#define COMPUTE_RHO_DIMENSIONS( c, CEID_order, N_coor_red ) FUNCTION_CHECK( compute_rho_dimensions( &(c), (CEID_order), (N_coor_red) ), compute_rho_dimensions )

//------------------------------------------

/* rho_index_inverse */

int rho_index_inverse( int, ivector_p, ivector_p );

#define RHO_INDEX_INVERSE( index, i, j ) FUNCTION_CHECK( rho_index_inverse( (index), &(i), &(j) ), rho_index_inverse )

//------------------------------------------

/* rho_index_inverse_aux */

int rho_index_inverse_aux( int, ivector_p );

#define RHO_INDEX_INVERSE_AUX( index, i ) FUNCTION_CHECK( rho_index_inverse_aux( (index), &(i)  ), rho_index_inverse_aux )

//------------------------------------------

/* compute_rho_index_table */

int compute_rho_index_table( unsigned short int );

#define COMPUTE_RHO_INDEX_TABLE( flag ) FUNCTION_CHECK( compute_rho_index_table( (flag) ), compute_rho_index_table ) 

//------------------------------------------

/* search_rho_index_table */

int search_rho_index_table( ivector_p );

#define SEARCH_RHO_INDEX_TABLE( ind ) search_rho_index_table( &(ind) )

//------------------------------------------

/* compute_rho_index_other_tables */

int compute_rho_index_other_tables( void );

#define COMPUTE_RHO_INDEX_OTHER_TABLES() FUNCTION_CHECK( compute_rho_index_other_tables(), compute_rho_index_other_tables ) 

//------------------------------------------

/* check indexing */

int check_indexing( void );

#define CHECK_INDEXING() FUNCTION_CHECK( check_indexing(), check_indexing )

//------------------------------------------

#endif /* __PolyCEID_INDEXING__ */
