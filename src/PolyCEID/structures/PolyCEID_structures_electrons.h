
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

#ifndef __PolyCEID_STRUCTURES_ELECTRONS__
#define __PolyCEID_STRUCTURES_ELECTRONS__


#include "../algebra/PolyCEID_constants.h"
#include "../algebra/PolyCEID_linear_algebra.h"
#include <stdio.h>


/* parameters */
#define NP_ELECTRONS  4


/*********************
  ELECTRONS
*********************/

typedef struct{

  matrix        H_matrix;         //
  matrix_array  F_matrix;         // [N_coor]
  matrix_array  delta_F_matrix;   // [N_coor]
  matrix_array  K_matrix;         // [N_coor *N_coor] vectorised!
  //
  matrix        dipole_many;      //
  //
  matrix        Ehrenfest_frame;  //[N_levels_single,N_levels_single]  
  //
  matrix_array  rho;              // [max_rho_index+1] vectorised!
  //
  matrix        mu00;             // reduced electronic density matrix
  matrix_array  mu01;             // [N_coor], P type
  matrix_array  mu10;             // [N_coor], R type
  matrix_array  mu02;             // [N_coor *N_coor] vectorised!, PP type
  matrix_array  mu11;             // [N_coor *N_coor] vectorised!, PR type (symmetrised)
  matrix_array  mu20;             // [N_coor *N_coor] vectorised!, RR type

} PolyCEID_electrons;


/*********************
  POINTERS
*********************/

typedef PolyCEID_electrons  electrons;

typedef PolyCEID_electrons* electrons_p;


/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_electrons_allocate( int, int*, electrons_p );

#define ELECTRONS_ALLOCATE( np, pp, e )  FUNCTION_CHECK( PolyCEID_electrons_allocate( (np), (pp), &(e) ), PolyCEID_electrons_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_electrons_free( electrons_p );

#define ELECTRONS_FREE( e )  FUNCTION_CHECK( PolyCEID_electrons_free( &(e) ), PolyCEID_electrons_free )

//------------------------------------------

/* READ */

int PolyCEID_electrons_read( FILE*, electrons_p );

#define ELECTRONS_READ( fp, e )  FUNCTION_CHECK( PolyCEID_electrons_read( fp, &(e) ), PolyCEID_electrons_read )

//------------------------------------------

/* PRINT */

int PolyCEID_electrons_print( FILE*, const electrons );

#define ELECTRONS_PRINT( fp, e )  FUNCTION_CHECK( PolyCEID_electrons_print( (fp), (e) ), PolyCEID_electrons_print )

//------------------------------------------

/* PRINT VERBOSE */

int PolyCEID_electrons_verbose_print( FILE*, const electrons );

#define ELECTRONS_VERBOSE_PRINT( fp, e )  FUNCTION_CHECK( PolyCEID_electrons_verbose_print( (fp), (e) ), PolyCEID_electrons_verbose_print )

//------------------------------------------

/* COPY [ e1 = e2 ] */

int PolyCEID_electrons_copy( electrons_p, const electrons );

#define ELECTRONS_COPY( e1, e2 )  FUNCTION_CHECK( PolyCEID_electrons_copy( &(e1), (e2) ), PolyCEID_electrons_copy )

//------------------------------------------

/* COMPARE [ e1 = e2 ] */

int PolyCEID_electrons_compare( const electrons, const electrons );

#define ELECTRONS_COMPARE( e1, e2 )  FUNCTION_CHECK( PolyCEID_electrons_compare( (e1), (e2) ), PolyCEID_electrons_compare )

//------------------------------------------

#endif /* __PolyCEID_STRUCTURES_ELECTRONS__ */

