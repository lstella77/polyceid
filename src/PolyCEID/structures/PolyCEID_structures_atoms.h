
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

#ifndef __PolyCEID_STRUCTURES_ATOMS__
#define __PolyCEID_STRUCTURES_ATOMS__


#include "../algebra/PolyCEID_constants.h"
#include "../algebra/PolyCEID_linear_algebra.h"
#include <string.h>
#include <stdio.h>


/* parameters */
#define NP_ATOMS  2

#define ATOM_NAME_LENGTH  8


/*********************
  ATOMS
*********************/

typedef struct{

  int spacial_dimension;
  int N_atoms;

  char**   names;           // [ATOM_NAME_LEGTH][N_coor]           vectorised!
  rvector  masses;          // [N_coor]           vectorised!
  rvector  masses_aux;      // [N_coor]
  double   mass_tot;        // total mass
  double   mass_ave;        // averaged mass
  rvector  positions;       // [N_coor]           vectorised!
  rvector  momenta;         // [N_coor]           vectorised!
  rvector  forces;          // [N_coor]           vectorised!
  rvector  forces_cons;     // [N_coor]           vectorised!
  rvector  distances;       // [N_atoms *N_atoms] vectorised, as a matrix.
  rvector  centre_of_mass;  // [sdim]
  ivector  mask;            // [N_atom]           vectorised!

} PolyCEID_atoms;

/*********************
  POINTERS
*********************/

typedef PolyCEID_atoms  atoms;

typedef PolyCEID_atoms* atoms_p;


/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_atoms_allocate( int, int*, atoms_p );

#define ATOMS_ALLOCATE( np, pp, a )  FUNCTION_CHECK( PolyCEID_atoms_allocate( (np), (pp), &(a) ), PolyCEID_atoms_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_atoms_free( atoms_p );

#define ATOMS_FREE( a )  FUNCTION_CHECK( PolyCEID_atoms_free( &(a) ), PolyCEID_atoms_free )

//------------------------------------------

/* READ */

int PolyCEID_atoms_read( FILE*, atoms_p );

#define ATOMS_READ( fp, a )  FUNCTION_CHECK( PolyCEID_atoms_read( fp, &(a) ), PolyCEID_atoms_read )

//------------------------------------------

/* PRINT */

int PolyCEID_atoms_print( FILE*, const atoms );

#define ATOMS_PRINT( fp, a )  FUNCTION_CHECK( PolyCEID_atoms_print( (fp), (a) ), PolyCEID_atoms_print )

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_atoms_verbose_print( FILE*, const atoms );

#define ATOMS_VERBOSE_PRINT( fp, a )  FUNCTION_CHECK( PolyCEID_atoms_verbose_print( (fp), (a) ), PolyCEID_atoms_verbose_print )

//------------------------------------------


/* COPY [ a1 = a2 ] */

int PolyCEID_atoms_copy( atoms_p, const atoms );

#define ATOMS_COPY( a1, a2 )  FUNCTION_CHECK( PolyCEID_atoms_copy( &(a1), (a2) ), PolyCEID_atoms_copy )

//------------------------------------------

/* COMPARE [ a1 = a2 ] */

int PolyCEID_atoms_compare( const atoms, const atoms );

#define ATOMS_COMPARE( a1, a2 )  FUNCTION_CHECK( PolyCEID_atoms_compare( (a1), (a2) ), PolyCEID_atoms_compare )

//------------------------------------------

#endif /* __PolyCEID_STRUCTURES_ATOMS__ */

