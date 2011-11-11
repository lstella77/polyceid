
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

#ifndef __PolyCEID_STRUCTURES_PHONONS__
#define __PolyCEID_STRUCTURES_PHONONS__


#include "../algebra/PolyCEID_constants.h"
#include "../input/PolyCEID_parsing.h"
#include "../algebra/PolyCEID_linear_algebra.h"


#define NP_PHONONS  2


/*********************
  PHONONS
*********************/

typedef struct{

  /* phonons */
  vector              initial_condition_atoms;           // [sqrt_max_rho_index] 
  matrix              original_hessian;                  // [N_coor,N_coor]
  matrix              canonical_transform;               // [N_coor,N_coor]
  rvector             canonical_frequencies;             // [N_coor]
  rvector             omega;                             // [N_coor] Are the last two different?
  rvector             ar;                                // [N_coor]
  rvector             ap;                                // [N_coor]
  rvector             phononic_dos;                      // [N_coor]

  matrix              symplectic_transform_R;            //  [N_coor,N_coor]
  matrix              conjugate_symplectic_transform_R;  //  [N_coor,N_coor]
  matrix              symplectic_transform_P;            //  [N_coor,N_coor]
  matrix              conjugate_symplectic_transform_P;  //  [N_coor,N_coor]

} PolyCEID_phonons;


/*********************
  POINTERS
*********************/

typedef PolyCEID_phonons phonons;

typedef PolyCEID_phonons* phonons_p;

/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_phonons_allocate( int, int*, phonons_p );

#define PHONONS_ALLOCATE( np, pp, c )  FUNCTION_CHECK( PolyCEID_phonons_allocate( (np), (pp), &(c) ), PolyCEID_phonons_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_phonons_free( phonons_p );

#define PHONONS_FREE( c )  FUNCTION_CHECK( PolyCEID_phonons_free( &(c) ), PolyCEID_phonons_free )

//------------------------------------------

/* READ */

int PolyCEID_phonons_read( FILE*, phonons_p );

#define PHONONS_READ( fp, c )  FUNCTION_CHECK( PolyCEID_phonons_read( fp, &(c) ), PolyCEID_phonons_read )

//------------------------------------------

/* PRINT */

int PolyCEID_phonons_print( FILE*, const phonons );

#define PHONONS_PRINT( fp, c )  FUNCTION_CHECK( PolyCEID_phonons_print( (fp), (c) ), PolyCEID_phonons_print )

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_phonons_verbose_print( FILE*, const phonons );

#define PHONONS_VERBOSE_PRINT( fp, c )  FUNCTION_CHECK( PolyCEID_phonons_verbose_print( (fp), (c) ), PolyCEID_phonons_verbose_print )

//------------------------------------------


/* COPY [ c1 = c2 ] */

int PolyCEID_phonons_copy( phonons_p, const phonons );

#define PHONONS_COPY( c1, c2 )  FUNCTION_CHECK( PolyCEID_phonons_copy( &(c1), (c2) ), PolyCEID_phonons_copy )

//------------------------------------------

#endif /* __PolyCEID_STRUCTURES_PHONONS__ */
