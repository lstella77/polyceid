
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

#ifndef __PolyCEID_STRUCTURES_HAMILTONIAN__
#define __PolyCEID_STRUCTURES_HAMILTONIAN__


#include "../algebra/PolyCEID_constants.h"
#include "../algebra/PolyCEID_linear_algebra.h"
#include <string.h>
#include <stdio.h>


/* parameters */
#define NP_HAMILTONIAN         0  //bugfix: temporary!


/*********************
  HAMILTONIAN
*********************/

typedef struct{

  char      class[ MAX_STRING_LENGTH ];
  
  int       N_nodes;

  int       N_links;

  int*      neighbour_index;
     
  int*      neighbour_list;

  int*      neighbour_link;

  int       N_par_node;

  int       N_par_link;

  int       N_par_extra;

  double**  par_node;

  double**  par_link;

  double*   par_extra;

} PolyCEID_hamiltonian;

/*********************
  POINTERS
*********************/

typedef PolyCEID_hamiltonian  hamiltonian;

typedef PolyCEID_hamiltonian* hamiltonian_p;


/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_hamiltonian_allocate( int, int*, hamiltonian_p );

#define HAMILTONIAN_ALLOCATE( np, pp, h )  FUNCTION_CHECK( PolyCEID_hamiltonian_allocate( (np), (pp), &(h) ), PolyCEID_hamiltonian_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_hamiltonian_free( hamiltonian_p );

#define HAMILTONIAN_FREE( h )  FUNCTION_CHECK( PolyCEID_hamiltonian_free( &(h) ), PolyCEID_hamiltonian_free )

//------------------------------------------

/* READ */

int PolyCEID_hamiltonian_read( FILE*, hamiltonian_p );

#define HAMILTONIAN_READ( fp, h )  FUNCTION_CHECK( PolyCEID_hamiltonian_read( fp, &(h) ), PolyCEID_hamiltonian_read )

//------------------------------------------

/* PRINT */

int PolyCEID_hamiltonian_print( FILE*, const hamiltonian );

#define HAMILTONIAN_PRINT( fp, h )  FUNCTION_CHECK( PolyCEID_hamiltonian_print( (fp), (h) ), PolyCEID_hamiltonian_print )

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_hamiltonian_verbose_print( FILE*, const hamiltonian );

#define HAMILTONIAN_VERBOSE_PRINT( fp, h )  FUNCTION_CHECK( PolyCEID_hamiltonian_verbose_print( (fp), (h) ), PolyCEID_hamiltonian_verbose_print )

//------------------------------------------


/* COPY [ h1 = h2 ] */

int PolyCEID_hamiltonian_copy( hamiltonian_p, const hamiltonian );

#define HAMILTONIAN_COPY( h1, h2 )  FUNCTION_CHECK( PolyCEID_hamiltonian_copy( &(h1), (h2) ), PolyCEID_hamiltonian_copy )

//------------------------------------------

/* COMPARE [ h1 = h2 ] */

int PolyCEID_hamiltonian_compare( const hamiltonian, const hamiltonian );

#define HAMILTONIAN_COMPARE( h1, h2 )  FUNCTION_CHECK( PolyCEID_hamiltonian_compare( (h1), (h2) ), PolyCEID_hamiltonian_compare )

//------------------------------------------

#endif /* __PolyCEID_STRUCTURES_HAMILTONIAN__ */

