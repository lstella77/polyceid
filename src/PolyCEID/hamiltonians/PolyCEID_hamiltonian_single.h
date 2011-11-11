
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

#ifndef __PolyCEID_HAMILTONIAN_SINGLE__
#define __PolyCEID_HAMILTONIAN_SINGLE__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "PolyCEID_distances.h"
#include "../initial/PolyCEID_indexing.h"
#include "PolyCEID_forces.h"
#include "PolyCEID_hamiltonian_single_CHAIN_MOD2.h"
#include "PolyCEID_hamiltonian_single_SSH.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* H_matrix_single update */

int H_matrix_single_update( const constants, state_p, config_p, matrix_p );

#define H_MATRIX_SINGLE_UPDATE( co, st, cf, mat ) FUNCTION_CHECK( H_matrix_single_update( (co), &(st), &(cf), &(mat) ), H_matrix_single_update ) 

//------------------------------------------

/* F_matrix_single update */

int F_matrix_single_update( const constants, state_p, config_p, int, matrix_p );

#define F_MATRIX_SINGLE_UPDATE( co, st, cf, i, mat ) FUNCTION_CHECK( F_matrix_single_update( (co), &(st), &(cf), (i), &(mat) ), F_matrix_single_update ) 

//------------------------------------------

/* K_matrix_single update */

int K_matrix_single_update( const constants, state_p, config_p, int, int, matrix_p );

#define K_MATRIX_SINGLE_UPDATE( co, st, cf, i, j, mat ) FUNCTION_CHECK( K_matrix_single_update( (co), &(st), &(cf), (i), (j), &(mat) ), K_matrix_single_update )

//------------------------------------------

/* utilities */

//------------------------------------------

/* initialise_hamiltonian_single */

int initialise_hamiltonian_single( constants );

#define INITIALISE_HAMILTONIAN_SINGLE( co )  FUNCTION_CHECK( initialise_hamiltonian_single( (co) ), initialise_hamiltonian_single ) 

//------------------------------------------
//------------------------------------------

#ifdef __DEBUG__

//------------------------------------------

int Hamiltonian_single_check( const constants, state_p, config_p );

#define HAMILTONIAN_SINGLE_CHECK( co, st, cf ) FUNCTION_CHECK( Hamiltonian_single_check( (co), &(st), &(cf) ), Hamiltonian_single_check )

//------------------------------------------
//------------------------------------------

int matrix_compare_entries( matrix_p, matrix_p, double );

#define MATRIX_COMPARE_ENTRIES( mat1, mat2, eps ) FUNCTION_CHECK( matrix_compare_entries( &(mat1), &(mat2), (eps) ), matrix_compare_entries )

//------------------------------------------

#endif /* __DEBUG__ */

/* H_matrix_single_check */

int H_matrix_single_check( const constants, rvector_p, matrix_p );

#define H_MATRIX_SINGLE_CHECK( co, vec, mat ) FUNCTION_CHECK( H_matrix_single_check( (co), &(vec), &(mat) ), H_matrix_single_check )

//------------------------------------------
//------------------------------------------

/* classical_dipole_single update */

int classical_dipole_single_update( const constants, state_p, config_p, int, matrix_p );

#define CLASSICAL_DIPOLE_SINGLE_UPDATE( co, st, cf, i, mat ) FUNCTION_CHECK( classical_dipole_single_update( (co), &(st), &(cf), (i), &(mat) ), classical_dipole_single_update ) 

//------------------------------------------

#endif /* __PolyCEID_HAMILTONIAN_SINGLE__ */
