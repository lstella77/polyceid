
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

#ifndef __PolyCEID_HAMILTONIAN_SINGLE_CHAIN_MOD2__
#define __PolyCEID_HAMILTONIAN_SINGLE_CHAIN_MOD2__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "PolyCEID_distances.h"
#include "../initial/PolyCEID_indexing.h"
#include "PolyCEID_forces.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* H_matrix_single_CHAIN_MOD2 update */

int H_matrix_single_CHAIN_MOD2_update( const constants, state_p, config_p, matrix_p );

#define H_MATRIX_SINGLE_CHAIN_MOD2_UPDATE( co, st, cf, mat ) FUNCTION_CHECK( H_matrix_single_CHAIN_MOD2_update( (co), &(st), &(cf), &(mat) ), H_matrix_single_CHAIN_MOD2_update ) 

//------------------------------------------

/* F_matrix_single_CHAIN_MOD2 update */

int F_matrix_single_CHAIN_MOD2_update( const constants, state_p, config_p, int, matrix_p );

#define F_MATRIX_SINGLE_CHAIN_MOD2_UPDATE( co, st, cf, i, mat ) FUNCTION_CHECK( F_matrix_single_CHAIN_MOD2_update( (co), &(st), &(cf), (i), &(mat) ), F_matrix_single_CHAIN_MOD2_update ) 

//------------------------------------------

/* K_matrix_single_CHAIN_MOD2 update */

int K_matrix_single_CHAIN_MOD2_update( const constants, state_p, config_p, int, int, matrix_p );

#define K_MATRIX_SINGLE_CHAIN_MOD2_UPDATE( co, st, cf, i, j, mat ) FUNCTION_CHECK( K_matrix_single_CHAIN_MOD2_update( (co), &(st), &(cf), (i), (j), &(mat) ), K_matrix_single_CHAIN_MOD2_update )

//------------------------------------------
//------------------------------------------

/* utilities */

int Hamiltonian_single_CHAIN_MOD2_parameters_check( const constants );

#define HAMILTONIAN_SINGLE_CHAIN_MOD2_PARAMETERS_CHECK( co ) FUNCTION_CHECK( Hamiltonian_single_CHAIN_MOD2_parameters_check( (co) ), Hamiltonian_single_CHAIN_MOD2_parameters_check )

//------------------------------------------
//------------------------------------------

int H_matrix_single_CHAIN_MOD2_check( const constants, rvector_p, matrix_p  );

#define H_MATRIX_SINGLE_CHAIN_MOD2_CHECK( co, pos, mat ) FUNCTION_CHECK( H_matrix_single_CHAIN_MOD2_check( co, &(pos), &(mat) ), H_matrix_single_CHAIN_MOD2_check )

//------------------------------------------

/* classical_dipole_single_CHAIN_MOD2_update */

int classical_dipole_single_CHAIN_MOD2_update( const constants, state_p, config_p, int, matrix_p );

#define CLASSICAL_DIPOLE_SINGLE_CHAIN_MOD2_UPDATE( co, st, cf, comp, mat ) FUNCTION_CHECK( classical_dipole_single_CHAIN_MOD2_update( (co), &(st), &(cf), (comp), &(mat) ), classical_dipole_single_CHAIN_MOD2_update )


//------------------------------------------

#endif /* __PolyCEID_HAMILTONIAN_SINGLE_CHAIN_MOD2__ */
