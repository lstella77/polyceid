
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

#ifndef __PolyCEID_HAMILTONIAN_SINGLE_SSH__
#define __PolyCEID_HAMILTONIAN_SINGLE_SSH__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "PolyCEID_distances.h"
#include "../initial/PolyCEID_indexing.h"
#include "PolyCEID_forces.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* H_matrix_single_SSH update */

int H_matrix_single_SSH_update( const constants, state_p, config_p, matrix_p );

#define H_MATRIX_SINGLE_SSH_UPDATE( co, st, cf, mat ) FUNCTION_CHECK( H_matrix_single_SSH_update( (co), &(st), &(cf), &(mat) ), H_matrix_single_SSH_update ) 

//------------------------------------------

/* F_matrix_single_SSH update */

int F_matrix_single_SSH_update( const constants, state_p, config_p, int, matrix_p );

#define F_MATRIX_SINGLE_SSH_UPDATE( co, st, cf, i, mat ) FUNCTION_CHECK( F_matrix_single_SSH_update( (co), &(st), &(cf), (i), &(mat) ), F_matrix_single_SSH_update ) 

//------------------------------------------

/* K_matrix_single_SSH update */

int K_matrix_single_SSH_update( const constants, state_p, config_p, int, int, matrix_p );

#define K_MATRIX_SINGLE_SSH_UPDATE( co, st, cf, i, j, mat ) FUNCTION_CHECK( K_matrix_single_SSH_update( (co), &(st), &(cf), (i), (j), &(mat) ), K_matrix_single_SSH_update )

//------------------------------------------
//------------------------------------------

/* utilities */

int Hamiltonian_single_SSH_parameters_check( const constants );

#define HAMILTONIAN_SINGLE_SSH_PARAMETERS_CHECK( c ) FUNCTION_CHECK( Hamiltonian_single_SSH_parameters_check( (c) ), Hamiltonian_single_SSH_parameters_check )

//------------------------------------------
//------------------------------------------

int H_matrix_single_SSH_check( const constants, rvector_p, matrix_p  );

#define H_MATRIX_SINGLE_SSH_CHECK( co, pos, mat ) FUNCTION_CHECK( H_matrix_single_SSH_check( co, &(pos), &(mat) ), H_matrix_single_SSH_check )

//------------------------------------------
//------------------------------------------

double distance_loc( const constants, int, int, rvector_p );

#define DISTANCE_LOC( co, a1, a2, pos ) distance_loc( (co), a1, a2, &(pos) )

//------------------------------------------

/* classical_dipole_single_SSH_update */

int classical_dipole_single_SSH_update( const constants, state_p, config_p, int, matrix_p );

#define CLASSICAL_DIPOLE_SINGLE_SSH_UPDATE( co, st, cf, i, mat ) FUNCTION_CHECK( classical_dipole_single_SSH_update( (co), &(st), &(cf), (i), &(mat) ), classical_dipole_single_SSH_update )

//------------------------------------------

#endif /* __PolyCEID_HAMILTONIAN_SINGLE_SSH__ */

