
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

#ifndef __PolyCEID_HAMILTONIAN_MANY__
#define __PolyCEID_HAMILTONIAN_MANY__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../hamiltonians/PolyCEID_distances.h"
#include "../hamiltonians/PolyCEID_centre_of_mass.h"
#include "../initial/PolyCEID_indexing.h"
#include "../hamiltonians/PolyCEID_forces.h"
#include "../hamiltonians/PolyCEID_hamiltonian_single.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* hamiltonian_many_update */

int hamiltonian_many_update( const constants, state_p, config_p );

#define HAMILTONIAN_MANY_UPDATE( co, st, cf ) FUNCTION_CHECK( hamiltonian_many_update( (co), &(st), &(cf) ), hamiltonian_many_update ) 

//------------------------------------------
//------------------------------------------

/* H_matrix_many update */

int H_matrix_many_update( const constants, state_p, config_p );

#define H_MATRIX_MANY_UPDATE( co, st, cf ) FUNCTION_CHECK( H_matrix_many_update( (co), &(st), &(cf) ), H_matrix_many_update )

//------------------------------------------

/* F_matrix_many update */

int F_matrix_many_update( const constants, state_p, config_p );

#define F_MATRIX_MANY_UPDATE( co, st, cf ) FUNCTION_CHECK( F_matrix_many_update( (co), &(st), &(cf) ), F_matrix_many_update ) 
//------------------------------------------

/* delta_F_matrix_many update */

int delta_F_matrix_many_update( const constants, state_p, config_p );

#define DELTA_F_MATRIX_MANY_UPDATE( co, st, cf ) FUNCTION_CHECK( delta_F_matrix_many_update( (co), &(st), &(cf) ), delta_F_matrix_many_update ) 

//------------------------------------------

/* K_matrix_many update */

int K_matrix_many_update( const constants, state_p, config_p );

#define K_MATRIX_MANY_UPDATE( co, st, cf ) FUNCTION_CHECK( K_matrix_many_update( (co), &(st), &(cf) ), K_matrix_many_update )

//------------------------------------------

/* classical_dipole_many update */

int classical_dipole_many_update( const constants, state_p, config_p, const int );

#define CLASSICAL_DIPOLE_MANY_UPDATE( co, st, cf, comp ) FUNCTION_CHECK( classical_dipole_many_update( (co), &(st), &(cf), (comp) ), classical_dipole_many_update)

//------------------------------------------

/* utilities */

//------------------------------------------

/* initialise_hamiltonian_many */

int initialise_hamiltonian_many( const constants );

#define INITIALISE_HAMILTONIAN_MANY( co )  FUNCTION_CHECK( initialise_hamiltonian_many( (co) ), initialise_hamiltonian_many ) 
//------------------------------------------

#ifdef __DEBUG__

//------------------------------------------

/* H_matrix_many_check */

int H_matrix_many_check( const constants, state_p, config_p, rvector_p, matrix_p );

#define H_MATRIX_MANY_CHECK( co, st, cf, vec, mat ) FUNCTION_CHECK( H_matrix_many_check( (co), &(st), &(cf), &(vec), &(mat) ), H_matrix_many_check )

//------------------------------------------
//------------------------------------------

int Hamiltonian_many_check( const constants, state_p, config_p );

#define HAMILTONIAN_MANY_CHECK( co, st, cf ) FUNCTION_CHECK( Hamiltonian_many_check( (co), &(st), &(cf) ), Hamiltonian_many_check )

//------------------------------------------

#endif /* __DEBUG__ */

//------------------------------------------

/* check_particle_hole_matrix */

int check_particle_hole_matrix( matrix_p, char*  );

#ifdef __CHECK_PARTICLE_HOLE__

#define CHECK_PARTICLE_HOLE_MATRIX( mat )   FUNCTION_CHECK( check_particle_hole_matrix( &(mat), (#mat) ), check_particle_hole_matrix )

#else /* __CHECK_PARTICLE_HOLE__ */

#define CHECK_PARTICLE_HOLE_MATRIX( mat )   0

#endif /* __CHECK_PARTICLE_HOLE__ */

//------------------------------------------

#endif /* __PolyCEID_HAMILTONIAN_MANY__ */
