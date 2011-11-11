
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

#ifndef __PolyCEID_INITIAL_RHO_ELECTRON__
#define __PolyCEID_INITIAL_RHO_ELECTRON__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../initial/PolyCEID_indexing.h"
#include "../hamiltonians/PolyCEID_hamiltonian_many.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* compute_initial_rho_electron */

int compute_initial_rho_electron( const constants, state_p );

#define COMPUTE_INITIAL_RHO_ELECTRON( co, st ) FUNCTION_CHECK( compute_initial_rho_electron( (co), &(st) ), compute_initial_rho_electron )

//------------------------------------------

/* compute_excited_rho_electron */

int compute_excited_rho_electron( const constants, state_p );

#define COMPUTE_EXCITED_RHO_ELECTRON( co, st ) FUNCTION_CHECK( compute_excited_rho_electron( (co), &(st) ), compute_excited_rho_electron )

//------------------------------------------
//------------------------------------------

int allocate_symmetry_multiplicity( constants_p );

#define ALLOCATE_SYMMETRY_MULTIPLICITY( co ) FUNCTION_CHECK( allocate_symmetry_multiplicity( &(co) ), allocate_symmetry_multiplicity )

//------------------------------------------

int allocate_many_body_hybridisation_table( constants_p );

#define ALLOCATE_MANY_BODY_HYBRIDISATION_TABLE( co ) FUNCTION_CHECK( allocate_many_body_hybridisation_table( &(co) ), allocate_many_body_hybridisation_table )

//------------------------------------------

int allocate_many_body_hybridisation_table_aux( const constants, const ivector, ivector, int*, int*, int*, int*, int* );

#define ALLOCATE_MANY_BODY_HYBRIDISATION_TABLE_AUX( c, ivec1, ivec2, orb1, spin1, orb2, spin2, sign ) FUNCTION_CHECK( \
                 allocate_many_body_hybridisation_table_aux( (c), (ivec1), (ivec2), &(orb1), &(spin1), &(orb2), &(spin2), &(sign) ), allocate_many_body_hybridisation_table_aux )


//------------------------------------------

/* utilities */

//------------------------------------------

int allocate_compute_initial_rho_electron_workspace( const constants );

#define ALLOCATE_COMPUTE_INITIAL_RHO_ELECTRON_WORKSPACE( co ) FUNCTION_CHECK( allocate_compute_initial_rho_electron_workspace( (co) ), allocate_compute_initial_rho_electron_workspace )

//------------------------------------------

int free_compute_initial_rho_electron_workspace( void );

#define FREE_COMPUTE_INITIAL_RHO_ELECTRON_WORKSPACE() FUNCTION_CHECK( free_compute_initial_rho_electron_workspace(), \
								      free_compute_initial_rho_electron_workspace )

//------------------------------------------

#endif /* __PolyCEID_INITIAL_RHO_ELECTRON__ */
