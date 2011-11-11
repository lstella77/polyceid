
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

#ifndef __PolyCEID_EFFECTIVE_HAMILTONIAN__
#define __PolyCEID_EFFECTIVE_HAMILTONIAN__


#include "../evolution/PolyCEID_transform_mu.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian_aux.h"
#include "../initial/PolyCEID_initial_rho_electron.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* effective_hamiltonian_update */

int effective_hamiltonian_update( const constants, state_p, config_p, matrix_p );

#define EFFECTIVE_HAMILTONIAN_UPDATE( co, st, cf, mat ) FUNCTION_CHECK( effective_hamiltonian_update( (co), &(st), &(cf), &(mat) ), effective_hamiltonian_update )

//------------------------------------------

/* effective_hamiltonian_update_aux */

int effective_hamiltonian_update_aux( const constants, state_p, config_p, vector, matrix_p );

#define EFFECTIVE_HAMILTONIAN_UPDATE_AUX( co, st, cf, vec, mat ) FUNCTION_CHECK( effective_hamiltonian_update_aux( (co), &(st), &(cf), (vec), &(mat) ), effective_hamiltonian_update_aux )

//------------------------------------------

/* utilities */

//------------------------------------------

int allocate_effective_hamiltonian_update_workspace( const constants );

#define ALLOCATE_EFFECTIVE_HAMILTONIAN_UPDATE_WORKSPACE( co ) FUNCTION_CHECK( allocate_effective_hamiltonian_update_workspace( (co) ), allocate_effective_hamiltonian_update_workspace )

//------------------------------------------

int free_effective_hamiltonian_update_workspace( void );

#define FREE_EFFECTIVE_HAMILTONIAN_UPDATE_WORKSPACE() FUNCTION_CHECK( free_effective_hamiltonian_update_workspace(), \
								      free_effective_hamiltonian_update_workspace )

//------------------------------------------

#endif /* __PolyCEID_EFFECTIVE_HAMILTONIAN__ */
