
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

#ifndef __PolyCEID_DIPOLE__
#define __PolyCEID_DIPOLE__


#include "../algebra/PolyCEID_linear_algebra.h"
#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../initial/PolyCEID_indexing.h"
#include "../hamiltonians/PolyCEID_hamiltonian_single.h"
#include "../hamiltonians/PolyCEID_centre_of_mass.h"
#include "../output/PolyCEID_energies.h"
#include "../initial/PolyCEID_initial_conditions.h"


/*********************
  FUNCTIONS & MACROS
*********************/


/* compute superpositions */

int compute_superpositions( const constants, state_p, config_p );

#define COMPUTE_SUPERPOSITIONS( co, st, cf ) FUNCTION_CHECK( compute_superpositions( (co), &(st), &(cf) ), compute_superpositions ) 

//------------------------------------------

/* compute classical_dipole */

int compute_classical_dipole( const constants, state_p, config_p );

#define COMPUTE_CLASSICAL_DIPOLE( co, st, cf ) FUNCTION_CHECK( compute_classical_dipole( (co), &(st), &(cf) ), compute_classical_dipole ) 

//------------------------------------------

#endif /* __PolyCEID_DIPOLE__ */
