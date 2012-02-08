
/******************************************************************************

  Copyright (C) 2011-2012 by Lorenzo Stella <lorenzo DOT stella77 AT gmail.com>

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

#ifndef __PolyCEID_ENRGY_CORRECTIONS__
#define __PolyCEID_ENRGY_CORRECTIONS__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../initial/PolyCEID_indexing.h"
#include "../initial/PolyCEID_rho_extra_indexing.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian_aux.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int compute_kinetic_energy_corrections( const constants, state_p, config_p );

#define COMPUTE_KINETIC_ENERGY_CORRECTIONS( co, st, cf ) FUNCTION_CHECK( compute_kinetic_energy_corrections( (co), &(st), &(cf) ), compute_kinetic_energy_corrections )

//------------------------------------------

int compute_potential_energy_corrections( const constants, state_p, config_p );

#define COMPUTE_POTENTIAL_ENERGY_CORRECTIONS( co, st, cf ) FUNCTION_CHECK( compute_potential_energy_corrections( (co), &(st), &(cf) ), compute_potential_energy_corrections )

//------------------------------------------

#endif /* __PolyCEID_ENREGY_CORRECTIONS__ */
