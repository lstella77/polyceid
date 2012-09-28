
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

#ifndef __PolyCEID_NONADIABATIC_COUPLING__
#define __PolyCEID_NONADIABATIC_COUPLING__


#include "../hamiltonians/PolyCEID_hamiltonian_many.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian_aux.h"
#include "../algebra/PolyCEID_combinatorics.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int compute_nonadiabaticity( constants, state_p, config_p );

#define COMPUTE_NONADIABATICITY( co, st, cf ) FUNCTION_CHECK( compute_nonadiabaticity( (co), &(st), &(cf) ), compute_nonadiabaticity )

//------------------------------------------

/* utilities */

//------------------------------------------

int compute_adiabatic_populations( const constants, state_p, config_p, matrix_p, rvector_p );

#define COMPUTE_ADIABATIC_POPULATIONS( co, st, cf, astates, populations ) FUNCTION_CHECK(  compute_adiabatic_populations( (co), &(st), &(cf), &(astates), &(populations) ),  compute_adiabatic_populations )

//------------------------------------------

int compute_nonadiabatic_forces( const constants, state_p, config_p, matrix_p, int, matrix_p );

#define COMPUTE_NONADIABATIC_FORCES( co, st, cf, astates, mode, forces ) FUNCTION_CHECK(  compute_nonadiabatic_forces( (co), &(st), &(cf), &(astates), (mode), &(forces) ),  compute_nonadiabatic_forces )

//------------------------------------------

#endif /* __PolyCEID_NONADIABATIC_COUPLING__ */
