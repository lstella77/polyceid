
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

#ifndef __PolyCEID_INITIAL_CONDITIONS__
#define __PolyCEID_INITIAL_CONDITIONS__


#include "../initial/PolyCEID_normal_modes.h"
#include "../hamiltonians/PolyCEID_effective_hamiltonian.h"
#include "../evolution/PolyCEID_transform_mu.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian.h"
#include "../output/PolyCEID_one_body_electronic_density_matrix.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int initial_condition_rho( const constants, state_p, config_p );

#define INITIAL_CONDITION_RHO( co, st, cf ) FUNCTION_CHECK( initial_condition_rho( (co), &(st), &(cf) ), initial_condition_rho )

//------------------------------------------

/* utilities */

//------------------------------------------

/*  compute_initial_condition_atoms */

int compute_initial_condition_atoms( const constants, state_p, config_p );

#define COMPUTE_INITIAL_CONDITION_ATOMS( co, st, cf ) FUNCTION_CHECK( compute_initial_condition_atoms( (co), &(st), &(cf) ), compute_initial_condition_atoms )

//------------------------------------------

#endif /*  __PolyCEID_INITIAL_CONDITIONS__ */
