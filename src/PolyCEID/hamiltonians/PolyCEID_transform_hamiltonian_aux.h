
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

#ifndef __PolyCEID_TRANSFORM_HAMILTONIAN_AUX__
#define __PolyCEID_TRANSFORM_HAMILTONIAN_AUX__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../structures/PolyCEID_structures_config.h"
#include "../initial/PolyCEID_indexing.h"
#include "../algebra/PolyCEID_linear_algebra.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* transform_hamiltonian_aux */

int transform_hamiltonian_aux( constants, state_p, config_p, matrix_array_p, rvector_p );

#define TRANSFORM_HAMILTONIAN_AUX( co, st, cf, forces, momenta ) FUNCTION_CHECK( transform_hamiltonian_aux( (co), &(st), &(cf), &(forces), &(momenta) ), transform_hamiltonian_aux )

//------------------------------------------
//------------------------------------------

/* transform_F_matrix */

int transform_F_matrix( constants, state_p, config_p, matrix_array_p );

#define TRANSFORM_F_MATRIX( co, st, cf, forces ) FUNCTION_CHECK( transform_F_matrix( (co), &(st), &(cf), &(forces) ), transform_F_matrix )

//------------------------------------------

/* transform_momenta */

int transform_momenta( constants, state_p, config_p, rvector_p );

#define TRANSFORM_MOMENTA( co, st, cf, momenta ) FUNCTION_CHECK( transform_momenta( (co), &(st), &(cf), &(momenta) ), transform_momenta )

//------------------------------------------

/* transform_forces */

int transform_forces( constants, state_p, config_p, rvector_p );

#define TRANSFORM_FORCES( co, st, cf, forces ) FUNCTION_CHECK( transform_forces( (co), &(st), &(cf), &(forces) ), transform_forces )

//------------------------------------------

#endif /* __PolyCEID_TRANSFORM_HAMILTONIAN_AUX__ */
