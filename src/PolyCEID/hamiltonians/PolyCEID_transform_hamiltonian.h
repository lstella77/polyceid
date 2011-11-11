
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

#ifndef __PolyCEID_TRANSFORM_HAMILTONIAN__
#define __PolyCEID_TRANSFORM_HAMILTONIAN__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../structures/PolyCEID_structures_config.h"
#include "../initial/PolyCEID_indexing.h"
#include "../algebra/PolyCEID_linear_algebra.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* transform_hamiltonian */

int transform_hamiltonian( constants, state_p, config_p, matrix_array_p, matrix_array_p );

#define TRANSFORM_HAMILTONIAN( co, st, cf, fors, pots ) FUNCTION_CHECK( transform_hamiltonian( (co), &(st), &(cf), &(fors), &(pots) ), transform_hamiltonian )

//------------------------------------------
//------------------------------------------

/* transform_delta_F_matrix */

int transform_delta_F_matrix( constants, state_p, config_p, matrix_array_p );

#define TRANSFORM_DELTA_F_MATRIX( co, st, cf, fors ) FUNCTION_CHECK( transform_delta_F_matrix( (co), &(st), &(cf), &(fors) ), transform_delta_F_matrix )

//------------------------------------------

/* transform_K_matrix */

int transform_K_matrix( constants, state_p, config_p, matrix_array_p );

#define TRANSFORM_K_MATRIX( co, st, cf, pots ) FUNCTION_CHECK( transform_K_matrix( (co), &(st), &(cf), &(pots) ), transform_K_matrix )

//------------------------------------------

#endif /* __PolyCEID_TRANSFORM_HAMILTONIAN__ */
