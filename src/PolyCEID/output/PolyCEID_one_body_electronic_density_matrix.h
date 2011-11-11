
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

#ifndef __PolyCEID_ONE_BODY_DENSITY_MATRIX__
#define __PolyCEID_ONE_BODY_DENSITY_MATRIX__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../structures/PolyCEID_structures_config.h"
#include "../initial/PolyCEID_indexing.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int compute_one_body_electronic_density_matrix( const constants, state_p, config_p );

#define COMPUTE_ONE_BODY_ELECTRONIC_DENSITY_MATRIX( co, st, cf ) FUNCTION_CHECK( compute_one_body_electronic_density_matrix( (co), &(st), &(cf) ), compute_one_body_electronic_density_matrix )

//------------------------------------------

/* utilities */

//------------------------------------------

int compute_one_body_electronic_density_matrix_Ehrenfest( const constants, state_p, config_p );

#define COMPUTE_ONE_BODY_ELECTRONIC_DENSITY_MATRIX_EHRENFEST( co, st, cf ) FUNCTION_CHECK( compute_one_body_electronic_density_matrix_Ehrenfest( (co), &(st), &(cf) ), compute_one_body_electronic_density_matrix_Ehrenfest )

//------------------------------------------

#endif /* __PolyCEID_ONE_BODY_DENSITY_MATRIX__ */
