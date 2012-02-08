
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

#ifndef __PolyCEID_EHRENFEST__
#define __PolyCEID_EHRENFEST__


#include "../algebra/PolyCEID_linear_algebra.h"
#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../hamiltonians/PolyCEID_hamiltonian_single.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int initialise_Ehrenfest_frame( const constants, state_p, config_p );

#define INITIALISE_EHRENFEST_FRAME( co, st, cf )  FUNCTION_CHECK( initialise_Ehrenfest_frame( (co), &(st), &(cf) ), initialise_Ehrenfest_frame )

//------------------------------------------
//------------------------------------------

int Ehrenfest_frame_update( const constants, state_p, config_p, config_p, double );

#define EHRENFEST_FRAME_UPDATE( co, st, cf_def, cf_tmp, ratio ) FUNCTION_CHECK( Ehrenfest_frame_update( (co), &(st), &(cf_def), &(cf_tmp), (ratio) ), Ehrenfest_frame_update )

//------------------------------------------
//------------------------------------------

/* utilities */

//------------------------------------------

int exact_Ehrenfest( const constants, state_p, config_p, config_p, double );

#define EXACT_EHRENFEST( co, st, cf_def, cf_tmp, ratio )  FUNCTION_CHECK( exact_Ehrenfest( (co), &(st), &(cf_def), &(cf_tmp), (ratio) ), exact_Ehrenfest )

//------------------------------------------

int orthonormalise_Ehrenfest_frame( const constants, state_p, config_p );

#define ORTHONORMALISE_EHRENFEST_FRAME( co, st, cf ) FUNCTION_CHECK( orthonormalise_Ehrenfest_frame( (co), &(st), &(cf) ), orthonormalise_Ehrenfest_frame )

//------------------------------------------

#endif /* __PolyCEID_EHRENFEST__ */
