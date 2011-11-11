
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

#ifndef __PolyCEID_RK2__
#define __PolyCEID_RK2__


#include "PolyCEID_rho_dot.h"
#include "PolyCEID_transform_mu.h"
#include "../hamiltonians/PolyCEID_hamiltonian_many_aux.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian.h"
#include "../output/PolyCEID_rho_norm_and_projection.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int RK2_electrons( const constants, state_p, config_p, config_p, double );

#define RK2_ELECTRONS( co, st, cf_tmp, cf_def, ratio ) FUNCTION_CHECK( RK2_electrons( (co), &(st), &(cf_tmp), &(cf_def), (ratio) ), RK2_electrons )

//------------------------------------------
//------------------------------------------

/* utilities */

//------------------------------------------

int single_step_electrons_RK2( const constants, state_p, config_p, config_p, double );

#define SINGLE_STEP_ELECTRONS_RK2( co, st, cf_tmp, cf_def, ratio ) FUNCTION_CHECK( single_step_electrons_RK2( (co), &(st), &(cf_tmp), &(cf_def), (ratio) ), single_step_electrons_RK2 )

//------------------------------------------

#endif /* __PolyCEID_RK2__ */
