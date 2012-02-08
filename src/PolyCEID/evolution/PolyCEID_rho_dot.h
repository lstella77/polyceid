
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

#ifndef __PolyCEID_RHO_DOT__
#define __PolyCEID_RHO_DOT__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../initial/PolyCEID_indexing.h"
#include "../hamiltonians/PolyCEID_hamiltonian_many.h"
#include "../algebra/PolyCEID_linear_algebra.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* rho_dot_update */

int rho_dot_update( const constants, state_p, config_p );

#define RHO_DOT_UPDATE( co, st, cf ) FUNCTION_CHECK(  rho_dot_update( (co), &(st), &(cf) ), rho_dot_update  )

//------------------------------------------

/* utilites */

//------------------------------------------

#ifdef __DEBUG__

/* rho_dot_is_zero */

int rho_dot_is_zero( const constants, state_p, config_p );

#define RHO_DOT_IS_ZERO( co, st, cf ) FUNCTION_CHECK( rho_dot_is_zero( (co), &(st), &(cf) ), rho_dot_is_zero )

//------------------------------------------

#endif /* __DEBUG__ */

//------------------------------------------

#endif /* __PolyCEID_RHO_DOT__ */
