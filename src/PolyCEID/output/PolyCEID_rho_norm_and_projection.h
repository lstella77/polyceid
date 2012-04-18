
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

#ifndef __PolyCEID_RHO_NORM_AND_PROJECTION__
#define __PolyCEID_RHO_NORM_AND_PROJECTION__


#include "../evolution/PolyCEID_rho_dot.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int rho_norm_and_projection_update( const constants, state_p, config_p );

#define RHO_NORM_AND_PROJECTION_UPDATE( co, st, cf ) FUNCTION_CHECK( rho_norm_and_projection_update( (co), &(st), &(cf) ), rho_norm_and_projection_update )

//------------------------------------------

/* utilities */

//------------------------------------------

#ifdef __DEBUG__

/* rho_hermiticity_test */

int rho_hermiticity_test( const constants, state_p, config_p );

#define RHO_HERMITICITY_TEST( co, st, cf ) FUNCTION_CHECK( rho_hermiticity_test( (co), &(st), &(cf) ), rho_hermiticity_test ) 

//------------------------------------------

/* rho_dot_hermiticity_test */

int rho_dot_hermiticity_test( const constants, state_p, config_p );

#define RHO_DOT_HERMITICITY_TEST( co, st, cf ) FUNCTION_CHECK( rho_dot_hermiticity_test( (co), &(st), &(cf) ), rho_dot_hermiticity_test ) 

//------------------------------------------

#endif /* __DEBUG__ */

//------------------------------------------

#endif /* __PolyCEID_RHO_NORM_AND_PROJECTION__ */
