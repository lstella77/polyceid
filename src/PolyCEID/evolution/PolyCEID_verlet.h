
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

#ifndef __PolyCEID_VERLET__
#define __PolyCEID_VERLET__


#include "../hamiltonians/PolyCEID_hamiltonian_many_aux.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int velocity_verlet_ions( const constants, state_p, config_p, config_p, double );

#define VELOCITY_VERLET_IONS( co, st, cf_tmp, cf_def, ratio ) FUNCTION_CHECK( velocity_verlet_ions( (co), &(st), &(cf_tmp), &(cf_def), (ratio) ), velocity_verlet_ions )

//------------------------------------------

int nose_hoover_chain( const constants, state_p, config_p, double );

#define NOSE_HOOVER_CHAIN( co, st, cf, ratio ) FUNCTION_CHECK( nose_hoover_chain( (co), &(st), &(cf), (ratio) ), nose_hoover_chain )

//------------------------------------------

#endif /* __PolyCEID_VERLET__ */
