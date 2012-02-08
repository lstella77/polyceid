
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

#ifndef __PolyCEID_FORCES__
#define __PolyCEID_FORCES__


#include "PolyCEID_hamiltonian_many.h"
#include "../evolution/PolyCEID_mu.h"
#include "../algebra/PolyCEID_finite_rotations.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int forces_update( const constants, state_p, config_p );

#define FORCES_UPDATE( co, st, cf ) FUNCTION_CHECK( forces_update( (co), &(st), &(cf) ), forces_update )

//------------------------------------------
//------------------------------------------

int forces_thermostat_update( const constants, state_p, config_p );

#define FORCES_THERMOSTAT_UPDATE( co, st, cf ) FUNCTION_CHECK( forces_thermostat_update( (co), &(st), &(cf) ), forces_thermostat_update )

//------------------------------------------

#endif /* __PolyCEID_FORCES__ */
