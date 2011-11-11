
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

#ifndef __PolyCEID_HAMILTONIAN_MANY_AUX__
#define __PolyCEID_HAMILTONIAN_MANY_AUX__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../hamiltonians/PolyCEID_distances.h"
#include "../initial/PolyCEID_indexing.h"
#include "../hamiltonians/PolyCEID_forces.h"
#include "../hamiltonians/PolyCEID_hamiltonian_single.h"
#include "../hamiltonians/PolyCEID_hamiltonian_many.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* hamiltonian_many_aux update */

int hamiltonian_many_aux_update( const constants, state_p, config_p );

#define HAMILTONIAN_MANY_AUX_UPDATE( co, st, cf ) FUNCTION_CHECK( hamiltonian_many_aux_update( (co), &(st), &(cf) ), hamiltonian_many_aux_update ) 

//------------------------------------------

#endif /* __PolyCEID_HAMILTONIAN_MANY_AUX__ */
