
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

#ifndef __PolyCEID_MU__
#define __PolyCEID_MU__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../structures/PolyCEID_structures_config.h"
#include "../initial/PolyCEID_indexing.h"
#include "../algebra/PolyCEID_linear_algebra.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* mu_def_update */

int mu_update( const constants, state_p, config_p );

#define MU_UPDATE( co, st, cf ) FUNCTION_CHECK( mu_update( (co), &(st), &(cf) ), mu_update )

//------------------------------------------
//------------------------------------------

/* mu00_update */

int mu00_update( const constants, state_p, config_p );

#define MU00_UPDATE( co, st, cf ) FUNCTION_CHECK( mu00_update( (co), &(st), &(cf) ), mu00_update )

//------------------------------------------

/* mu01_update */

int mu01_update( const constants, state_p, config_p );

#define MU01_UPDATE( co, st, cf ) FUNCTION_CHECK( mu01_update( (co), &(st), &(cf) ), mu01_update )

//------------------------------------------

/* mu10_update */

int mu10_update( const constants, state_p, config_p );

#define MU10_UPDATE( co, st, cf ) FUNCTION_CHECK( mu10_update( (co), &(st), &(cf) ), mu10_update )

//------------------------------------------

/* mu02_update */

int mu02_update( const constants, state_p, config_p );

#define MU02_UPDATE( co, st, cf ) FUNCTION_CHECK( mu02_update( (co), &(st), &(cf) ), mu02_update )

//------------------------------------------

/* mu11_update */

int mu11_update( const constants, state_p, config_p );

#define MU11_UPDATE( co, st, cf ) FUNCTION_CHECK( mu11_update( (co), &(st), &(cf) ), mu11_update )

//------------------------------------------
/* mu20_update */

int mu20_update( const constants, state_p, config_p );

#define MU20_UPDATE( co, st, cf ) FUNCTION_CHECK( mu20_update( (co), &(st), &(cf) ), mu20_update )

//------------------------------------------

/* utilities */

//------------------------------------------

/* irrelevant_mu_update */
int irrelevant_mu_update( const constants, state_p, config_p );

#define IRRELEVANT_MU_UPDATE( co, st, cf ) FUNCTION_CHECK( irrelevant_mu_update( (co), &(st), &(cf) ), irrelevant_mu_update )

//------------------------------------------

#endif /* __PolyCEID_MU__ */
