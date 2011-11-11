
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

#ifndef __PolyCEID_TRANSFORM_MU__
#define __PolyCEID_TRANSFORM_MU__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../structures/PolyCEID_structures_config.h"
#include "../initial/PolyCEID_indexing.h"
#include "../algebra/PolyCEID_linear_algebra.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* transform_mu */

int transform_mu( constants, state_p, config_p );

#define TRANSFORM_MU( co, st, cf ) FUNCTION_CHECK( transform_mu( (co), &(st), &(cf) ), transform_mu )

//------------------------------------------
//------------------------------------------

/* transform_mu01 */

int transform_mu01( constants, state_p, config_p );

#define TRANSFORM_MU01( co, st, cf ) FUNCTION_CHECK( transform_mu01( (co), &(st), &(cf) ), transform_mu01 )

//------------------------------------------

/* transform_mu10 */

int transform_mu10( constants, state_p, config_p );

#define TRANSFORM_MU10( co, st, cf ) FUNCTION_CHECK( transform_mu10( (co), &(st), &(cf) ), transform_mu10f )

//------------------------------------------

/* transform_mu02 */

int transform_mu02( constants, state_p, config_p );

#define TRANSFORM_MU02( co, st, cf ) FUNCTION_CHECK( transform_mu02( (co), &(st), &(cf) ), transform_mu02 )

//------------------------------------------

/* transform_mu11 */

int transform_mu11( constants, state_p, config_p );

#define TRANSFORM_MU11( co, st, cf ) FUNCTION_CHECK( transform_mu11( (co), &(st), &(cf) ), transform_mu11 )

//------------------------------------------

/* transform_mu20 */

int transform_mu20( constants, state_p, config_p );

#define TRANSFORM_MU20( co, st, cf ) FUNCTION_CHECK( transform_mu20( (co), &(st), &(cf) ), transform_mu20 )

//------------------------------------------

#endif /* __PolyCEID_TRANSFORM_MU__ */
