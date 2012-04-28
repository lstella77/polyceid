
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

#ifndef  __PolyCEID_START_FILE__
#define  __PolyCEID_START_FILE__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../structures/PolyCEID_structures_config.h"

//------------------------------------------

/* READ */

int read_start_file( constants_p, state_p, config_p );

#define READ_START_FILE( co, st, cf ) FUNCTION_CHECK( read_start_file( &(co), &(st), &(cf) ), read_start_file )


//------------------------------------------

/* PRINT */

int print_start_file( const constants, const state, const config );

#define PRINT_START_FILE( co, st, cf ) FUNCTION_CHECK( print_start_file( (co), (st), (cf) ), print_start_file )

//------------------------------------------


#endif /* __PolyCEID_START_FILE__  */
