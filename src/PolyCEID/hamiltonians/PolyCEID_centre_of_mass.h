
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

#ifndef __PolyCEID_CENTRE_OF_MASS__
#define __PolyCEID_CENTRE_OF_MASS__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../structures/PolyCEID_structures_config.h"
#include "../utils/PolyCEID_utilities.h"


/* compute_centre_of_mass */

int compute_centre_of_mass( const constants, atoms_p );

#define COMPUTE_CENTRE_OF_MASS( co, at ) FUNCTION_CHECK( compute_centre_of_mass( (co), &(at) ), compute_centre_of_mass )

//------------------------------------------

/* centre_of_mass_transform */

int centre_of_mass_transform( const constants, atoms_p );

#define CENTRE_OF_MASS_TRANSFORM( co, at )   FUNCTION_CHECK( centre_of_mass_transform( (co), &(at) ), centre_of_mass_transform )

//------------------------------------------

#endif /* __PolyCEDID_CENTRE_OF_MASS__ */
