
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

#ifndef __PolyCEID_STRUCTURES_THERMOSTAT__
#define __PolyCEID_STRUCTURES_THERMOSTAT__


#include "../algebra/PolyCEID_constants.h"
#include "../algebra/PolyCEID_linear_algebra.h"
#include <stdio.h>


/* parameters */
#define NP_THERMOSTAT  1


/*********************
  THERMOSTAT
*********************/

typedef struct{

  rvector  positions;       // [N_chain]           vectorised!
  rvector  momenta;         // [N_chain]           vectorised!
  rvector  forces;          // [N_chain]           vectorised!

} PolyCEID_thermostat;

/*********************
  POINTERS
*********************/

typedef PolyCEID_thermostat  thermostat;

typedef PolyCEID_thermostat* thermostat_p;


/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_thermostat_allocate( int, int*, thermostat_p );

#define THERMOSTAT_ALLOCATE( np, pp, a )  FUNCTION_CHECK( PolyCEID_thermostat_allocate( (np), (pp), &(a) ), PolyCEID_thermostat_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_thermostat_free( thermostat_p );

#define THERMOSTAT_FREE( a )  FUNCTION_CHECK( PolyCEID_thermostat_free( &(a) ), PolyCEID_thermostat_free )

//------------------------------------------

/* READ */

int PolyCEID_thermostat_read( FILE*, thermostat_p );

#define THERMOSTAT_READ( fp, a )  FUNCTION_CHECK( PolyCEID_thermostat_read( fp, &(a) ), PolyCEID_thermostat_read )

//------------------------------------------

/* PRINT */

int PolyCEID_thermostat_print( FILE*, const thermostat );

#define THERMOSTAT_PRINT( fp, a )  FUNCTION_CHECK( PolyCEID_thermostat_print( (fp), (a) ), PolyCEID_thermostat_print )

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_thermostat_verbose_print( FILE*, const thermostat );

#define THERMOSTAT_VERBOSE_PRINT( fp, a )  FUNCTION_CHECK( PolyCEID_thermostat_verbose_print( (fp), (a) ), PolyCEID_thermostat_verbose_print )

//------------------------------------------


/* COPY [ a1 = a2 ] */

int PolyCEID_thermostat_copy( thermostat_p, const thermostat );

#define THERMOSTAT_COPY( a1, a2 )  FUNCTION_CHECK( PolyCEID_thermostat_copy( &(a1), (a2) ), PolyCEID_thermostat_copy )

//------------------------------------------

/* COMPARE [ a1 = a2 ] */

int PolyCEID_thermostat_compare( const thermostat, const thermostat );

#define THERMOSTAT_COMPARE( a1, a2 )  FUNCTION_CHECK( PolyCEID_thermostat_compare( (a1), (a2) ), PolyCEID_thermostat_compare )

//------------------------------------------

#endif /* __PolyCEID_STRUCTURES_THERMOSTAT__ */

