
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

#ifndef __PolyCEID_STRUCTURES_CONFIG__
#define __PolyCEID_STRUCTURES_CONFIG__


#include <stdio.h>
#include "../structures/PolyCEID_structures_electrons.h"
#include "../structures/PolyCEID_structures_atoms.h"
#include "../structures/PolyCEID_structures_thermostat.h"


/* parameters */
#define NP_CONFIG  7


/*********************
  CONFIG
*********************/

typedef struct{

  atoms        atoms; 
  //
  thermostat   thermostat; 
  //
  electrons    electrons;
  //

} PolyCEID_config;


/*********************
  POINTERS
*********************/

typedef PolyCEID_config  config;

typedef PolyCEID_config* config_p;


/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_config_allocate( int, int*, config_p );

#define CONFIG_ALLOCATE( np, pp, cf )  FUNCTION_CHECK( PolyCEID_config_allocate( (np), (pp), &(cf) ), PolyCEID_config_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_config_free( config_p );

#define CONFIG_FREE( cf )  FUNCTION_CHECK( PolyCEID_config_free( &(cf) ), PolyCEID_config_free )

//------------------------------------------

/* READ */

int PolyCEID_config_read( FILE*, config_p );

#define CONFIG_READ( fp, cf )  FUNCTION_CHECK( PolyCEID_config_read( fp, &(cf) ), PolyCEID_config_read )

//------------------------------------------

/* PRINT */

int PolyCEID_config_print( FILE*, const config );

#define CONFIG_PRINT( fp, cf )  FUNCTION_CHECK( PolyCEID_config_print( (fp), (cf) ), PolyCEID_config_print )

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_config_verbose_print( FILE*, const config );

#define CONFIG_VERBOSE_PRINT( fp, cf )  FUNCTION_CHECK( PolyCEID_config_verbose_print( (fp), (cf) ), PolyCEID_config_verbose_print )

//------------------------------------------


/* COPY [ cf1 = cf2 ] */

int PolyCEID_config_copy( config_p, const config );

#define CONFIG_COPY( cf1, cf2 )  FUNCTION_CHECK( PolyCEID_config_copy( &(cf1), (cf2) ), PolyCEID_config_copy )

//------------------------------------------

/* COMPARE [ a1 = a2 ] */

int PolyCEID_config_compare( const config, const config );

#define CONFIG_COMPARE( cf1, cf2 )  FUNCTION_CHECK( PolyCEID_config_compare( (cf1), (cf2) ), PolyCEID_config_compare )

//------------------------------------------

#endif /* __PolyCEID_STRUCTURES_CONFIG__ */

