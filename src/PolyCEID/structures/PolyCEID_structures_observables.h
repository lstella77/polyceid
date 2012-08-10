
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

#ifndef __PolyCEID_STRUCTURES_OBSERVABLES__
#define __PolyCEID_STRUCTURES_OBSERVABLES__


#include "../algebra/PolyCEID_linear_algebra.h"
#include "../structures/PolyCEID_structures_constants.h"
#include <stdio.h>


/* parameters */
#define NP_OBSERVABLES  5


/*********************
  OBSERVABLES
*********************/

typedef struct{

  rvector  kinetic_energy_atoms;        //[3*N_atoms]
  
  rvector  kinetic_energy_system;       //[4]
  rvector  potential_energy_system;     //[4]
  
  rvector  kinetic_energy_thermostat;   //[N_chain+1]
  rvector  potential_energy_thermostat; //[N_chain+1]
  
  double   total_energy_system;
  
  double   pseudo_energy_system;

  double   superposition_instantaneous_and_initial_state;
  double   superposition_instantaneous_and_excited_state;

  rvector  oscillator_strength_many;    // [N_levels_many *(N_levels_many-1)/2] 
  rvector  oscillator_frequency_many;   // [N_levels_many *(N_levels_many -1)/2]

  rvector  oscillator_strength_single;  // [N_levels_single]
  rvector  oscillator_frequency_single; // [N_levels_single]

} PolyCEID_observables;


/*********************
  POINTERS
*********************/

typedef PolyCEID_observables  observables;

typedef PolyCEID_observables* observables_p;


/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_observables_allocate( int, int*, observables_p );

#define OBSERVABLES_ALLOCATE( np, pp, o )  FUNCTION_CHECK( PolyCEID_observables_allocate( (np), (pp), &(o) ), PolyCEID_observables_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_observables_free( observables_p );

#define OBSERVABLES_FREE( o )  FUNCTION_CHECK( PolyCEID_observables_free( &(o) ), PolyCEID_observables_free )

//------------------------------------------

/* READ */

int PolyCEID_observables_read( FILE*, observables_p );

#define OBSERVABLES_READ( fp, o )  FUNCTION_CHECK( PolyCEID_observables_read( fp, &(o) ), PolyCEID_observables_read )

//------------------------------------------

/* PRINT */

int PolyCEID_observables_print( FILE*, const observables );

#define OBSERVABLES_PRINT( fp, o )  FUNCTION_CHECK( PolyCEID_observables_print( (fp), (o) ), PolyCEID_observables_print )

//------------------------------------------

/* PRINT VERBOSE */

int PolyCEID_observables_verbose_print( FILE*, const observables );

#define OBSERVABLES_VERBOSE_PRINT( fp, o )  FUNCTION_CHECK( PolyCEID_observables_verbose_print( (fp), (o) ), PolyCEID_observables_verbose_print )

//------------------------------------------

/* COPY [ o1 = o2 ] */

int PolyCEID_observables_copy( observables_p, const observables );

#define OBSERVABLES_COPY( o1, o2 )  FUNCTION_CHECK( PolyCEID_observables_copy( &(o1), (o2) ), PolyCEID_observables_copy )

//------------------------------------------


#endif /* __PolyCEID_STRUCTURES_OBSERVABLES__ */

