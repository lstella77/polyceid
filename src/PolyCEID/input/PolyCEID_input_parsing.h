
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

#ifndef __PolyCEID_INPUT_PARSING__
#define __PolyCEID_INPUT_PARSING__


#include "../input/PolyCEID_parsing.h"
#include "../input/PolyCEID_list.h"
#include "../input/PolyCEID_assign_variable.h"
#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_atoms.h"
#include "../initial/PolyCEID_indexing.h"
#include "../initial/PolyCEID_initial_rho_electron.h"
#include "../initial/PolyCEID_CI_table.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* INPUT PARSING */

int input_parsing( FILE*, constants_p constants_p, int, rvector );

#define INPUT_PARSING( fp, co, counter, rvec ) FUNCTION_CHECK( input_parsing( (fp), &(co), (counter), (rvec) ), input_parsing )

//------------------------------------------

int hamiltonian_parsing( list_p list_p, constants_p constants_p );

#define HAMILTONIAN_PARSING( lp, co ) FUNCTION_CHECK( hamiltonian_parsing( (lp), &(co) ), hamiltonian_parsing )

//------------------------------------------

int assign_initial_condition_atoms( list_p, constants_p );

#define ASSIGN_INITIAL_CONDITION_ATOMS( lp, co ) FUNCTION_CHECK( assign_initial_condition_atoms( (lp), &(co) ), assign_initial_condition_atoms )

//------------------------------------------

int initial_condition_atoms_parsing( constants_p, int, rvector );

#define INITIAL_CONDITION_ATOMS_PARSING( co,ocounter, rvec ) FUNCTION_CHECK( initial_condition_atoms_parsing( &(co), (counter), (rvec) ), initial_condition_atoms_parsing )

//------------------------------------------

int initial_condition_electrons_parsing( constants_p constants_p );

#define INITIAL_CONDITION_ELECTRONS_PARSING( co ) FUNCTION_CHECK( initial_condition_electrons_parsing( &(co) ), initial_condition_electrons_parsing )

//------------------------------------------

int construct_transition( const constants constants, ivector, rvector_p transition_p );

#define CONSTRUCT_TRANSITION( co, occ, rvec ) FUNCTION_CHECK( construct_transition( (co), (occ), &(rvec) ), construct_transition ) 

//------------------------------------------

int output_parsing( list_p, constants_p );

#define OUTPUT_PARSING( lp, co ) FUNCTION_CHECK( output_parsing( (lp), &(co) ), output_parsing )

//------------------------------------------

/* utilites */

//------------------------------------------

#endif /* __PolyCEID_INPUT_PARSING__ */
