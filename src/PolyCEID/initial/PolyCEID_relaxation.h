
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

#ifndef __PolyCEID_RELAXATION__
#define __PolyCEID_RELAXATION__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../hamiltonians/PolyCEID_hamiltonian_many.h"
#include "../initial/PolyCEID_initial_rho_electron.h"
#include "../initial/PolyCEID_normal_modes.h"
#include "../evolution/PolyCEID_Ehrenfest.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int relaxation( const constants, state_p, config_p );

#define RELAXATION( co, st, cf ) FUNCTION_CHECK( relaxation( (co), &(st), &(cf) ), relaxation )

//------------------------------------------

int one_step_relaxation( const constants, state_p, config_p, unsigned short int );

#define ONE_STEP_RELAXATION( co, st, cf, flag ) FUNCTION_CHECK( one_step_relaxation( (co), &(st), &(cf), (flag) ), one_step_relaxation )

//------------------------------------------

/* utilities */

int allocate_one_step_relaxation_workspace( const constants );

#define ALLOCATE_ONE_STEP_RELAXATION_WORKSPACE( co ) FUNCTION_CHECK( allocate_one_step_relaxation_workspace( (co) ), allocate_one_step_relaxation_workspace )

//------------------------------------------

int free_one_step_relaxation_workspace( void );

#define FREE_ONE_STEP_RELAXATION_WORKSPACE() FUNCTION_CHECK( free_one_step_relaxation_workspace(), free_one_step_relaxation_workspace )

//------------------------------------------
//------------------------------------------

int compute_F_partial( const constants, state_p, config_p, vector_p );

#define COMPUTE_F_PARTIAL( co, st, cf, vec ) FUNCTION_CHECK( compute_F_partial( (co), &(st), &(cf), &(vec) ), compute_F_partial )

//------------------------------------------

int compute_K_partial( const constants, state_p, config_p, matrix_p );

#define COMPUTE_K_PARTIAL( co, st, cf, mat ) FUNCTION_CHECK( compute_K_partial( (co), &(st), &(cf), &(mat) ), compute_K_partial )

//------------------------------------------

#endif /* __PolyCEID_RELAXATION__  */

