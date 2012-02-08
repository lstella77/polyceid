
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

#ifndef __PolyCEID_NORMAL_MODES__
#define __PolyCEID_NORMAL_MODES__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../initial/PolyCEID_indexing.h"
#include "../hamiltonians/PolyCEID_hamiltonian_many.h"
#include "../initial/PolyCEID_relaxation.h"
#include "../algebra/PolyCEID_finite_rotations.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* compute phonons */

int compute_phonons( const constants, state_p, config_p, config_p );

#define COMPUTE_PHONONS( co, st, cf_tmp, cf_def ) FUNCTION_CHECK( compute_phonons( (co), &(st), &(cf_tmp), &(cf_def) ), compute_phonons )

//------------------------------------------

/* utilities */

//------------------------------------------

/* compute original_hessian */

int compute_original_hessian( const constants, state_p, config_p );

#define COMPUTE_ORIGINAL_HESSIAN( co, st, cf ) FUNCTION_CHECK( compute_original_hessian( (co), &(st), &(cf) ), compute_original_hessian )

//------------------------------------------

int compute_normal_modes( const constants, state_p, config_p );

#define COMPUTE_NORMAL_MODES( co, st, cf ) FUNCTION_CHECK( compute_normal_modes( (co), &(st), &(cf) ), compute_normal_modes )

//------------------------------------------
//------------------------------------------

/*  check_symmetries */

int check_symmetries( const constants constants, state_p state_p, config_p );

#define CHECK_SYMMETRIES( co, st, cf  ) FUNCTION_CHECK( check_symmetries( (co), &(st), &(cf) ), check_symmetries )

//------------------------------------------

/* hessian_corrections */

int hessian_corrections( const constants, state_p, config_p, matrix_p );

#define HESSIAN_CORRECTIONS( co, st, cf, mat )   FUNCTION_CHECK( hessian_corrections( (co), &(st), &(cf), &(mat) ), hessian_corrections )

//------------------------------------------

/* hessian_corrections_aux */

int hessian_corrections_aux( const constants, state_p, config_p, vector, matrix_p );

#define HESSIAN_CORRECTIONS_AUX( co, st, cf, vec, mat )   FUNCTION_CHECK( hessian_corrections_aux( (co), &(st), &(cf), (vec), &(mat) ), hessian_corrections_aux )

//------------------------------------------

/* compute_phononic_dos */

int compute_phononic_dos( const constants, state_p );

#define COMPUTE_PHONONIC_DOS( co, st )   FUNCTION_CHECK( compute_phononic_dos( (co), &(st) ), compute_phononic_dos )

//------------------------------------------

#endif /* __PolyCEID_NORMAL_MODES__  */
