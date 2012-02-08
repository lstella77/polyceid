
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

#ifndef __PolyCEID_FINITE_ROTATIONS__
#define __PolyCEID_FINITE_ROTATIONS__

#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../initial/PolyCEID_indexing.h"
#include "../initial/PolyCEID_relaxation.h"
#include "../utils/PolyCEID_utilities.h"
#include "../algebra/PolyCEID_linear_algebra.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* check_finite_rotation_hamiltonian */

int check_finite_rotation_hamiltonian( const constants, state_p, config_p, int );

#define CHECK_FINITE_ROTATION_HAMILTONIAN( co, st, cf, angle ) FUNCTION_CHECK( check_finite_rotation_hamiltonian( (co), &(st), &(cf), (angle) ), check_finite_rotation_hamiltonian )

//------------------------------------------

/* check_finite_rotation_forces */

int check_finite_rotation_forces( const constants, state_p, config_p, int );

#define CHECK_FINITE_ROTATION_FORCES( co, st, cf, angle ) FUNCTION_CHECK( check_finite_rotation_forces( (co), &(st), &(cf), (angle) ), check_finite_rotation_forces )

//------------------------------------------

/* check_finite_rotation_positions */

int check_finite_rotation_positions( const constants, state_p, config_p, int );

#define CHECK_FINITE_ROTATION_POSITIONS( co, st, cf, angle ) FUNCTION_CHECK( check_finite_rotation_positions( (co), &(st), &(cf), (angle) ), check_finite_rotation_positions )

//------------------------------------------

/* check_finite_rotation_momenta */

int check_finite_rotation_momenta( const constants, state_p, config_p, int );

#define CHECK_FINITE_ROTATION_MOMENTA( co, st, cf, angle ) FUNCTION_CHECK( check_finite_rotation_momenta( (co), &(st), &(cf), (angle) ), check_finite_rotation_momenta )


//------------------------------------------

/* utilities */

//------------------------------------------

/* check_finite_rotation_vector */

int check_finite_rotation_vector( const constants, int, vector_p );

#define CHECK_FINITE_ROTATION_VECTOR( c, angle, vec )  FUNCTION_CHECK( check_finite_rotation_vector( (c), (angle), &(vec) ), check_finite_rotation_vector )

//------------------------------------------

/* check_finite_double_rotation_vector */

int check_finite_double_rotation_vector( const constants, int, vector_p );

#define CHECK_FINITE_DOUBLE_ROTATION_VECTOR( c, angle, vec )  FUNCTION_CHECK( check_finite_double_rotation_vector( (c), (angle), &(vec) ), check_finite_double_rotation_vector )

//------------------------------------------

/* check_finite_rotation_matrix */

int check_finite_rotation_matrix( const constants, int, matrix_p );

#define CHECK_FINITE_ROTATION_MATRIX( c, angle, mat )   FUNCTION_CHECK( check_finite_rotation_matrix( (c), (angle), &(mat) ), check_finite_rotation_matrix )

//------------------------------------------
//------------------------------------------

/* check_finite_rotation_vector_comp */

int check_finite_rotation_vector_comp( const constants, int, vector_p );

#define CHECK_FINITE_ROTATION_VECTOR_COMP( c, angle, vec )  FUNCTION_CHECK( check_finite_rotation_vector_comp( (c), (angle), &(vec) ), check_finite_rotation_vector_comp )

//------------------------------------------

/* check_finite_double_rotation_vector_comp */

int check_finite_double_rotation_vector_comp( const constants, int, vector_p );

#define CHECK_FINITE_DOUBLE_ROTATION_VECTOR_COMP( c, angle, vec )  FUNCTION_CHECK( check_finite_double_rotation_vector_comp( (c), (angle), &(vec) ), check_finite_double_rotation_vector_comp )

//------------------------------------------

/* check_finite_rotation_matrix_comp */

int check_finite_rotation_matrix_comp( const constants, int, matrix_p );

#define CHECK_FINITE_ROTATION_MATRIX_COMP( c, angle, mat )   FUNCTION_CHECK( check_finite_rotation_matrix_comp( (c), (angle), &(mat) ), check_finite_rotation_matrix_comp )

//------------------------------------------
//------------------------------------------

/* compute_finite_rotation */

int compute_finite_rotation( const constants, int , matrix_p );

#define COMPUTE_FINITE_ROTATION( c, angle, mat )   FUNCTION_CHECK( compute_finite_rotation( (c), (angle) , &(mat) ), compute_finite_rotation )

//------------------------------------------

/* compute_finite_rotation_comp */

int compute_finite_rotation_comp( const constants, int , matrix_p );

#define COMPUTE_FINITE_ROTATION_COMP( c, angle, mat )   FUNCTION_CHECK( compute_finite_rotation_comp( (c), (angle) , &(mat) ), compute_finite_rotation_comp )

//------------------------------------------
//------------------------------------------

/* centre_of_mass_transform_vector */

int centre_of_mass_transform_vector( const constants, vector_p );

#define CENTRE_OF_MASS_TRANSFORM_VECTOR( c, vec )   FUNCTION_CHECK( centre_of_mass_transform_vector( (c), &(vec) ), centre_of_mass_transform_vector )

//------------------------------------------
//------------------------------------------

/* symmetry_adapted_unitary_transform  */

int symmetry_adapted_unitary_transform( constants, double, matrix_p );

#define SYMMETRY_ADAPTED_UNITARY_TRANSFORM( c, angle, mat )  FUNCTION_CHECK( symmetry_adapted_unitary_transform( (c), (angle), &(mat) ), symmetry_adapted_unitary_transform )

//------------------------------------------

#endif /* __PolyCEID_FINITE_ROTATIONS__ */
