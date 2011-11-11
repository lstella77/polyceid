
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

#ifndef __PolyCEID_INTEGRATOR__
#define __PolyCEID_INTEGRATOR__


#include "PolyCEID_verlet.h"

#ifndef __RK4__

#include "PolyCEID_RK2.h"

#else /* __RK4__ */

#include "PolyCEID_RK4.h"

#endif /* __RK4__ */

#include "../output/PolyCEID_output.h"
#include "../output/PolyCEID_energies.h"
#include "../evolution/PolyCEID_rho_dot.h"
#include "../evolution/PolyCEID_transform_mu.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian.h"
#include "../output/PolyCEID_adiabatic_projection.h"
#include "../algebra/PolyCEID_finite_rotations.h"
#include "../hamiltonians/PolyCEID_hamiltonian_many_aux.h"
#include "../output/PolyCEID_rho_norm_and_projection.h"
#include "../output/PolyCEID_single_level_populations.h"
#include "../initial/PolyCEID_initialise_state.h"
#include "../output/PolyCEID_one_body_electronic_density_matrix.h"
#include "../output/PolyCEID_one_body_electronic_transition_matrix.h"
#include "../output/PolyCEID_dipole.h"
#include "../output/PolyCEID_dipole_single.h"
#include "../output/PolyCEID_ionic_density.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int global_CEID_integrator( const constants, state_p, config_p, config_p );

#define GLOBAL_CEID_INTEGRATOR( co, st, cf_tmp, cf_def ) FUNCTION_CHECK( global_CEID_integrator( (co), &(st), &(cf_tmp), &(cf_def) ), global_CEID_integrator )

//------------------------------------------

/* utilites */

//------------------------------------------

int single_step_CEID_integrator( const constants, state_p, config_p, config_p );

#define SINGLE_STEP_CEID_INTEGRATOR( co, st, cf_tmp, cf_def ) FUNCTION_CHECK( single_step_CEID_integrator( (co), &(st), &(cf_tmp), &(cf_def) ), single_step_CEID_integrator )

//------------------------------------------

int compute_observables( const constants, state_p, config_p );
#define COMPUTE_OBSERVABLES( co, st, cf ) FUNCTION_CHECK( compute_observables( (co), &(st), &(cf) ), compute_observables )

//------------------------------------------

int global_update( const constants, state_p, config_p, config_p );
#define GLOBAL_UPDATE( co, st, cf_def, cf_tmp ) FUNCTION_CHECK( global_update( (co), &(st), &(cf_def), &(cf_tmp) ), global_update )

//------------------------------------------

#endif /* __PolyCEID_INTEGRATOR__ */
