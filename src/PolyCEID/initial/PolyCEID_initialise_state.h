
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

#ifndef __PolyCEID_INITIALISE_STATE__
#define __PolyCEID_INITIALISE_STATE__


#include "../hamiltonians/PolyCEID_hamiltonian_many.h"
#include "../initial/PolyCEID_initial_conditions.h"
#include "../initial/PolyCEID_initialise_mu.h"
#include "../evolution/PolyCEID_transform_mu.h"
#include "../hamiltonians/PolyCEID_transform_hamiltonian.h"
#include "../evolution/PolyCEID_rho_dot.h"
#include "../output/PolyCEID_energies.h"
#include "../output/PolyCEID_adiabatic_projection.h"
#include "../evolution/PolyCEID_Ehrenfest.h"
#include "../output/PolyCEID_single_level_populations.h"
#include "../output/PolyCEID_rho_norm_and_projection.h"
#include "../output/PolyCEID_one_body_electronic_density_matrix.h"
#include "../output/PolyCEID_one_body_electronic_transition_matrix.h"
#include "../output/PolyCEID_dipole.h"
#include "../output/PolyCEID_dipole_single.h"
#include "../output/PolyCEID_ionic_density.h"
#include "../output/PolyCEID_nonadiabatic_coupling.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int initialise_state( const constants, state_p, config_p, config_p );

#define INITIALISE_STATE( co, st, cf_tmp, cf_def ) FUNCTION_CHECK( initialise_state( (co), &(st), &(cf_tmp), &(cf_def) ), initialise_state )

//------------------------------------------

#endif /* __PolyCEID_INITIALISE_STATE__ */
