
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

#ifndef __PolyCEID_OUTPUT__
#define __PolyCEID_OUTPUT__


#include "../structures/PolyCEID_structures_constants.h"
#include "../structures/PolyCEID_structures_state.h"
#include "../structures/PolyCEID_structures_config.h"
#include "../initial/PolyCEID_indexing.h"
#include "../output/PolyCEID_energies.h"
#include <string.h>


/*********************
  FUNCTIONS & MACROS
*********************/

int print_observables( const constants, const state, const config );

#define PRINT_OBSERVABLES( co, st, cf ) FUNCTION_CHECK(  print_observables( (co), (st), (cf) ), print_observables )

//------------------------------------------

int print_positions( const constants, const state, const config );

#define PRINT_POSITIONS( co, st, cf ) FUNCTION_CHECK(  print_positions( (co), (st), (cf) ), print_positions )

//------------------------------------------

int print_momenta( const constants, const state, const config );

#define PRINT_MOMENTA( co, st, cf ) FUNCTION_CHECK(  print_momenta( (co), (st), (cf) ), print_momenta )

//------------------------------------------

int print_forces( const constants, const state, const config );

#define PRINT_FORCES( co, st, cf ) FUNCTION_CHECK(  print_forces( (co), (st), (cf) ), print_forces )

//------------------------------------------

int print_positions_thermostat( const constants, const state, const config );

#define PRINT_POSITIONS_THERMOSTAT( co, st, cf ) FUNCTION_CHECK(  print_positions_thermostat( (co), (st), (cf) ), print_positions )

//------------------------------------------

int print_momenta_thermostat( const constants, const state, const config );

#define PRINT_MOMENTA_THERMOSTAT( co, st, cf ) FUNCTION_CHECK(  print_momenta_thermostat( (co), (st), (cf) ), print_momenta )

//------------------------------------------

int print_forces_thermostat( const constants, const state, const config );

#define PRINT_FORCES_THERMOSTAT( co, st, cf ) FUNCTION_CHECK(  print_forces_thermostat( (co), (st), (cf) ), print_forces )

//------------------------------------------

int print_populations( const constants, const state, const config );

#define PRINT_POPULATIONS( co, st, cf ) FUNCTION_CHECK(  print_populations( (co), (st), (cf) ), print_populations )

//------------------------------------------

int print_mu_trace( const constants, const state, const config );

#define PRINT_MU_TRACE( co, st, cf ) FUNCTION_CHECK(  print_mu_trace( (co), (st), (cf) ), print_mu_trace ) 

//------------------------------------------

int print_mu_norm( const constants, const state, const config );

#define PRINT_MU_NORM( co, st, cf ) FUNCTION_CHECK(  print_mu_norm( (co), (st), (cf) ), print_mu_norm ) 

//------------------------------------------

#ifdef __DEBUG_PLUS__

int print_rho_trace( const constants, const state, const config );

#define PRINT_RHO_TRACE( co, st, cf ) FUNCTION_CHECK(  print_rho_trace( (co), (st), (cf) ), print_rho_trace ) 

//------------------------------------------

int print_rho_norm( const constants, const state, const config );

#define PRINT_RHO_NORM( co, st, cf ) FUNCTION_CHECK(  print_rho_norm( (co), (st), (cf) ), print_rho_norm ) 

#endif /* __DEBUG_PLUS__ */

//------------------------------------------

int print_adiabatic_populations( const constants, const state );

#define PRINT_ADIABATIC_POPULATIONS( co, st ) FUNCTION_CHECK(  print_adiabatic_populations( (co), (st) ), print_adiabatic_populationsn ) 

//------------------------------------------

int print_adiabatic_projection( const constants, const state );

#define PRINT_ADIABATIC_PROJECTION( co, st ) FUNCTION_CHECK(  print_adiabatic_projection( (co), (st) ), print_adiabatic_projection ) 

//------------------------------------------

int print_projection( const state );

#define PRINT_PROJECTION( st ) FUNCTION_CHECK(  print_projection( (st) ), print_projection ) 

//------------------------------------------

int print_adiabatic_PES_many( const constants, const state );

#define PRINT_ADIABATIC_PES_MANY( co, st ) FUNCTION_CHECK(  print_adiabatic_PES_many( (co), (st) ), print_adiabatic_PES_many ) 

//------------------------------------------

int print_adiabatic_PES_single( const constants, const state );

#define PRINT_ADIABATIC_PES_SINGLE( co, st ) FUNCTION_CHECK(  print_adiabatic_PES_single( (co), (st) ), print_adiabatic_PES_single ) 

//------------------------------------------

int print_single_level_populations_Ehrenfest( const constants, const state );

#define PRINT_SINGLE_LEVEL_POPULATIONS_EHRENFEST( co, st ) FUNCTION_CHECK(  print_single_level_populations_Ehrenfest( (co), (st) ), print_single_level_populations_Ehrenfest ) 

//------------------------------------------

int print_single_level_populations_adiabatic( const constants, const state );

#define PRINT_SINGLE_LEVEL_POPULATIONS_ADIABATIC( co, st ) FUNCTION_CHECK(  print_single_level_populations_adiabatic( (co), (st) ), print_single_level_populations_adiabatic ) 

//------------------------------------------

int print_nonadiabatic_coupling( const constants, const state );

#define PRINT_NONADIABATIC_COUPLING( co, st ) FUNCTION_CHECK(  print_nonadiabatic_coupling( (co), (st) ), print_nonadiabatic_coupling ) 

//------------------------------------------

int print_nonadiabatic_rate( const constants, const state );

#define PRINT_NONADIABATIC_RATE( co, st ) FUNCTION_CHECK(  print_nonadiabatic_rate( (co), (st) ), print_nonadiabatic_rate ) 

//------------------------------------------

int print_one_body_electronic_density_matrix( const constants, const state );

#define PRINT_ONE_BODY_ELECTRONIC_DENSITY_MATRIX( co, st ) FUNCTION_CHECK(  print_one_body_electronic_density_matrix( (co), (st) ), print_one_body_electronic_density_matrix ) 

//------------------------------------------

int print_one_body_electronic_hole_matrix( const constants, const state );

#define PRINT_ONE_BODY_ELECTRONIC_HOLE_MATRIX( co, st ) FUNCTION_CHECK(  print_one_body_electronic_hole_matrix( (co), (st) ), print_one_body_electronic_hole_matrix ) 

//------------------------------------------

int print_one_body_electronic_particle_matrix( const constants, const state );

#define PRINT_ONE_BODY_ELECTRONIC_PARTICLE_MATRIX( co, st ) FUNCTION_CHECK(  print_one_body_electronic_particle_matrix( (co), (st) ), print_one_body_electronic_particle_matrix ) 

//------------------------------------------

int print_natural_orbitals( const constants, state );

#define PRINT_NATURAL_ORBITALS( co, st ) FUNCTION_CHECK(  print_natural_orbitals( (co), (st) ), print_natural_orbitals ) 

//------------------------------------------

int print_hole_orbitals( const constants, const state );

#define PRINT_HOLE_ORBITALS( co, st ) FUNCTION_CHECK(  print_hole_orbitals( (co), (st) ), print_hole_orbitals ) 

//------------------------------------------

int print_particle_orbitals( const constants, const state );

#define PRINT_PARTICLE_ORBITALS( co, st ) FUNCTION_CHECK(  print_particle_orbitals( (co), (st) ), print_particle_orbitals ) 

//------------------------------------------

int print_adiabatic_states( const constants, const state );

#define PRINT_ADIABATIC_STATES( co, st ) FUNCTION_CHECK(  print_adiabatic_states( (co), (st) ), print_adiabatic_states ) 

//------------------------------------------

int print_electronic_density_states( const constants, const state );

#define PRINT_ELECTRONIC_DENSITY_STATES( co, st ) FUNCTION_CHECK(  print_electronic_density_states( (co), (st) ), print_electronic_density_states ) 

//------------------------------------------

int print_ionic_density_states( const constants, const state );

#define PRINT_IONIC_DENSITY_STATES( co, st ) FUNCTION_CHECK(  print_ionic_density_states( (co), (st) ), print_ionic_density_states ) 

//------------------------------------------

int print_dipole_many( const constants, const state );

#define PRINT_DIPOLE_MANY( co, st ) FUNCTION_CHECK(  print_dipole_many( (co), (st) ), print_dipole_many ) 

//------------------------------------------

int print_dipole_single( const constants, const state );

#define PRINT_DIPOLE_SINGLE( co, st ) FUNCTION_CHECK(  print_dipole_single( (co), (st) ), print_dipole_single ) 

//------------------------------------------

int print_position_ion_frame( const constants, const state, const config );

#define PRINT_FRAME( co, st, cf )   FUNCTION_CHECK(  print_position_ion_frame( (co), (st), (cf) ), print_position_ion_frame  )

//------------------------------------------

int print_energies( const constants, const state );

#define PRINT_ENERGIES( co, st ) FUNCTION_CHECK(  print_energies( (co), (st) ), print_energies ) 

//------------------------------------------

/* utilities */

//------------------------------------------

int output_files_opening( const constants );

#define OUTPUT_OPENING( co ) FUNCTION_CHECK(  output_files_opening( (co) ), output_files_opening ) 

//------------------------------------------

int output_files_closing( const constants );

#define OUTPUT_CLOSING( co ) FUNCTION_CHECK(  output_files_closing( (co) ), output_files_closing ) 

//------------------------------------------

#endif /* __PolyCEID_OUTPUT__  */
