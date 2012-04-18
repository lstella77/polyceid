
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

#include "config.h"
#include "PolyCEID_initialise_state.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int initialise_state( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  int  sdim;
  int  N_coor;
  int  N_levels_single;
  int  N_levels_many;
  int  max_rho_index;
  int  sqrt_max_rho_index;
  int  N_atoms;
  int  N_chain;
  /* dummies */
  int  ppc[ NP_CONFIG ];
  int  pps[ NP_STATE ];
  int  info=0;


  flag_normal_mode_expansion  =  constants.flag_normal_mode_expansion;
  N_coor                      = constants.N_coor;
  N_levels_single             = constants.N_levels_single;
  N_levels_many               = constants.N_levels_many;
  max_rho_index               = constants.max_rho_index;
  sqrt_max_rho_index          = constants.sqrt_max_rho_index;
  N_atoms                     = constants.N_atoms;
  N_chain                     = constants.N_chain;
  sdim                        = constants.spacial_dimension;

  /* parameters */
  ppc[0] = N_coor;
  ppc[1] = N_levels_single;
  ppc[2] = N_levels_many;
  ppc[3] = max_rho_index;
  ppc[4] = N_atoms;
  ppc[5] = N_chain;
  ppc[6] = sdim;

  /* config_tmp allocation */
  if( CONFIG_ALLOCATE( NP_CONFIG, ppc, *config_tmp_p ) ) info=1;

  /* config_def allocation */
  if( CONFIG_ALLOCATE( NP_CONFIG, ppc, *config_def_p ) ) info=1;


  /* parameters */
  pps[0] = N_coor;
  pps[1] = N_levels_single;
  pps[2] = N_levels_many;
  pps[3] = max_rho_index;
  pps[4] = sqrt_max_rho_index;
  pps[5] = N_atoms;
  pps[6] = N_chain;
  pps[7] = sdim;

  /* state allocation */
  if( STATE_ALLOCATE( NP_STATE, pps, *state_p ) ) info=1;


  // CONFIG_VERBOSE_PRINT( stdout, *config_tmp_p );

  // CONFIG_VERBOSE_PRINT( stdout, *config_def_p );

  // STATE_VERBOSE_PRINT( stdout, *state_p );


  /* other state variables */
  if( !(constants.flag_restart) ){

    /* step_counter already initialised */

    /* time */
    config_def_p->time = 0;

  }
  else{

    fprintf( stderr, "ERROR: restarting is not coded yet!\n" );
    fflush( stderr );

  }

  /* copy initial configurations */

  if( !info ){

    /* initial positions etc. */
    if( ATOMS_COPY( config_def_p->atoms, constants.initial_atoms ) ) info=1;

    if( ATOMS_COPY( config_tmp_p->atoms, constants.initial_atoms ) ) info=1;

  }


  if( !info ){

    /* compute_phonos */
    if( COMPUTE_PHONONS( constants, *state_p, *config_tmp_p, *config_def_p ) ) info=1; 

  }


  if( !info ){

    /* initial condition for the rho matrix */
    if( INITIAL_CONDITION_RHO( constants, *state_p, *config_def_p ) ) info=1; 

  }


  if( !info ){

    /* initialise mu */
    if( INITIALISE_MU( constants, *state_p, *config_def_p ) ) info=1;

  }


  if( !info ){

    /* initialise distances */
    if( DISTANCES_UPDATE( constants, *state_p, *config_def_p ) ) info=1;
  
  }


  if( !info ){

    /* initilise_Ehrenfest_frame_def */
    if( INITIALISE_EHRENFEST_FRAME( constants, *state_p, *config_def_p ) ) info=1;

  }


  if( !info ){

    /* initialise hamiltonian and forces */
    if( HAMILTONIAN_MANY_UPDATE( constants, *state_p, *config_def_p ) ) info=1;

  }


  if( !info && N_chain ){

    /* initialise thermostat forces */
    if( FORCES_THERMOSTAT_UPDATE( constants, *state_p, *config_def_p ) ) info=1;

  }


  /* initialise energy correction */
  state_p->observables.kinetic_energy_system_correction   = 0.0e0;
  state_p->observables.potential_energy_system_correction = 0.0e0;
  

  if( !info ){

    /* initialise energies */
    if( ENERGIES_UPDATE( constants, *state_p, *config_def_p ) ) info=1;

  }

  if( !info ){

    /* initialise adiabatic_projection */
    if( ADIABATIC_PROJECTION_UPDATE( constants, *state_p, *config_def_p ) ) info=1;

  }

#ifndef __NO_DIPOLE__

  if( !info ){

    /* compute superpositions */
    if( COMPUTE_SUPERPOSITIONS( constants, *state_p, *config_def_p ) ) info=1;

  }

  if( !info ){

    /* compute classical dipole */
    if( COMPUTE_CLASSICAL_DIPOLE( constants, *state_p, *config_def_p ) ) info=1;

  }

#endif /* __NO_DIPOLE__ */

  if( !info ){

    /* initialise ionic_density */
    if( IONIC_DENSITY_UPDATE( constants, *state_p, *config_def_p ) ) info=1;

  }


  if( !info ){

    /* initialise one_body_electronic_density_matrix */
    if( COMPUTE_ONE_BODY_ELECTRONIC_DENSITY_MATRIX( constants, *state_p, *config_def_p ) ) info=1;

  }


  if( !info ){

    /* initialise single_level_populations */
    if( COMPUTE_SINGLE_LEVEL_POPULATIONS( constants, *state_p, *config_def_p ) ) info=1;

  }


  if( !info ){

    /* initialise one_body_electronic_transition_matrix */
    if( COMPUTE_ONE_BODY_ELECTRONIC_TRANSITION_MATRIX( constants, *state_p, *config_def_p ) ) info=1;

  }


#ifndef __NO_DIPOLE__

  if( !info ){

    /* compute dipole single*/
    if( COMPUTE_DIPOLE_SINGLE( constants, *state_p, *config_def_p ) ) info=1;

  }

#endif /* __NO_DIPOLE__ */

  /* copy def variables to tmp ones */
  if( !info ){

    if( CONFIG_COPY( *config_tmp_p, *config_def_p ) ) info=1;

  }

#ifdef __DEBUG__

  if( !info ){

    if( CONFIG_COMPARE( *config_def_p, *config_tmp_p ) ) info=1;
    
  }
  
  
  /* update tmp values... */
  if( !info ){

    if( MU_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

    if( flag_normal_mode_expansion ){

      if( TRANSFORM_MU( constants, *state_p, *config_tmp_p ) ) info=1;
    
    }
    
    if( DISTANCES_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

    if( INITIALISE_EHRENFEST_FRAME( constants, *state_p, *config_tmp_p ) ) info=1;
    
    if( HAMILTONIAN_MANY_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

    if( N_chain ){

      if( FORCES_THERMOSTAT_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

    }

  }
  
  
  /* ...and check whether they are same  */
  if( !info ){

    if( CONFIG_COMPARE( *config_def_p, *config_tmp_p ) ) info=1;

  }
  
#endif /* __DEBUG__ */


  //BUGFIX: why TMP?
  /* initialise rho_dot */
  if( !info ){

    if( flag_normal_mode_expansion ){

      if( TRANSFORM_HAMILTONIAN( constants, *state_p, *config_tmp_p, config_tmp_p->electrons.delta_F_matrix,  config_tmp_p->electrons.K_matrix ) ) info=1;

    }

    /* initialise rho_dot */
    if( RHO_DOT_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

  }


  if( !info ){

    /* initialise rho_norm_and_projection */
    if( RHO_NORM_AND_PROJECTION_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

  }


  if( !info ){

    /* compute nonadiabaticity */
    // WARNING: it needs the transformed Hamiltonian!
    if( COMPUTE_NONADIABATICITY( constants, *state_p, *config_tmp_p ) ) info=1;

    // if( COMPUTE_NONADIABATIC_RATE( constants, *state_p, *config_tmp_p ) ) info=1;

  }


  return info;

}

//------------------------------------------
