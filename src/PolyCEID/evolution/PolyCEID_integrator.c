
/******************************************************************************

  Copyright (C) 2011 by Lorenzo Stella <lorenzo DOT stella77 AT gmail DOT com>

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
#include "PolyCEID_integrator.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int global_CEID_integrator( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p ){

  /* constants */
  int                 skip_write;
  double              dt;
  double              time_length;
  /* state */
  long unsigned int*  step_counter_p;
  double*             time_p;
  /* dummies */
  int  info=0;


  skip_write            =  constants.skip_write;
  dt                    =  constants.dt;
  time_length           =  constants.time_length;

  step_counter_p        = &(state_p->step_counter);
  time_p                = &(state_p->time);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: integrator\n" );

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: observables [start]\n" );

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    /* compute observables */
    if( COMPUTE_OBSERVABLES( constants, *state_p, *config_def_p ) ) info=1;
	  
    /* print observables */
    if( PRINT_OBSERVABLES( constants, *state_p, *config_def_p ) ) info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: initial copy\n" );

#endif /* __DEBUG_PLUS__ */


  /* initial copy */
  // BUGFIX: in principle they should equal after the initialisation...
  if( CONFIG_COPY( *config_tmp_p, *config_def_p ) ) info=1;


  /*--------------------------------
   -------MAIN INTEGRATOR LOOP------
   ---------------------------------*/
  while( ( time_length - (*time_p) ) > EPS ){


#ifdef __DEBUG__

    if( CONFIG_COMPARE( *config_def_p, *config_tmp_p ) ) info=1;
	
#endif /* __DEBUG__ */


    /* step_counter increase */
    (*step_counter_p)++;
      
      
#ifdef __DEBUG_PLUS__

    fprintf( stdout, "# ---> integration step %lu.\n", *step_counter_p );
    fprintf( stdout, "# ---> started at time = %e ...\n", *time_p );
    fflush( stdout );
      
#endif /* __DEBUG_PLUS__ */
      
      
#ifdef __DEBUG_PLUS__

    fprintf( stdout, "DOING: single step integrator\n" );

#endif /* __DEBUG_PLUS__ */


    /* single step integrator */
    if( SINGLE_STEP_CEID_INTEGRATOR( constants, *state_p, *config_tmp_p, *config_def_p ) ) info=1;
	
      
#ifdef __DEBUG_PLUS__

    fprintf( stdout, "DOING: final copy\n" );

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG__

    if( CONFIG_COMPARE( *config_def_p, *config_tmp_p ) ) info=1;
	
#endif /* __DEBUG__ */


    /* info check */
    if( info ){

      fprintf( stderr, "ERROR: Main integrations loop broken\n" );
      fflush( stderr );

      break;

    }      

      
    /* time increasing */
    *time_p += dt;
   

    /* energy corrections */
#ifndef __NO_ENERGY_CORRECTION__


    if( COMPUTE_KINETIC_ENERGY_CORRECTIONS( constants, *state_p, *config_def_p ) ) info=1;

    if( COMPUTE_POTENTIAL_ENERGY_CORRECTIONS( constants, *state_p, *config_def_p ) ) info=1;


#endif /* __NO_ENERGY_CORRECTION__ */


    /* write observables */
    if( 0 == (*step_counter_p) %skip_write ){

	  
#ifdef __DEBUG_PLUS__
	  
      fprintf( stdout, "DOING: observables\n" );

#endif /* __DEBUG_PLUS__ */


      /* compute observables */
      if( COMPUTE_OBSERVABLES( constants, *state_p, *config_def_p ) ) info=1;
	  
      /* print observables */
      if( PRINT_OBSERVABLES( constants, *state_p, *config_def_p ) ) info=1;


#ifdef __DEBUG__

      if( CONFIG_COMPARE( *config_def_p, *config_tmp_p ) ) info=1;
	
#endif /* __DEBUG__ */


    } /* end step_counter conditional */

      
#ifdef __DEBUG_PLUS__

    fprintf( stdout, "# ---> ... concluded at time = %e.\n\n", *time_p );
    fflush( stdout );

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG__

    if( 0 == ( (*step_counter_p) %skip_write) ){

      fprintf( stderr, "+" );
      fflush( stderr );

      }

#endif /* __DEBUG__ */


  } /* end while */


  /* final write observables */
  if( 0 != (*step_counter_p) %skip_write ){
      

#ifdef __DEBUG_PLUS__
      
    fprintf( stdout, "DOING: observables [end]\n" );

#endif /* __DEBUG_PLUS__ */


    /* compute observables */
    if( COMPUTE_OBSERVABLES( constants, *state_p, *config_def_p ) ) info=1;
	  
    /* print observables */
    if( PRINT_OBSERVABLES( constants, *state_p, *config_def_p ) ) info=1;


#ifdef __DEBUG__

    if( CONFIG_COMPARE( *config_def_p, *config_tmp_p ) ) info=1;
	
#endif /* __DEBUG__ */


  } /* end time_length conditional */
     

#ifdef __DEBUG__

  if( !info ){

    fprintf( stderr, ":-)\n\n" );

  }
  else{

    fprintf( stderr, ":-(\n\n" );

  }

  fflush( stderr );

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------

/* utilites */

//------------------------------------------

//BUGFIX: this was originally TMP
int single_step_CEID_integrator( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  unsigned short int  flag_no_Ehrenfest_frame;
  /* dummies */
  int  info=0;
  
  
  flag_normal_mode_expansion =  constants.flag_normal_mode_expansion;
  flag_no_Ehrenfest_frame    =  constants.flag_no_Ehrenfest_frame;


  //WARNING: no update a this stage: there have been a copy before

#ifdef __DEBUG_PLUS__

#ifndef __NO_VERLET__

  fprintf( stdout, "   DOING: verlet [1]\n" );

#endif /* __NO_VERLET__ */
    
#endif /* __DEBUG_PLUS__ */
  

#ifndef __NO_VERLET__

  /* first half step ions */
  if( VELOCITY_VERLET_IONS( constants, *state_p, *config_tmp_p, *config_def_p, 0.5e0 ) ) info=1;
    
#endif /* __NO_VERLET__ */


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "   DOING: update\n" );

#endif /* __DEBUG_PLUS__ */


  /* global update */
  if( GLOBAL_UPDATE( constants, *state_p, *config_def_p, *config_tmp_p ) ) info=1;

  /* copy back after a global update */
  if( CONFIG_COPY( *config_def_p, *config_tmp_p ) ) info=1;

  if( flag_normal_mode_expansion ){

    if( TRANSFORM_HAMILTONIAN( constants, *state_p, *config_tmp_p, config_tmp_p->electrons.delta_F_matrix, config_tmp_p->electrons.K_matrix ) ) info=1;

  }

  /* rho_dot_update */
  if( RHO_DOT_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;


  if( !flag_no_Ehrenfest_frame ){

    /* half-step Ehrenfest_frame_update */
    if( EHRENFEST_FRAME_UPDATE( constants, *state_p, *config_def_p, *config_tmp_p, 0.5e0 ) )  info=1;

  }


  /* the following part is useless in a pure Ehrenfest simulation */
#ifndef __NO_CEID__


  /* compute nonadiabaticity */
  if( 0 == (state_p->step_counter) %(constants.skip_write) ){

    // WARNING! It needs the transformed Hamiltonian!
    if( COMPUTE_NONADIABATICITY( constants, *state_p, *config_tmp_p ) ) info=1;
      
    // if( COMPUTE_NONADIABATIC_RATE( constants, *state_p, *config_tmp_p ) ) info=1;
      
  }


#ifndef __DEBUG__ 

  if( 0 == ( state_p->step_counter ) %( constants.skip_write ) ){

    if( RHO_NORM_AND_PROJECTION_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

  }

#else /* __DEBUG__ */

  if( RHO_NORM_AND_PROJECTION_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

#endif /* __DEBUG__ */


#ifdef __DEBUG_PLUS__

#ifndef __RK4__

  fprintf( stdout, "   DOING: RK2\n" );

#else /* __RK4__ */

  fprintf( stdout, "   DOING: RK4\n" );

#endif /* __RK4__ */

#endif /* __DEBUG_PLUS__ */


  /* single step electrons */
#ifndef __RK4__

  if( RK2_ELECTRONS( constants, *state_p, *config_tmp_p, *config_def_p, 1.0e0 ) ) info=1;

#else /* __RK4__ */

  if( RK4_ELECTRONS( constants, *state_p, *config_tmp_p, *config_def_p, 1.0e0 ) ) info=1;

#endif /* __RK4__ */

#endif /* __NO_CEID__ */


  if( !flag_no_Ehrenfest_frame ){

    /* half-step Ehrenfest_frame_update */
    if( EHRENFEST_FRAME_UPDATE( constants, *state_p, *config_def_p, *config_tmp_p, 0.5e0 ) )  info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "   DOING: update\n" );

#endif /* __DEBUG_PLUS__ */


  /* global update */
  if( GLOBAL_UPDATE( constants, *state_p, *config_def_p, *config_tmp_p ) ) info=1;

  /* copy back after a global update */
  if( CONFIG_COPY( *config_def_p, *config_tmp_p ) ) info=1;


#ifdef __DEBUG_PLUS__
    
#ifndef __NO_VERLET__

  fprintf( stdout, "   DOING: verlet [2]\n" );
    
#endif /* __NO_VERLET__ */

#endif /* __DEBUG_PLUS__ */


#ifndef __NO_VERLET__

  /* second half step ions */
  if( VELOCITY_VERLET_IONS( constants, *state_p, *config_tmp_p, *config_def_p, 0.5e0 ) ) info=1;
    
#endif /* __NO_VERLET__ */

  /* global update */
  if( GLOBAL_UPDATE( constants, *state_p, *config_def_p, *config_tmp_p ) ) info=1;

  /* copy backi after a global update */
  if( CONFIG_COPY( *config_def_p, *config_tmp_p ) ) info=1;


  return info;

}

//------------------------------------------
// WARNING: this function was originally DEF
int compute_observables( const constants constants, state_p state_p, config_p config_p ){

  /* dummy */
  int info=0;


  /* update energy */
  if( ENERGIES_UPDATE( constants, *state_p, *config_p ) ) info=1;

  /* compute adiabatic projection */
  if( ADIABATIC_PROJECTION_UPDATE( constants, *state_p, *config_p ) ) info=1;

#ifndef __NO_DIPOLE__

  /* compute superpositions */
  if( COMPUTE_SUPERPOSITIONS( constants, *state_p, *config_p ) ) info=1;

  /* compute classical dipole */
  if( COMPUTE_CLASSICAL_DIPOLE( constants, *state_p, *config_p ) ) info=1;

#endif /* __NO_DIPOLE__ */

  /* update one_body_electronic_density_matrix */
  if( COMPUTE_ONE_BODY_ELECTRONIC_DENSITY_MATRIX( constants, *state_p, *config_p ) ) info=1; 

  /* compute single level populations */
  if( COMPUTE_SINGLE_LEVEL_POPULATIONS( constants, *state_p, *config_p ) ) info=1;
	  
  /* compute ionic density */
  if( IONIC_DENSITY_UPDATE( constants, *state_p, *config_p ) ) info=1;
	  
  /* update one_body_electronic_transition_matrix */
  if( COMPUTE_ONE_BODY_ELECTRONIC_TRANSITION_MATRIX( constants, *state_p, *config_p ) ) info=1;

#ifndef __NO_DIPOLE__

  /* compute dipole single*/
  if( COMPUTE_DIPOLE_SINGLE( constants, *state_p, *config_p ) ) info=1;

#endif /* __NO_DIPOLE__ */


  return info;

}  
	
//------------------------------------------

int global_update( const constants constants, state_p state_p, config_p config_def_p, config_p config_tmp_p ){

  /* dummy */
  int info=0;


  if( MU_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

  if( constants.flag_normal_mode_expansion ){

    if( TRANSFORM_MU( constants, *state_p, *config_tmp_p ) ) info=1;

  }           

#ifdef __K_MATRIX_UPDATE__ 

  if( HAMILTONIAN_MANY_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1; 

#else /* __K_MATRIX_UPDATE__ */

  if( MATRIX_ARRAY_COPY( config_tmp_p->electrons.K_matrix, config_def_p->electrons.K_matrix ) ) info=1;

  if( HAMILTONIAN_MANY_AUX_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1; 

#endif /* __K_MATRIX_UPDATE__ */                     

  // BUGFIX: should be here?
  if( constants.N_chain ){

    if( FORCES_THERMOSTAT_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

  }        


  return info;

}

//------------------------------------------

