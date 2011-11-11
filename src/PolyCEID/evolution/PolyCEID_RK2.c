
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

#include "config.h"
#include "PolyCEID_RK2.h"



/*********************
  FUNCTIONS & MACROS
*********************/

int RK2_electrons( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p, double ratio ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  //unsigned short int  flag_no_Ehrenfest_frame;
  /* dummies */
  int info=0;

  
  flag_normal_mode_expansion =  constants.flag_normal_mode_expansion;
  //flag_no_Ehrenfest_frame    =  constants.flag_no_Ehrenfest_frame;


  /*-------------------*/
  /* half initial step */
  /*-------------------*/
  if( SINGLE_STEP_ELECTRONS_RK2( constants, *state_p, *config_tmp_p, *config_def_p, 0.5e0 *ratio ) ) info=1;

  //if( !flag_no_Ehrenfest_frame ){

  /* half-step Ehrenfest_frame_update */
    //if( EHRENFEST_FRAME_UPDATE( constants, *state_p, 0.5e0 *ratio ) )  info=1; //WARNING: this changes Ehrenfest_frame_tmp!

  //}


  /* update */
  if( !info ){

    if( MU_UPDATE( constants, *state_p, *config_tmp_p ) )       info=1;

    if( flag_normal_mode_expansion ){
                                                               
      if( TRANSFORM_MU( constants, *state_p, *config_tmp_p ) )  info=1;

    }

    // It's needed because of Delta F
#ifdef __K_MATRIX_UPDATE__ 

    if( HAMILTONIAN_MANY_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1; 

#else /* __K_MATRIX_UPDATE__ */

    if( MATRIX_ARRAY_COPY( config_tmp_p->electrons.K_matrix, config_def_p->electrons.K_matrix ) ) info=1;

    if( HAMILTONIAN_MANY_AUX_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1; 

#endif /* __K_MATRIX_UPDATE__ */


    if( flag_normal_mode_expansion ){

      if( TRANSFORM_HAMILTONIAN( constants, *state_p, *config_tmp_p, config_tmp_p->electrons.delta_F_matrix, config_tmp_p->electrons.K_matrix ) ) info=1;

    }

    if( RHO_DOT_UPDATE( constants, *state_p, *config_tmp_p ) )        info=1;

  }

#ifdef __DEBUG__

  if( RHO_NORM_AND_PROJECTION_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

#endif /* __DEBUG__ */


  /*-----------------------------------------*/
  /* whole step according to the new rho_dot */
  /*-----------------------------------------*/
  if( SINGLE_STEP_ELECTRONS_RK2( constants, *state_p, *config_tmp_p, *config_def_p, ratio ) ) info=1; /* WARNING: this term is different */


  //if( !flag_no_Ehrenfest_frame ){

  /* half-step Ehrenfest_frame_update */
    //if( EHRENFEST_FRAME_UPDATE( constants, *state_p, 0.5e0 *ratio ) )  info=1; //WARNING: this changes Ehrenfest_frame_tmp!

  //}


  return info;

}

//------------------------------------------
//------------------------------------------

/* utilities */

//------------------------------------------

int single_step_electrons_RK2( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p, double ratio ){


  /* constants*/
  int             N_chain;
  int             max_rho_index;
  double          dt;
  /* state */
  matrix_array_p  rho_def_p;
  matrix_array_p  rho_tmp_p;
  matrix_array_p  rho_dot_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  int             i;
  int             info=0;

  
  N_chain         =  constants.N_chain;
  max_rho_index   =  constants.max_rho_index;

  dt              =  ratio *( constants.dt );
  rho_def_p       = &(config_def_p->electrons.rho);
  rho_tmp_p       = &(config_tmp_p->electrons.rho);
  rho_dot_p       = &(state_p->rho_dot);
  dummy_matrix1_p = &state_p->dummy_matrix1;


  if( N_chain ){

#ifdef __TIME_SCALE__

    /* BUGFIX: why tmp? */
    dt *= exp( config_tmp_p->thermostat.positions.rvector[ 0 ] );

#endif /* __TIME_SCALE__ */

    // fprintf( stdout, "time_scale = %le [single_step_electrons_RK2]\n", exp( state_p->thermostat_tmp.positions.rvector[ 0 ] ) );

  }        

  for( i=0; i<max_rho_index; i++ ){

    if( MATRIX_SCALAR_PRODUCT( rho_dot_p->array[ i ], CMPLX( dt ), *dummy_matrix1_p ) ) info=1;

    if( MATRIX_SUM( rho_def_p->array[ i ], *dummy_matrix1_p, rho_tmp_p->array[ i ] ) ) info=1;

  }


  return info;

}

//------------------------------------------
