
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
#include "PolyCEID_RK4.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int RK4_electrons( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p, double ratio ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  int            N_chain;
  int            max_rho_index;
  int            N_levels_many;
  double         dt;
  /* state */
  matrix_array_p rho_tmp_p;
  matrix_array_p rho_dot_p;
  matrix_p       dummy_matrix2_p;
  /* dummies */
  int            i;
  matrix_array   rho_tmp2;
  int info=0;


  flag_normal_mode_expansion =  constants.flag_normal_mode_expansion;


  N_chain         =  constants.N_chain;
  max_rho_index   =  constants.max_rho_index;
  N_levels_many   =  constants.N_levels_many;
  dt              =  ratio *( constants.dt );

  rho_tmp_p       = &(config_tmp_p->electrons.rho);
  rho_dot_p       = &(state_p->rho_dot);
  dummy_matrix2_p = &(state_p->dummy_matrix2);


  if( N_chain ){

#ifdef __TIME_SCALE__

    /* BUGFIX: why tmp? */
    dt *= exp( config_tmp_p->thermostat.positions.rvector[ 0 ] );

#endif /* __TIME_SCALE__ */

  }
		    
  
  /* allocation */
  if( MATRIX_ARRAY_ALLOCATE( max_rho_index+1, N_levels_many, N_levels_many, rho_tmp2 ) ) info=1;

  /* initial copy */
  if( MATRIX_ARRAY_COPY( rho_tmp2, *rho_tmp_p) ) info=1;


  /*-------------------*/
  /* partial update */
  for( i=0; i<max_rho_index; i++ ){

    if( MATRIX_SCALAR_PRODUCT( rho_dot_p->array[ i ], CMPLX( ONEO6 *dt ), *dummy_matrix2_p ) ) info=1;
    if( MATRIX_SUM(  rho_tmp2.array[ i ],  *dummy_matrix2_p, rho_tmp2.array[ i ] ) ) info=1; /* WARNING: overwriting */

  }


  /* first step */
  if( SINGLE_STEP_ELECTRONS_RK4( constants, *state_p, *config_tmp_p, *config_def_p, 0.5e0 *ratio ) ) info=1;
  /*-------------------*/


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING 1st update\n" );

#endif /* __DEBUG_PLUS__ */


  /* update */
  if( !info ){

    if( MU_UPDATE( constants, *state_p, *config_tmp_p ) )             info=1;

    if( flag_normal_mode_expansion ){

      if( TRANSFORM_MU( constants, *state_p, *config_tmp_p ) )          info=1;

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


  /*-------------------*/
  /* partial update */
  for( i=0; i<max_rho_index; i++ ){

    if( MATRIX_SCALAR_PRODUCT( rho_dot_p->array[ i ], CMPLX( ONEO3 *dt ), *dummy_matrix2_p ) ) info=1;
    if( MATRIX_SUM( rho_tmp2.array[ i ],  *dummy_matrix2_p, rho_tmp2.array[ i ] ) ) info=1; /* WARNING: overwriting */

  }

  /* second step */
  if( SINGLE_STEP_ELECTRONS_RK4( constants, *state_p, *config_tmp_p, *config_def_p, 1.0e0 *ratio ) ) info=1;
  /*-------------------*/


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING 2nd update\n" );

#endif /* __DEBUG_PLUS__ */


  /* update */
  if( !info ){

    if( MU_UPDATE( constants, *state_p, *config_tmp_p ) )             info=1;

    if( flag_normal_mode_expansion ){

      if( TRANSFORM_MU( constants, *state_p, *config_tmp_p ) )          info=1;

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

  /*-------------------*/
  /* partial update */
  for( i=0; i<max_rho_index; i++ ){

    if( MATRIX_SCALAR_PRODUCT( rho_dot_p->array[ i ], CMPLX( ONEO3 *dt ), *dummy_matrix2_p ) ) info=1;
    if( MATRIX_SUM( rho_tmp2.array[ i ],  *dummy_matrix2_p, rho_tmp2.array[ i ] ) ) info=1; /* WARNING: overwriting */

  }

  /* third step */
  if( SINGLE_STEP_ELECTRONS_RK4( constants, *state_p, *config_tmp_p, *config_def_p, 0.5e0 *ratio ) ) info=1;
  /*-------------------*/


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING 3rd update\n" );

#endif /* __DEBUG_PLUS__ */


  /* update */
  if( !info ){

    if( MU_UPDATE( constants, *state_p, *config_tmp_p ) )             info=1;

    if( flag_normal_mode_expansion ){

      if( TRANSFORM_MU( constants, *state_p, *config_tmp_p ) )          info=1;

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

  /*-------------------*/
  /* partial update */
  for( i=0; i<max_rho_index; i++ ){

    if( MATRIX_SCALAR_PRODUCT( rho_dot_p->array[ i ], CMPLX( ONEO6 *dt ), *dummy_matrix2_p ) ) info=1;
    if( MATRIX_SUM( rho_tmp2.array[ i ],  *dummy_matrix2_p, rho_tmp2.array[ i ] ) ) info=1; /* WARNING: overwriting */

  }

  /* no update is needed */
  /*-------------------*/

  /* final copy */
  if( MATRIX_ARRAY_COPY( *rho_tmp_p, rho_tmp2 ) ) info=1;


  /* deallocation */
  if( MATRIX_ARRAY_FREE( rho_tmp2 ) ) info=1;


  return info;

}

//------------------------------------------
//------------------------------------------

/* utilities */

//------------------------------------------

int single_step_electrons_RK4( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p, double ratio ){


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

  }


  for( i=0; i<max_rho_index; i++ ){

    if( MATRIX_SCALAR_PRODUCT( rho_dot_p->array[ i ], CMPLX( dt ), *dummy_matrix1_p ) ) info=1;

    if( MATRIX_SUM( rho_def_p->array[ i ], *dummy_matrix1_p, rho_tmp_p->array[ i ] ) ) info=1;

  }


  return info;

}

//------------------------------------------
