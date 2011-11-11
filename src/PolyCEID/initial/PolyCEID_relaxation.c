
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
#include "PolyCEID_relaxation.h"


#define  TOLL         1000 *( DBL_EPSILON )
#define  MAX_COUNTER  100
#define  N_RATIO      1.0e0
#define  SD_RATIO     0.1e-1


/* private varaiables */

vector  delta_x1;
vector  delta_x2;
vector  F_partial;
vector  F_partial_trans;
matrix  K_partial;
rvector eigenvalues;
matrix  eigenvectors;
matrix  conjugate_eigenvectors;


/*********************
  FUNCTIONS & MACROS
*********************/

//BUGFIX: this was originally TMP
int relaxation( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  /* state*/
  /* dummies */
  unsigned short int SD_flag;
  int     counter=0;
  double  norm;
  int info=0;


#ifdef __DEBUG__

  fprintf( stdout, "DOING: relaxation\n" );

#endif /* __DEBUG__ */


  /* allocation */
  if( ALLOCATE_ONE_STEP_RELAXATION_WORKSPACE( constants ) ) info=1;


  /* compute initial_rho_electron */
  if( COMPUTE_INITIAL_RHO_ELECTRON( constants, *state_p ) )    info=1;
  
  if( DISTANCES_UPDATE( constants, *state_p, *config_p ) )           info=1; /* useless? */

  if( INITIALISE_EHRENFEST_FRAME( constants, *state_p, *config_p ) ) info=1; /* useless? */

  if( F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) )       info=1;

  if( K_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) )       info=1;

  
  if( COMPUTE_F_PARTIAL( constants, *state_p, *config_p, F_partial ) ) info=1;

  /*
    fprintf( stdout, "positions [at the beginning]\n" );
    if( RVECTOR_PRINT_PLUS( stdout, config_p->atoms.positions ) ) info=1;
    fprintf( stdout, "\n" );
  
    fprintf( stdout, "initial_rho_electron [at the beginning]\n" );
    if( MATRIX_PRINT_PLUS( stdout, state_p->initial_rho_electron ) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "Ehrenfest_frame [at the beginning]\n" );
    if( MATRIX_PRINT_PLUS( stdout, config_p->electrons.Ehrenfest_frame ) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "F_partial [at the beginning]\n" );
    if( VECTOR_PRINT_PLUS( stdout, F_partial ) ) info=1;
    fprintf( stdout, "\n" );
  */

  
  if( COMPUTE_K_PARTIAL( constants, *state_p, *config_p, K_partial ) ) info=1;

  /*
    fprintf( stdout, "K_partial\n" );
    if( MATRIX_PRINT_PLUS( stdout, K_partial ) ) info=1;
    fprintf( stdout, "\n" );
  */

  norm = VECTOR_NORM( F_partial );


#ifdef __DEBUG__

  fprintf( stdout, "counter = %d, norm = %le. Next step use 'Newton' [default]\n", counter, norm );

#endif /* __DEBUG__ */


#ifndef __NO_RELAXATION__

  /* relaxation loop */

  if( !info ){

    counter=0;

    SD_flag=0;

    while( norm > TOLL && counter < MAX_COUNTER ){
      
      counter++;


      /* one_step_relaxation */
      if( ONE_STEP_RELAXATION( constants, *state_p, *config_p, SD_flag ) ) info=1;
      
      /*
        fprintf( stdout, "positions [in the middle]\n" );
        if( RVECTOR_PRINT_PLUS( stdout, config_p->atoms.positions ) ) info=1;
        fprintf( stdout, "\n" );
      
	fprintf( stdout, "F_partial [after one-step relaxation]\n" );
	if( VECTOR_PRINT_PLUS( stdout, F_partial ) ) info=1;
	fprintf( stdout, "\n" );
      */
  
      norm     = VECTOR_NORM( F_partial );
      
      // switch
      /*
      if( norm > norm_old ){

	SD_flag=1;

      }
      */

      if( SD_flag ){

#ifdef __DEBUG__

	fprintf( stdout, "counter = %d, norm = %le. Next step use 'steepest descent'\n", counter, norm );

#endif /* __DEBUG__ */
      }
      else{

#ifdef __DEBUG__
	
	fprintf( stdout, "counter = %d, norm = %le. Next step use 'Newton'\n", counter, norm );

#endif /* __DEBUG__ */
      
      }      


      /* breaking instruction */
      if( info ) break;

    } /* end while loop*/

    /*
      fprintf( stdout, "positions [at the end]\n" );
      if( RVECTOR_PRINT_PLUS( stdout, config_p->atoms.positions ) ) info=1;
      fprintf( stdout, "\n" );
    
      fprintf( stdout, "initial_rho_electron [at the end]\n" );
      if( MATRIX_PRINT_PLUS( stdout, state_p->initial_rho_electron ) ) info=1;
      fprintf( stdout, "\n" );

      fprintf( stdout, "Ehrenfest_frame [at the end]\n" );
      if( MATRIX_PRINT_PLUS( stdout, config_p->electrons.Ehrenfest_frame ) ) info=1;
      fprintf( stdout, "\n" );

      fprintf( stdout, "F_partial [at the end]\n" );
      if( VECTOR_PRINT_PLUS( stdout, F_partial ) ) info=1;
      fprintf( stdout, "\n" );

    */

  }

#else /* __NO_RELAXATION__ */

#ifdef __DEBUG__

  fprintf( stdout, "WARNING: relaxation stopped here\n" );

#endif /* __DEBUG__ */

#endif /* __NO_RELAXATION__ */


#ifdef __DEBUG__

  fprintf( stdout, "END: relaxation\n" );

#endif /* __DEBUG__ */


  /* deallocation */
  if( FREE_ONE_STEP_RELAXATION_WORKSPACE() ) info=1;


  return info;

}

//------------------------------------------

//BUGFIX: this was originally TMP
int one_step_relaxation( const constants constants, state_p state_p, config_p config_p, unsigned short int SD_flag ){

  /* constants */
  int N_coor;
  int sdim;
  /* state */
  rvector_p       positions_p;
  /* dummies */
  int i;
  int info=0;


  N_coor          =  constants.N_coor;
  sdim            =  constants.spacial_dimension;
  positions_p     = &(config_p->atoms.positions);


  if( !info ){

    if( !SD_flag ){

      // Newton

      if( DIAGONALISATION( K_partial, eigenvectors, eigenvalues ) ) info=1;
    
      /*
	fprintf( stdout, "eigenvectors\n" );
	if( MATRIX_PRINT_PLUS( stdout, eigenvectors ) ) info=1;
	fprintf( stdout, "\n" );
	
	fprintf( stdout, "eigenvalues\n" );
	if( RVECTOR_PRINT_PLUS( stdout, eigenvalues ) ) info=1;
	fprintf( stdout, "\n" );
      */

      if( MATRIX_ADJOINT( eigenvectors, conjugate_eigenvectors ) ) info=1;

      /*
	fprintf( stdout, "conjugate_eigenvectors\n" );
	if( MATRIX_PRINT_PLUS( stdout, conjugate_eigenvectors ) ) info=1;
	fprintf( stdout, "\n" );
      */

      if( MATRIX_VECTOR_PRODUCT( conjugate_eigenvectors, F_partial, F_partial_trans ) ) info=1;

      /*
	fprintf( stdout, "F_partial_trans\n" );
	if( VECTOR_PRINT_PLUS( stdout, F_partial_trans ) ) info=1;
	fprintf( stdout, "\n" );
      
	fprintf( stdout, "delta_x1 [1]\n" );
	if( VECTOR_PRINT_PLUS( stdout, delta_x1 ) ) info=1;
	fprintf( stdout, "\n" );
      */

      if( !info ){

	for( i=0; i<N_coor; i++ ){
      
	  if( fabs( eigenvalues.rvector[ i ] ) > EPS ){

	    delta_x1.vector[ i ] = CMPLX_DIVISION( F_partial_trans.vector[ i ], CMPLX( fabs( eigenvalues.rvector[ i ] ) ) ); //WARNING: notice the absolute value!

	  }
	  else{

	    delta_x1.vector[ i ] = CMPLX_ZERO;

	  }
	  
	}

      }

      /*
	fprintf( stdout, "delta_x1 [2]\n" );
	if( VECTOR_PRINT_PLUS( stdout, delta_x1 ) ) info=1;
	fprintf( stdout, "\n" );
      */
      
      if( MATRIX_VECTOR_PRODUCT( eigenvectors, delta_x1, delta_x2 ) ) info=1;
      
      /*
	fprintf( stdout, "delta_x2 [2]\n" );
	if( VECTOR_PRINT_PLUS( stdout, delta_x2 ) ) info=1;
	fprintf( stdout, "\n" );
	
	fprintf( stdout, "positions [1]\n" );
	if( RVECTOR_PRINT_PLUS( stdout, *positions_p ) ) info=1;
	fprintf( stdout, "\n" );
      */


      for( i=0; i<N_coor; i++ ){
	
	positions_p->rvector[ i ] += N_RATIO *REAL( delta_x2.vector[ i ] );

	if( constants.flag_periodic_boundary_condition ){

	  positions_p->rvector[ i ] = PBC_DIFF( positions_p->rvector[ i ], constants.cell_dim.rvector[ i %sdim ] );

	}	
      
      }
      
    }
    else{

      // Steepest descent
    
      //      fprintf( stdout, " DOING: steepest descent\n" );
      
      for( i=0; i<N_coor; i++ ){
	
	positions_p->rvector[ i ] += SD_RATIO *REAL( F_partial.vector[ i ] );
      
	if( constants.flag_periodic_boundary_condition ){

	  positions_p->rvector[ i ] = PBC_DIFF( positions_p->rvector[ i ], constants.cell_dim.rvector[ i %sdim ] );

	}	
      
      }

    }

    /*
      fprintf( stdout, "positions [2]\n" );
      if( RVECTOR_PRINT_PLUS( stdout, *positions_p ) ) info=1;
      fprintf( stdout, "\n" );
    */

    /* compute initial_rho_electron */
    if( COMPUTE_INITIAL_RHO_ELECTRON( constants, *state_p ) ) info=1;
  
    if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;

    if( INITIALISE_EHRENFEST_FRAME( constants, *state_p, *config_p ) ) info=1;
    
    if( F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;
    
#ifdef __K_MATRIX_UPDATE__

    if( K_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;

#endif /* __K_MATRIX_UPDATE__ */ 

    if( COMPUTE_F_PARTIAL( constants, *state_p, *config_p, F_partial ) ) info=1;

    /*
      fprintf( stdout, "F_partial\n" );
      if( VECTOR_PRINT_PLUS( stdout, F_partial ) ) info=1;
      fprintf( stdout, "\n" );
    */
    
    if( COMPUTE_K_PARTIAL( constants, *state_p, *config_p, K_partial ) ) info=1;
    
    /*
      fprintf( stdout, "K_partial\n" );
      if( MATRIX_PRINT_PLUS( stdout, K_partial ) ) info=1;
      fprintf( stdout, "\n" );
    */

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "\n" );
  fprintf( stdout, "RELAXED positions \n" );
  if( RVECTOR_PRINT_PLUS( stdout, *positions_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int allocate_one_step_relaxation_workspace( const constants constants ){

  /* constants */
  int N_coor;
  /* dummies */
  int info=0;


  N_coor = constants.N_coor;

  
  if( !info ){

    /* allocations */
    if( ALLOCATE_COMPUTE_INITIAL_RHO_ELECTRON_WORKSPACE( constants ) ) info=1;

    if( VECTOR_ALLOCATE( N_coor, delta_x1 ) ) info=1;
    
    if( VECTOR_ALLOCATE( N_coor, delta_x2 ) ) info=1;
    
    if( VECTOR_ALLOCATE( N_coor, F_partial ) ) info=1;
    
    if( VECTOR_ALLOCATE( N_coor, F_partial_trans ) ) info=1;
    
    if( MATRIX_ALLOCATE( N_coor, N_coor, K_partial ) ) info=1;
    
    if( RVECTOR_ALLOCATE( N_coor, eigenvalues ) ) info=1;
    
    if( MATRIX_ALLOCATE( N_coor, N_coor, eigenvectors ) ) info=1;
    
    if( MATRIX_ALLOCATE( N_coor, N_coor, conjugate_eigenvectors ) ) info=1;

  }


  return info;

}

//------------------------------------------

int free_one_step_relaxation_workspace( void ){

  /* dummies */
  int info=0;


  if( !info ){

    /* deallocations */
    if( VECTOR_FREE( delta_x1 ) ) info=1;
    
    if( VECTOR_FREE( delta_x2 ) ) info=1;
    
    if( VECTOR_FREE( F_partial ) ) info=1;
    
    if( VECTOR_FREE( F_partial_trans ) ) info=1;
    
    if( MATRIX_FREE( K_partial ) ) info=1;
  
    if( RVECTOR_FREE( eigenvalues ) ) info=1;

    if( MATRIX_FREE( eigenvectors ) ) info=1;
    
    if( MATRIX_FREE( conjugate_eigenvectors ) ) info=1;

    if( FREE_COMPUTE_INITIAL_RHO_ELECTRON_WORKSPACE() ) info=1;

  }


  return info;

}

//------------------------------------------
//------------------------------------------

//BUGFIX: this was originally TMP
int compute_F_partial( const constants constants, state_p state_p, config_p config_p, vector_p F_partial_p ){

  /* constants */
  int             N_coor;
  /* state */
  matrix_p        initial_rho_electron_p;
  matrix_array_p  F_matrix_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  int i;
  int info=0;


  N_coor                   =  constants.N_coor;

  initial_rho_electron_p   = &(state_p->initial_rho_electron);
  F_matrix_p               = &(config_p->electrons.F_matrix);
  dummy_matrix1_p          = &state_p->dummy_matrix1;


  for( i=0; i<N_coor; i++ ){

    if( MATRIX_MATRIX_PRODUCT( F_matrix_p->array[ i ], *initial_rho_electron_p, *dummy_matrix1_p ) ) info=1;
  
    F_partial_p->vector[ i ] = MATRIX_TRACE( *dummy_matrix1_p );
    
  } /* i_loop */


  return info;

}

//------------------------------------------

//BUGFIX: this was originally TMP
int compute_K_partial( const constants constants, state_p state_p, config_p config_p, matrix_p K_partial_p ){

  /* constants */
  int             N_coor;
  /* state */
  matrix_p        initial_rho_electron_p;
  matrix_array_p  K_matrix_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  int index;
  int i, j;
  int info=0;


  N_coor                   =  constants.N_coor;

  initial_rho_electron_p   = &(state_p->initial_rho_electron);
  K_matrix_p               = &(config_p->electrons.K_matrix);
  dummy_matrix1_p          = &state_p->dummy_matrix1;


  for( i=0; i<N_coor; i++ ){

    for( j=0; j<N_coor; j++ ){

      index=COORDINATE_INDEX( i, j );

      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ index ], *initial_rho_electron_p, *dummy_matrix1_p ) ) info=1;

      K_partial_p->matrix[ index ] = MATRIX_TRACE( *dummy_matrix1_p );

    } /* j loop */

  } /* i_loop */


#ifndef __NO_HESSIAN_CORRECTIONS__

 /* Hamiltonian */
  if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;

  if( H_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;

  if( F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;


  /* correction to the bare K_partial */
  if( HESSIAN_CORRECTIONS( constants, *state_p, *config_p, *K_partial_p ) ) info=1;

#endif /* __NO_HESSIAN_CORRECTIONS__ */


  return info;

}

//------------------------------------------
