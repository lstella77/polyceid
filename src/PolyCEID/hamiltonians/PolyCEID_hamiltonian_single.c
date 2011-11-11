
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
#include "PolyCEID_hamiltonian_single.h"


#define EPS0     1.0e-15

#define DELTA1   1.0e-6
#define EPS1     1.0e-7

#define DELTA2   1.0e-4
#define EPS2     1.0e-5


/*********************
  FUNCTIONS & MACROS
*********************/


/* pointers to function */
int (* Hamiltonian_single_parameters_check)( const constants );

int (* H_matrix_single_update_p)( const constants, state_p, config_p, matrix_p );
int (* F_matrix_single_update_p)( const constants, state_p, config_p, int, matrix_p );
int (* K_matrix_single_update_p)( const constants, state_p, config_p, int, int, matrix_p );

int (* classical_dipole_single_update_p)( const constants, state_p, config_p, int, matrix_p );

int (* H_matrix_single_check_p)( const constants, rvector_p, matrix_p );


//------------------------------------------

/* H_matrix_single update */

int H_matrix_single_update( const constants constants, state_p state_p, config_p config_p, matrix_p matrix_p ){

  /* dummies */
  int        info=0;


  if( (* H_matrix_single_update_p)( constants, state_p, config_p, matrix_p ) ) info=1;


  return info;

}

//------------------------------------------

/* F_matrix_single update */

int F_matrix_single_update( const constants constants, state_p state_p, config_p config_p, int i, matrix_p matrix_p ){

  /* config */
  ivector_p  mask_p;
  /* dummies */
  int        info=0;


  mask_p = &config_p->atoms.mask;


  if( mask_p->ivector[ i ] ){

    if( (* F_matrix_single_update_p)( constants, state_p, config_p, i, matrix_p ) ) info=1;

  }  
  else{

    if( MATRIX_ZERO( *matrix_p ) ) info=1;

  }        

  return info;

}

//------------------------------------------

/* K_matrix_single update */

int K_matrix_single_update( const constants constants, state_p state_p, config_p config_p, int i, int j, matrix_p matrix_p){

  /* config */
  ivector_p  mask_p;
  /* dummies */
  int        info=0;


  mask_p = &config_p->atoms.mask;


  if( mask_p->ivector[ i ] && mask_p->ivector[ j ] ){

    if( (* K_matrix_single_update_p)( constants, state_p, config_p, i, j, matrix_p ) ) info=1;

  }
  else{

    if( MATRIX_ZERO( *matrix_p ) ) info=1;

  }        


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int initialise_hamiltonian_single( constants constants ){

  /* dummies */
  int info=0;


  /* which Hamiltonian_single? */

  if( !strcmp( constants.hamiltonian.class, "CHAIN_MOD2" ) ){

    Hamiltonian_single_parameters_check = &Hamiltonian_single_CHAIN_MOD2_parameters_check;

    H_matrix_single_update_p         = &H_matrix_single_CHAIN_MOD2_update;
    F_matrix_single_update_p         = &F_matrix_single_CHAIN_MOD2_update;
    K_matrix_single_update_p         = &K_matrix_single_CHAIN_MOD2_update;

    classical_dipole_single_update_p = &classical_dipole_single_CHAIN_MOD2_update;

    H_matrix_single_check_p          = &H_matrix_single_CHAIN_MOD2_check;

  }
  else if( !strcmp( constants.hamiltonian.class, "SSH" ) ){

    Hamiltonian_single_parameters_check = &Hamiltonian_single_SSH_parameters_check;

    H_matrix_single_update_p         = &H_matrix_single_SSH_update;
    F_matrix_single_update_p         = &F_matrix_single_SSH_update;
    K_matrix_single_update_p         = &K_matrix_single_SSH_update;

    classical_dipole_single_update_p = &classical_dipole_single_SSH_update;

    H_matrix_single_check_p          = &H_matrix_single_SSH_check;

  }
  else{

    fprintf( stderr, "ERROR: this Hamiltonian_single is not supported\n" );
    fflush( stderr );

    info=1;

  }

  if( !info ){

    /* parameters check */
    if( (*Hamiltonian_single_parameters_check)( constants ) ) info=1;

  }


  return info;

}

//------------------------------------------
//------------------------------------------

#ifdef __DEBUG__

//------------------------------------------

/* Hamiltonian_single_check */

int Hamiltonian_single_check( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  /* state */
  rvector_p       positions_p;
  rvector_p       dummy_positions_p;
  ivector_p       mask_p;
  matrix_p        dummy_matrix_single1_p;
  matrix_p        dummy_matrix_single2_p;
  matrix_p        dummy_matrix_single3_p;
  matrix_p        dummy_matrix_single4_p;
  /*dummies*/
  int i_coor1, i_coor2;
  int info=0;


  N_coor                    = constants.N_coor;

  positions_p               = &(config_p->atoms.positions);
  mask_p                    = &(config_p->atoms.mask);
  dummy_positions_p         = &(state_p->dummy_positions);
  dummy_matrix_single1_p    = &(state_p->dummy_matrix_single1);
  dummy_matrix_single2_p    = &(state_p->dummy_matrix_single2);
  dummy_matrix_single3_p    = &(state_p->dummy_matrix_single3);
  dummy_matrix_single4_p    = &(state_p->dummy_matrix_single4);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: Hamiltonian_single_check\n" );

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    /* H_matrix_single check */

    if( RVECTOR_COPY( *dummy_positions_p, *positions_p ) ) info=1;

    if( H_MATRIX_SINGLE_CHECK( constants, *dummy_positions_p, *dummy_matrix_single1_p ) ) info=1;

    if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, *dummy_matrix_single2_p ) ) info=1;

    if( MATRIX_COMPARE_ENTRIES( *dummy_matrix_single1_p, *dummy_matrix_single2_p, EPS0  ) ) info=1;

    if( info ){

      fprintf( stderr, "----------\n" );
      fprintf( stderr, "ERROR in H_matrix_single check [The first one is the numerical estimate]\n" );
      fprintf( stderr, "----------\n" );

    }

  }

  if( !info ){

    /* F_matrix_single check */

    for( i_coor1=0; i_coor1<N_coor; i_coor1++ ){

      if( info ) break;
      
      if( mask_p->ivector[ i_coor1 ] ){ 
      
        //    if( RVECTOR_COPY( *dummy_positions_p, *positions_p ) ) info=1;
      
        // first term;
      
        dummy_positions_p->rvector[ i_coor1 ] += DELTA1;
      
        if( H_MATRIX_SINGLE_CHECK( constants, *dummy_positions_p, *dummy_matrix_single1_p ) ) info=1;
      
#ifdef __DEBUG_PLUS__
      
        fprintf( stdout, "H_matrix_single[1]\n" );
        if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
        fprintf( stdout, "\n" );
	    
#endif /* __DEBUG_PLUS__ */
	    
        dummy_positions_p->rvector[ i_coor1 ] -= DELTA1;

        // second term;

        dummy_positions_p->rvector[ i_coor1 ] -= DELTA1;

        if( H_MATRIX_SINGLE_CHECK( constants, *dummy_positions_p, *dummy_matrix_single2_p ) ) info=1;
      
#ifdef __DEBUG_PLUS__
      
        fprintf( stdout, "H_matrix_single[2]\n" );
        if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
        fprintf( stdout, "\n" );
	    
#endif /* __DEBUG_PLUS__ */
	    
        dummy_positions_p->rvector[ i_coor1 ] += DELTA1;

        // final

        if( MATRIX_DIF( *dummy_matrix_single1_p, *dummy_matrix_single2_p, *dummy_matrix_single1_p ) ) info=1; // WARNING: overwriting
      
        /* WARNING: the F_matrix_single is MINUS the derivative of H_matrix_single */
        if( MATRIX_SCALAR_DIVISION( *dummy_matrix_single1_p, CMPLX( -2.0e0 *DELTA1 ), *dummy_matrix_single1_p ) ) info=1; // WARNING: overwriting


        if( F_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, i_coor1, *dummy_matrix_single2_p ) ) info=1;

        if( MATRIX_COMPARE_ENTRIES( *dummy_matrix_single1_p, *dummy_matrix_single2_p, EPS1  ) ) info=1;

        if( info ){

	  fprintf( stderr, "----------\n" );
	  fprintf( stderr, "ERROR in F_matrix_single check [The first one is the numerical estimate]\n" );
	  fprintf( stderr, "---> coor: (%d)\n", i_coor1 );
	  fprintf( stderr, "----------\n" );

        }

      }  

    } /* i_coor1 loop */
    
  }

  if( !info ){

    /* K_matrix_single check */

    for( i_coor1=0; i_coor1<N_coor; i_coor1++ ){
      
      for( i_coor2=0; i_coor2<N_coor; i_coor2++ ){

	if( info ) break;
	
        if( mask_p->ivector[ i_coor1 ] && mask_p->ivector[ i_coor2 ] ){ 
      
	  //      if( RVECTOR_COPY( *dummy_positions_p, *positions_p ) ) info=1;

	  // first term

	  dummy_positions_p->rvector[ i_coor1 ] += DELTA2;

	  dummy_positions_p->rvector[ i_coor2 ] += DELTA2;

	  if( H_MATRIX_SINGLE_CHECK( constants, *dummy_positions_p, *dummy_matrix_single1_p ) ) info=1;

#ifdef __DEBUG_PLUS__
      
          fprintf( stdout, "H_matrix_single[1]\n" );
          if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
          fprintf( stdout, "\n" );
	    
#endif /* __DEBUG_PLUS__ */
	    
	  dummy_positions_p->rvector[ i_coor2 ] -= DELTA2;
	
	  dummy_positions_p->rvector[ i_coor1 ] -= DELTA2;

	  // second term

	  dummy_positions_p->rvector[ i_coor1 ] += DELTA2;

	  dummy_positions_p->rvector[ i_coor2 ] -= DELTA2;

	  if( H_MATRIX_SINGLE_CHECK( constants, *dummy_positions_p, *dummy_matrix_single2_p ) ) info=1;

#ifdef __DEBUG_PLUS__
      
          fprintf( stdout, "H_matrix_single[2]\n" );
          if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
          fprintf( stdout, "\n" );
	    
#endif /* __DEBUG_PLUS__ */
	    
	  dummy_positions_p->rvector[ i_coor2 ] += DELTA2;

	  dummy_positions_p->rvector[ i_coor1 ] -= DELTA2;

	  // third term

          dummy_positions_p->rvector[ i_coor1 ] -= DELTA2;

	  dummy_positions_p->rvector[ i_coor2 ] += DELTA2;

	  if( H_MATRIX_SINGLE_CHECK( constants, *dummy_positions_p, *dummy_matrix_single3_p ) ) info=1;

#ifdef __DEBUG_PLUS__
      
          fprintf( stdout, "H_matrix_single[3]\n" );
          if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single3_p ) ) info=1;
          fprintf( stdout, "\n" );
	    
#endif /* __DEBUG_PLUS__ */
	    
	  dummy_positions_p->rvector[ i_coor2 ] -= DELTA2;

	  dummy_positions_p->rvector[ i_coor1 ] += DELTA2;

	  // fourth term

	  dummy_positions_p->rvector[ i_coor1 ] -= DELTA2;

	  dummy_positions_p->rvector[ i_coor2 ] -= DELTA2;

	  if( H_MATRIX_SINGLE_CHECK( constants, *dummy_positions_p, *dummy_matrix_single4_p ) ) info=1;

#ifdef __DEBUG_PLUS__
      
          fprintf( stdout, "H_matrix_single[4]\n" );
          if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single4_p ) ) info=1;
          fprintf( stdout, "\n" );
	    
#endif /* __DEBUG_PLUS__ */
	    
	  dummy_positions_p->rvector[ i_coor2 ] += DELTA2;

	  dummy_positions_p->rvector[ i_coor1 ] += DELTA2;

	  // final

	  if( MATRIX_DIF( *dummy_matrix_single1_p, *dummy_matrix_single2_p, *dummy_matrix_single1_p ) ) info=1; // WARNING: overwriting

	  if( MATRIX_DIF( *dummy_matrix_single1_p, *dummy_matrix_single3_p, *dummy_matrix_single1_p ) ) info=1; // WARNING: overwriting

	  if( MATRIX_SUM( *dummy_matrix_single1_p, *dummy_matrix_single4_p, *dummy_matrix_single1_p ) ) info=1; // WARNING: overwriting

	  if( MATRIX_SCALAR_DIVISION( *dummy_matrix_single1_p, CMPLX( 4.0e0 *DELTA2 *DELTA2 ), *dummy_matrix_single1_p ) ) info=1; // WARNING: overwriting

	  if( K_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, i_coor1, i_coor2, *dummy_matrix_single2_p ) ) info=1;

	  if( MATRIX_COMPARE_ENTRIES( *dummy_matrix_single1_p, *dummy_matrix_single2_p, EPS2  ) ) info=1;

	  if( info ){

	    fprintf( stderr, "----------\n" );
	    fprintf( stderr, "ERROR in K_matrix_single check [The first one is the numerical estimate]\n" );
	    fprintf( stderr, "---> coor: (%d, %d)\n", i_coor1, i_coor2 );
	    fprintf( stderr, "----------\n" );
	  
	  }

	}  

      } /* i_coor2 loop */

    } /* i_coor1 loop */

  }


  if( info ){

    fprintf( stdout, "WARNING: Hamiltonian_single after a consistency error occurred\n\n" );

    if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, *dummy_matrix_single2_p ) ) info=1;

    fprintf( stdout, "H_matrix_single\n" );
    if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
    fprintf( stdout, "\n" );


    for( i_coor1=0; i_coor1<N_coor; i_coor1++ ){
     
      if( F_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, i_coor1, *dummy_matrix_single2_p ) ) info=1;
 
      fprintf( stdout, "F_matrix_single[%d]\n", i_coor1 );
      if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
      fprintf( stdout, "\n" );

    }


    for( i_coor1=0; i_coor1<N_coor; i_coor1++ ){
    
      for( i_coor2=0; i_coor2<N_coor; i_coor2++ ){
 
	if( K_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, i_coor1, i_coor2, *dummy_matrix_single2_p ) ) info=1;
 
	fprintf( stdout, "K_matrix_single[%d,%d]\n", i_coor1,i_coor2 );
	if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
	fprintf( stdout, "\n" );

      }

    }

    fflush( stdout );

  }


  return info;

}

//------------------------------------------
//------------------------------------------

int matrix_compare_entries( matrix_p mat1_p, matrix_p mat2_p, double eps ){

  /* dummies */
  int     i, j;
  int     dim_row, dim_column;
  int     index;
  double  norm;
  double  norm2;
  complex dummy;
  complex dummy2;
  int     info=0;


  if( MATRIX_CONSINSTENCY( *mat1_p, *mat2_p ) ) info=1;

  if( !info ){

    dim_row    = mat1_p->matrix_dim_row;

    dim_column = mat1_p->matrix_dim_column;

    for( i=0; i<dim_row; i++ ){

      for( j=0; j<dim_column; j++ ){

	if( info ) break;

	index = MATRIX_ARG( i, j, dim_row, dim_column );

	dummy = CMPLX_DIF( mat1_p->matrix[ index ], mat2_p->matrix[ index ] );

	dummy2 = CMPLX_SUM( mat1_p->matrix[ index ], mat2_p->matrix[ index ] );

	norm  = CMPLX_NORM( dummy );

	norm2 = 0.5e0 *CMPLX_NORM( dummy2 );

	if( norm > eps *norm2 && norm > eps ){

	  fprintf( stderr, "ERROR in entry ( %d, %d ): ( %20.15le, %20.15le ) and ( %20.15le, %20.15le ) are different\n", i, j, 
		   mat1_p->matrix[ index ].z[0], mat1_p->matrix[ index ].z[1], mat2_p->matrix[ index ].z[0], mat2_p->matrix[ index ].z[1] );
	  fprintf( stderr, "[ norm  is %12.5le, threshold are %12.5le and %12.5le ]\n", norm, eps, eps *norm2 );
	  fprintf( stderr, "[ norm2 is %12.5le ]\n", norm2 );

	  fflush( stderr );

	  info=1;

	}

      } /* j loop */

    } /* i loop */


  }


  return info;

}

//------------------------------------------

#endif /* __DEBUG__ */

//------------------------------------------

/* H_matrix_single_check */

int H_matrix_single_check( const constants constants, rvector_p positions_p, matrix_p matrix_p ){

  /* dummies */
  int info=0;


  if( (* H_matrix_single_check_p)( constants, positions_p, matrix_p ) ) info=1;


  return info;

}

//------------------------------------------
//------------------------------------------

/* classical_dipole_single update */

int classical_dipole_single_update( const constants constants, state_p state_p, config_p config_p, int comp, matrix_p matrix_p ){

  /* dummies */
  int        info=0;


  if( (* classical_dipole_single_update_p)( constants, state_p, config_p, comp, matrix_p ) ) info=1;


  return info;

}

//------------------------------------------
