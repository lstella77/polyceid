
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
#include "PolyCEID_Ehrenfest.h"



/*********************
  FUNCTIONS & MACROS
*********************/

int initialise_Ehrenfest_frame( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        N_levels_single;
  /* state */
  matrix_p   Ehrenfest_frame_p;
  rvector_p  dummy_rvector_single_p;
  matrix_p   dummy_matrix_single1_p;
  /* dummies */
  //  double energy_shift;
  int    i, j;
  int    info=0;


  N_levels_single         =  constants.N_levels_single;

  Ehrenfest_frame_p       = &(config_p->electrons.Ehrenfest_frame);
  dummy_rvector_single_p  = &(state_p->dummy_rvector_single);
  dummy_matrix_single1_p  = &(state_p->dummy_matrix_single1);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: initialise_Ehrenfest_frame\n" );

#endif /* __DEBUG_PLUS__ */


  if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, *dummy_matrix_single1_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "H_matrix_single\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  if( DIAGONALISATION( *dummy_matrix_single1_p, *Ehrenfest_frame_p, *dummy_rvector_single_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "eigenvectors [initial Ehrenfest_frame]\n" );
  if( MATRIX_PRINT_PLUS( stdout, *Ehrenfest_frame_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "eigenvalues\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *dummy_rvector_single_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  for( i=0; i<N_levels_single; i++ ){

    if( Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( 0, i )].z[0] < 0.0e0 ){

      for( j=0; j<N_levels_single; j++ ){

	Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( j, i ) ].z[0] = -Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( j, i ) ].z[0];

	Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( j, i ) ].z[1] = -Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( j, i ) ].z[1];


      } /* end j loop */

    } /* end conditional */

  } /* end i loop */

  /*
  energy_shift = 0.0e0;

  for( i=0; i< N_levels_single; i++ ){

    energy_shift += dummy_rvector_single_p->rvector[ i ];

  }

  energy_shift /= (double) N_levels_single;

  fprintf( stdout, "energy_shift = %le\n", energy_shift );
  */


  return info;

}

//------------------------------------------
//------------------------------------------

int Ehrenfest_frame_update( const constants constants, state_p state_p, config_p config_def_p, config_p config_tmp_p, double ratio ){

  /* dummies */
  int  info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOIND: Ehrenfest_frame_update\n" );

#endif /* __DEBUG_PLUS__ */


  if(  EXACT_EHRENFEST( constants, *state_p, *config_def_p, *config_tmp_p, ratio ) ) info=1;


#ifdef __ORTHONORMALISE__

  if( ORTHONORMALISE_EHRENFEST_FRAME( constants, *state_p, *config_tmp_p ) ) info=1;

#endif /* __ORTHONORMALISE__ */


  return info;

}

//------------------------------------------
//------------------------------------------

/* utilities */

//------------------------------------------

int exact_Ehrenfest( const constants constants, state_p state_p, config_p config_def_p, config_p config_tmp_p, double ratio ){

  /* constants*/
  int       N_chain;
  int       N_levels_single;
  double    dt;
  /* state */
  matrix_p  Ehrenfest_frame_p;
  rvector_p dummy_rvector_single_p;
  matrix_p  dummy_matrix_single1_p;
  matrix_p  dummy_matrix_single2_p;
  matrix_p  dummy_matrix_single3_p;
  matrix_p  dummy_matrix_single4_p;
  matrix_p  dummy_Ehrenfest_frame_p;
  /* dummies */
  double    energy_shift;
  complex   shift;
  int       i;
  int       info=0;


  N_chain                 =  constants.N_chain;
  N_levels_single         =  constants.N_levels_single;
  dt                      =  ratio *( constants.dt ); // WARNING: rescaling!
  Ehrenfest_frame_p       = &(config_tmp_p->electrons.Ehrenfest_frame);
  dummy_rvector_single_p  = &state_p->dummy_rvector_single;
  dummy_matrix_single1_p  = &state_p->dummy_matrix_single1;
  dummy_matrix_single2_p  = &state_p->dummy_matrix_single2;
  dummy_matrix_single3_p  = &state_p->dummy_matrix_single3;
  dummy_matrix_single4_p  = &state_p->dummy_matrix_single4;
  dummy_Ehrenfest_frame_p = &(state_p->dummy_Ehrenfest_frame);


  if( N_chain ){

#ifdef __TIME_SCALE__

    dt *= exp( config_def_p->thermostat.positions.rvector[ 0 ] );

#endif /* __TIME_SCALE__ */

    // fprintf( stdout, "time_scale = %le [exact_Eherenfest]\n", exp( config_p->thermostat.positions.rvector[ 0 ] ) ); 

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: exact_Ehrenfest\n" );

#endif /* __DEBUG_PLUS__ */


  /* initial copy */
  if( MATRIX_COPY( *dummy_Ehrenfest_frame_p, *Ehrenfest_frame_p ) ) info=1;


  /*
    for( i=0; i<constants.N_atoms; i++ ){
    
    if( fabs( config_p->atoms.forces.rvector[ i ] +config_p->atoms.forces.rvector[ constants.N_atoms -i -1 ] ) > FLT_EPSILON ){
      
    fprintf( stderr, "ERROR: forces[%d] %le\n", i, config_p->atoms.forces.rvector[ i ] );
    fprintf( stderr, "ERROR: forces[%d] %le\n", constants.N_atoms -i -1, config_p->atoms.forces.rvector[ constants.N_atoms -i -1 ] );
    fprintf( stderr, "\n" );
      
    info=1;
    
    }
    
    }
  */

  /* compute H_matrix_single */
  if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_def_p, *dummy_matrix_single1_p ) ) info=1;

  /*
    for( i=0; i<N_levels_single; i++ ){
    
    for( j=0; j<N_levels_single; j++ ){
    
    if( fabs(dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( i, j ) ].z[0] - 
    dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( N_levels_single -i -1, N_levels_single -j -1  ) ].z[0] ) > FLT_EPSILON ){
	       
    fprintf( stderr, "ERROR: H_matrix_single[%d, %d] %le\n", i, j, dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( i, j ) ].z[0] );
    fprintf( stderr, "ERROR: H_matrix_single[%d, %d] %le\n", i, j, dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( N_levels_single -i -1, N_levels_single -j -1  )].z[0] );
    fprintf( stderr, "\n" );
	       
    info=1;
    
    }
    
    }
    
    }
	       
  */

  /* diagonalisation */
  if( DIAGONALISATION( *dummy_matrix_single1_p, *dummy_matrix_single2_p, *dummy_rvector_single_p  ) ) info=1;

  /*  
    fprintf( stdout, "eigenvectors [H_matrix]\n" );
    if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "eigenvalues [H_matrix]\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *dummy_rvector_single_p ) ) info=1;
    fprintf( stdout, "\n" );
  */


  /* compute energy shift */
  energy_shift = 0.0e0;

  for( i=0; i< N_levels_single; i++ ){

    energy_shift += dummy_rvector_single_p->rvector[ i ];
    
  }

  energy_shift /= (double) N_levels_single;

  //  fprintf( stdout, "energy_shift = %le\n", energy_shift );


#ifdef __DEBUG__

  /* diagonalisation check */
  if( MATRIX_ZERO( *dummy_matrix_single1_p ) ) info=1;

  for( i=0; i<N_levels_single; i++ ){
    
    dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( i, i ) ].z[ 0 ] = dummy_rvector_single_p->rvector[ i ];
    
  }

  if( MATRIX_ADJOINT( *dummy_matrix_single2_p, *dummy_matrix_single3_p ) ) info=1;
  
  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_matrix_single3_p, *dummy_matrix_single4_p ) ) info=1;
  
  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single2_p, *dummy_matrix_single4_p, *dummy_matrix_single1_p ) ) info=1;

  /* compute H_matrix_single */
  if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_def_p, *dummy_matrix_single3_p ) ) info=1;

  /* compare */
  if( MATRIX_COMPARE( *dummy_matrix_single1_p, *dummy_matrix_single3_p ) ){

    fprintf( stderr, "ERROR: inconsistency in H_matrix_single diagonalisation\n" );

    info=1;

  }

#endif /* __DEBUG__ */


  /* eigenvalue rescaling */
  for( i=0; i< N_levels_single; i++ ){

    dummy_rvector_single_p->rvector[ i ] -= energy_shift;
    
  }


  /* unitary transform */
  if( MATRIX_ZERO( *dummy_matrix_single1_p ) ) info=1;

  for( i=0; i<N_levels_single; i++ ){
    
    dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( i, i ) ].z[ 0 ] = cos( PIO4 -( dummy_rvector_single_p->rvector[ i ] ) *dt /HBAR );
    dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( i, i ) ].z[ 1 ] = sin( PIO4-( dummy_rvector_single_p->rvector[ i ] ) *dt /HBAR );

  }

  if( MATRIX_ADJOINT( *dummy_matrix_single2_p, *dummy_matrix_single3_p ) ) info=1;
  
  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_matrix_single3_p, *dummy_matrix_single4_p ) ) info=1;
  
  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single2_p, *dummy_matrix_single4_p, *dummy_matrix_single1_p ) ) info=1;

  shift = CMPLX_INIT( OSQRT2, -OSQRT2 );

  /*
    fprintf( stdout, "complex shift\n" );
    if( CMPLX_PRINT_PLUS( stdout, shift ) ) info=1;
    fprintf( stdout, "\n" );
  */

  if( MATRIX_SCALAR_PRODUCT( *dummy_matrix_single1_p, shift, *dummy_matrix_single1_p) ) info=1;


  /*
    fprintf( stdout, "Ehrenfest evolution\n" );
    if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
    fprintf( stdout, "\n" );
  */


#ifdef __DEBUG__

  /* unitarity check */
  if( MATRIX_ADJOINT( *dummy_matrix_single1_p, *dummy_matrix_single2_p ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_matrix_single2_p, *dummy_matrix_single3_p ) ) info=1;

  if( MATRIX_UNIT( *dummy_matrix_single4_p ) ) info=1;

  if( MATRIX_COMPARE( *dummy_matrix_single3_p, *dummy_matrix_single4_p ) ){

    fprintf( stderr, "ERROR: the Ehrenfest evolution matrix is not unitary\n");

    info=1;

  }

#endif /* __DEBUG__ */

  /*
    for( i=0; i<N_levels_single; i++ ){

    for( j=0; j<N_levels_single; j++ ){
    
    if( fabs(dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( i, j ) ].z[0] - 
    dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( N_levels_single -i -1, N_levels_single -j -1  ) ].z[0] ) > FLT_EPSILON ){
    
    fprintf( stderr, "ERROR: Ehrenfest_propagator[%d, %d] %le\n", i, j, dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( i, j ) ].z[0] );
    fprintf( stderr, "ERROR: Ehrenfest_propagator[%d, %d] %le\n", i, j, dummy_matrix_single1_p->matrix[ ELECTRON_SINGLE_INDEX( N_levels_single -i -1, N_levels_single -j -1  )].z[0] );
    fprintf( stderr, "\n" );
	
    info=1;

    }
	  
    }

    }
  */

  /* evolution */
  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_Ehrenfest_frame_p, *Ehrenfest_frame_p ) ) info=1;
  //  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *Ehrenfest_frame_p, *Ehrenfest_frame_p ) ) info=1;


#ifdef __DEBUG__

  /* unitarity check */
  if( MATRIX_ADJOINT( *Ehrenfest_frame_p, *dummy_matrix_single2_p ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( *Ehrenfest_frame_p, *dummy_matrix_single2_p, *dummy_matrix_single3_p ) ) info=1;

  //  if( MATRIX_UNIT( *dummy_matrix_single4_p ) ) info=1;

  if( MATRIX_COMPARE( *dummy_matrix_single3_p, *dummy_matrix_single4_p ) ){

    fprintf( stderr, "ERROR: Ehrenfest frame is not orthonormal\n");

    info=1;

  }

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------

int orthonormalise_Ehrenfest_frame( const constants constants, state_p state_p, config_p config_p ){

  /* constants*/
  int       N_levels_single;
  /* state */
  //  matrix_p  Ehrenfest_frame_p;
  matrix_p  Ehrenfest_frame_p;
  rvector_p dummy_rvector_single_p;
  matrix_p  dummy_matrix_single1_p;
  matrix_p  dummy_matrix_single2_p;
  matrix_p  dummy_matrix_single3_p;
  matrix_p  dummy_matrix_single4_p;
  matrix_p  dummy_Ehrenfest_frame_p;
  /* dummies */
  int       i;
  double    dummy;
  int       info=0;


  N_levels_single         =  constants.N_levels_single;
  Ehrenfest_frame_p       = &(config_p->electrons.Ehrenfest_frame);
  dummy_rvector_single_p  = &state_p->dummy_rvector_single;
  dummy_matrix_single1_p  = &state_p->dummy_matrix_single1;
  dummy_matrix_single2_p  = &state_p->dummy_matrix_single2;
  dummy_matrix_single3_p  = &state_p->dummy_matrix_single3;
  dummy_matrix_single4_p  = &state_p->dummy_matrix_single4;
  dummy_Ehrenfest_frame_p = &(state_p->dummy_Ehrenfest_frame);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: orthonormalise_Ehrenfest\n" );


  fprintf( stdout, "Ehrenfest_frame [before]\n" );
  if( MATRIX_PRINT_PLUS( stdout, *Ehrenfest_frame_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* initial copy */
  if( MATRIX_COPY( *dummy_Ehrenfest_frame_p, *Ehrenfest_frame_p ) ) info=1;


  if( MATRIX_ADJOINT( *dummy_Ehrenfest_frame_p, *dummy_matrix_single1_p ) ) info=1;


  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_Ehrenfest_frame_p, *dummy_matrix_single2_p ) ) info=1; // S matrix

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "S  matrix\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* diagonalisation */
  if( DIAGONALISATION( *dummy_matrix_single2_p, *dummy_matrix_single1_p, *dummy_rvector_single_p  ) ) info=1;


  /* orthonormalisation transform */
  if( MATRIX_ZERO( *dummy_matrix_single2_p ) ) info=1;

  for( i=0; i<N_levels_single; i++ ){
    
    dummy = dummy_rvector_single_p->rvector[ i ];

    if( dummy < EPS*EPS ){

      fprintf( stderr, "ERROR: diagonal entries of the S matrix must be positive\n" );

      info=1;

    }


    dummy_matrix_single2_p->matrix[ ELECTRON_SINGLE_INDEX( i, i ) ] = CMPLX( 1.0e0/sqrt( dummy ) );

  }

  if( MATRIX_ADJOINT( *dummy_matrix_single1_p, *dummy_matrix_single3_p ) ) info=1;
  
  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single2_p, *dummy_matrix_single3_p, *dummy_matrix_single4_p ) ) info=1;
  
  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_matrix_single4_p, *dummy_matrix_single2_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "orthonormalisation matrix\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* orthonormalisation */
  if( MATRIX_MATRIX_PRODUCT( *dummy_Ehrenfest_frame_p, *dummy_matrix_single2_p, *Ehrenfest_frame_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "Ehrenfest_frame [after]\n" );
  if( MATRIX_PRINT_PLUS( stdout, *Ehrenfest_frame_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------
