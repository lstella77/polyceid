
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
#include "PolyCEID_finite_rotations.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* check_finite_rotation_hamiltonian */
//BUGFIX: this was originally TMP
int check_finite_rotation_hamiltonian( const constants constants, state_p state_p, config_p config_p, int angle ){

  /* constants */
  int     N_coor;
  /* dummies */
  vector  F_partial_rot;
  matrix  K_partial_rot;
  int info=0;


  N_coor = constants.N_coor;


  /* allocations */
  if( VECTOR_ALLOCATE( N_coor, F_partial_rot ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, K_partial_rot ) ) info=1;


  /* F_partial */
  if( COMPUTE_F_PARTIAL( constants, *state_p, *config_p, F_partial_rot ) ) info=1;

  if( CHECK_FINITE_ROTATION_VECTOR( constants, angle, F_partial_rot ) ) info=1;

  /* K_partial */
  if( COMPUTE_K_PARTIAL( constants, *state_p, *config_p, K_partial_rot ) ) info=1;

  if( CHECK_FINITE_ROTATION_MATRIX( constants, angle, K_partial_rot ) ) info=1;


  /* dellocation */
  if( VECTOR_FREE( F_partial_rot ) ) info=1;
    
  if( MATRIX_FREE( K_partial_rot ) ) info=1;
  

  return info;

}

//------------------------------------------

/* check_finite_rotation_forces */
//BUGFIX: this was originally TMP
int check_finite_rotation_forces( const constants constants, state_p state_p, config_p config_p, int angle ){

  /* constants */
  int     N_coor;
  /* state */
  rvector_p  forces_p;
  /* dummies */
  vector  forces_rot;
  int i;
  int info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: check rotation forces\n" );

#endif /* __DEBUG_PLUS__ */


  N_coor    = constants.N_coor;

  forces_p  = &(config_p->atoms.forces);


  /* allocations */
  if( VECTOR_ALLOCATE( N_coor, forces_rot ) ) info=1;

  for( i=0; i<N_coor; i++ ){

    forces_rot.vector[ i ] = CMPLX( forces_p->rvector[ i ] );

  }

  if( CENTRE_OF_MASS_TRANSFORM_VECTOR( constants, forces_rot ) ) info=1;

  if( CHECK_FINITE_ROTATION_VECTOR( constants, angle, forces_rot ) ) info=1;


  /* dellocation */
  if( VECTOR_FREE( forces_rot ) ) info=1;
    

  return info;

}

//------------------------------------------

/* check_finite_rotation_positions */
//BUGFIX: this was originally TMP
int check_finite_rotation_positions( const constants constants, state_p state_p, config_p config_p, int angle ){

  /* constants */
  int     N_coor;
  /* state */
  rvector_p  positions_p;
  /* dummies */
  vector  positions_rot;
  int i;
  int info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: check rotation positions\n" );

#endif /* __DEBUG_PLUS__ */


  N_coor       = constants.N_coor;

  positions_p  = &(config_p->atoms.positions);


  /* allocations */
  if( VECTOR_ALLOCATE( N_coor, positions_rot ) ) info=1;

  for( i=0; i<N_coor; i++ ){

    positions_rot.vector[ i ] = CMPLX( positions_p->rvector[ i ] );

  }

  if( CENTRE_OF_MASS_TRANSFORM_VECTOR( constants, positions_rot ) ) info=1;

  if( CHECK_FINITE_ROTATION_VECTOR( constants, angle, positions_rot ) ) info=1;


  /* dellocation */
  if( VECTOR_FREE( positions_rot ) ) info=1;
    

  return info;

}

//------------------------------------------

/* check_finite_rotation_momenta */
//BUGFIX: this was originally TMP
int check_finite_rotation_momenta( const constants constants, state_p state_p, config_p config_p, int angle ){

  /* constants */
  int     N_coor;
  /* state */
  rvector_p  momenta_p;
  /* dummies */
  vector  momenta_rot;
  int i;
  int info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: check rotation momenta\n" );

#endif /* __DEBUG_PLUS__ */


  N_coor     = constants.N_coor;

  momenta_p  = &(config_p->atoms.momenta);


  /* allocations */
  if( VECTOR_ALLOCATE( N_coor, momenta_rot ) ) info=1;

  for( i=0; i<N_coor; i++ ){

    momenta_rot.vector[ i ] = CMPLX( momenta_p->rvector[ i ] );

  }

  if( CENTRE_OF_MASS_TRANSFORM_VECTOR( constants, momenta_rot ) ) info=1;

  if( CHECK_FINITE_ROTATION_VECTOR( constants, angle, momenta_rot ) ) info=1;


  /* dellocation */
  if( VECTOR_FREE( momenta_rot ) ) info=1;
    

  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

/* check_finite_rotation_vector */

int check_finite_rotation_vector( const constants constants, int angle, vector_p vector_p ){

  /* constants */
  int    N_coor;
  /* dummies */
  matrix rotation;
  vector dummy_vector1;
  int info=0;


  N_coor = constants.N_coor;

  /* allocation */
  if( MATRIX_ALLOCATE( N_coor, N_coor, rotation ) ) info=1;

  if( VECTOR_ALLOCATE( N_coor, dummy_vector1 ) ) info=1;


  /* compute_finite_rotation*/
  if( COMPUTE_FINITE_ROTATION( constants, angle, rotation) ) info=1;


  if( MATRIX_VECTOR_PRODUCT( rotation, *vector_p, dummy_vector1 ) ) info=1;


  if( VECTOR_COMPARE( *vector_p, dummy_vector1) ){

    fprintf( stderr, "ERROR: the vector is not invariant with respect to finite rotation\n");
    fflush( stderr );

    info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "finite rotation\n" );
  if( MATRIX_PRINT_PLUS( stdout, rotation ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "original vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, *vector_p ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "transformed vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, dummy_vector1 ) ) info=1;
  fprintf( stdout, "\n" );

  if( info ){

    fprintf( stdout, "FAILURE!\n" );
    fprintf( stdout, "\n" );

  }

#endif /* __DEBUG_PLUS__ */


  /* deallocation */
  if( VECTOR_FREE( dummy_vector1 ) ) info=1;

  if( MATRIX_FREE( rotation ) ) info=1;


  return info;

}

//------------------------------------------

/* check_finite_double_rotation_vector */

int check_finite_double_rotation_vector( const constants constants, int angle, vector_p vector_p ){

  /* constants */
  int    N_coor;
  /* dummies */
  matrix rotation;
  vector dummy_vector1;
  vector dummy_vector2;
  int info=0;


  N_coor = constants.N_coor;

  /* allocation */
  if( MATRIX_ALLOCATE( N_coor, N_coor, rotation ) ) info=1;

  if( VECTOR_ALLOCATE( N_coor, dummy_vector1 ) ) info=1;

  if( VECTOR_ALLOCATE( N_coor, dummy_vector2 ) ) info=1;


  /* compute_finite_rotation*/
  if( COMPUTE_FINITE_ROTATION( constants, angle, rotation) ) info=1;


  if( MATRIX_VECTOR_PRODUCT( rotation, *vector_p, dummy_vector1 ) ) info=1;

  if( MATRIX_VECTOR_PRODUCT( rotation, dummy_vector1, dummy_vector2 ) ) info=1;


  if( VECTOR_COMPARE( *vector_p, dummy_vector2 ) ){

    fprintf( stderr, "ERROR: the vector is not invariant with respect to finite rotation\n");
    fflush( stderr );

    info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "finite rotation\n" );
  if( MATRIX_PRINT_PLUS( stdout, rotation ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "original vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, *vector_p ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "transformed vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, dummy_vector2 ) ) info=1;
  fprintf( stdout, "\n" );

  if( info ){

    fprintf( stdout, "FAILURE!\n" );
    fprintf( stdout, "\n" );
    
  }

#endif /* __DEBUG_PLUS__ */


  /* deallocation */
  if( VECTOR_FREE( dummy_vector2 ) ) info=1;

  if( VECTOR_FREE( dummy_vector1 ) ) info=1;

  if( MATRIX_FREE( rotation ) ) info=1;


  return info;

}

//------------------------------------------

/* check_finite_rotation_matrix */

int check_finite_rotation_matrix( const constants constants, int angle, matrix_p matrix_p ){

  /* constants */
  int    N_coor;
  /* dummies */
  matrix rotation;
  matrix dummy_matrix1;
  matrix dummy_matrix2;
  matrix dummy_matrix3;
  int info=0;


  N_coor = constants.N_coor;

  /* allocation */
  if( MATRIX_ALLOCATE( N_coor, N_coor, rotation ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix1 ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix2 ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix3 ) ) info=1;


  /* compute_finite_rotation*/
  if( COMPUTE_FINITE_ROTATION( constants, angle, rotation) ) info=1;

  if( MATRIX_ADJOINT( rotation, dummy_matrix2 ) ) info=1;

  //  if( MATRIX_MATRIX_PRODUCT( rotation, dummy_matrix2, dummy_matrix1 ) ) info=1;

  /*
    fprintf( stdout, "unit matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, dummy_matrix1 ) ) info=1;
    fprintf( stdout, "\n" );
  */

  //
  if( MATRIX_MATRIX_PRODUCT( rotation, *matrix_p, dummy_matrix1 ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( dummy_matrix1, dummy_matrix2, dummy_matrix3  ) ) info=1;

  /*
    if( MATRIX_MATRIX_PRODUCT( *matrix_p, rotation, dummy_matrix1 ) ) info=1;
    
    if( MATRIX_MATRIX_PRODUCT( dummy_matrix2, dummy_matrix1, dummy_matrix3  ) ) info=1;
  */

  if( MATRIX_COMPARE( *matrix_p, dummy_matrix3 ) ){

    fprintf( stderr, "ERROR: the matrix is not invariant with respect to finite rotation\n");
    fflush( stderr );

    info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "finite rotation\n" );
  if( MATRIX_PRINT_PLUS( stdout, rotation ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "original matrix\n" );
  if( MATRIX_PRINT_PLUS( stdout, *matrix_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "transformed matrix\n" );
  if( MATRIX_PRINT_PLUS( stdout, dummy_matrix3 ) ) info=1;
  fprintf( stdout, "\n" );

  if( info ){

    fprintf( stdout, "FAILURE!\n" );
    fprintf( stdout, "\n" );

  }

#endif /* __DEBUG_PLUS__ */


  /* deallocation */
  if( MATRIX_FREE( dummy_matrix3 ) ) info=1;

  if( MATRIX_FREE( dummy_matrix2 ) ) info=1;

  if( MATRIX_FREE( dummy_matrix1 ) ) info=1;

  if( MATRIX_FREE( rotation ) ) info=1;


  return info;

}

//------------------------------------------
//------------------------------------------

/* check_finite_rotation_vector_comp */

int check_finite_rotation_vector_comp( const constants constants, int angle, vector_p vector_p ){

  /* constants */
  int    N_coor;
  /* dummies */
  matrix rotation;
  vector dummy_vector1;
  int info=0;


  N_coor = constants.N_coor;

  /* allocation */
  if( MATRIX_ALLOCATE( N_coor, N_coor, rotation ) ) info=1;

  if( VECTOR_ALLOCATE( N_coor, dummy_vector1 ) ) info=1;


  /* compute_finite_rotation_comp */
  if( COMPUTE_FINITE_ROTATION_COMP( constants, angle, rotation) ) info=1;


  if( MATRIX_VECTOR_PRODUCT( rotation, *vector_p, dummy_vector1 ) ) info=1;


  if( VECTOR_COMPARE( *vector_p, dummy_vector1) ){

    fprintf( stderr, "ERROR: the vector is not invariant with respect to finite rotation\n");
    fflush( stderr );

    info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "finite rotation_comp\n" );
  if( MATRIX_PRINT_PLUS( stdout, rotation ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "original vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, *vector_p ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "transformed vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, dummy_vector1 ) ) info=1;
  fprintf( stdout, "\n" );

  if( info ){

    fprintf( stdout, "FAILURE!\n" );
    fprintf( stdout, "\n" );

  }

#endif /* __DEBUG_PLUS__ */


  /* deallocation */
  if( VECTOR_FREE( dummy_vector1 ) ) info=1;

  if( MATRIX_FREE( rotation ) ) info=1;


  return info;

}

//------------------------------------------

/* check_finite_double_rotation_vector_comp */

int check_finite_double_rotation_vector_comp( const constants constants, int angle, vector_p vector_p ){

  /* constants */
  int    N_coor;
  /* dummies */
  matrix rotation;
  vector dummy_vector1;
  vector dummy_vector2;
  int info=0;


  N_coor = constants.N_coor;

  /* allocation */
  if( MATRIX_ALLOCATE( N_coor, N_coor, rotation ) ) info=1;

  if( VECTOR_ALLOCATE( N_coor, dummy_vector1 ) ) info=1;

  if( VECTOR_ALLOCATE( N_coor, dummy_vector2 ) ) info=1;


  /* compute_finite_rotation_comp */
  if( COMPUTE_FINITE_ROTATION_COMP( constants, angle, rotation) ) info=1;


  if( MATRIX_VECTOR_PRODUCT( rotation, *vector_p, dummy_vector1 ) ) info=1;

  if( MATRIX_VECTOR_PRODUCT( rotation, dummy_vector1, dummy_vector2 ) ) info=1;


  if( VECTOR_COMPARE( *vector_p, dummy_vector2 ) ){

    fprintf( stderr, "ERROR: the vector is not invariant with respect to finite rotation\n");
    fflush( stderr );

    info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "finite rotation_comp\n" );
  if( MATRIX_PRINT_PLUS( stdout, rotation ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "original vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, *vector_p ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "transformed vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, dummy_vector2 ) ) info=1;
  fprintf( stdout, "\n" );
  
  if( info ){

    fprintf( stdout, "FAILURE!\n" );
    fprintf( stdout, "\n" );

  }

#endif /* __DEBUG_PLUS__ */


  /* deallocation */
  if( VECTOR_FREE( dummy_vector2 ) ) info=1;

  if( VECTOR_FREE( dummy_vector1 ) ) info=1;

  if( MATRIX_FREE( rotation ) ) info=1;


  return info;

}

//------------------------------------------

/* check_finite_rotation_matrix_comp */

int check_finite_rotation_matrix_comp( const constants constants, int angle, matrix_p matrix_p ){

  /* constants */
  int    N_coor;
  /* dummies */
  matrix rotation;
  matrix dummy_matrix1;
  matrix dummy_matrix2;
  matrix dummy_matrix3;
  int info=0;


  N_coor = constants.N_coor;

  /* allocation */
  if( MATRIX_ALLOCATE( N_coor, N_coor, rotation ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix1 ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix2 ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix3 ) ) info=1;


  /* compute_finite_rotation_comp*/
  if( COMPUTE_FINITE_ROTATION_COMP( constants, angle, rotation) ) info=1;

  if( MATRIX_ADJOINT( rotation, dummy_matrix2 ) ) info=1;

  //  if( MATRIX_MATRIX_PRODUCT( rotation, dummy_matrix2, dummy_matrix1 ) ) info=1;

  /*
    fprintf( stdout, "unit matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, dummy_matrix1 ) ) info=1;
    fprintf( stdout, "\n" );
  */

  //
  if( MATRIX_MATRIX_PRODUCT( rotation, *matrix_p, dummy_matrix1 ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( dummy_matrix1, dummy_matrix2, dummy_matrix3  ) ) info=1;

  /*
    if( MATRIX_MATRIX_PRODUCT( *matrix_p, rotation, dummy_matrix1 ) ) info=1;
    
    if( MATRIX_MATRIX_PRODUCT( dummy_matrix2, dummy_matrix1, dummy_matrix3  ) ) info=1;
  */

  if( MATRIX_COMPARE( *matrix_p, dummy_matrix3 ) ){

    fprintf( stderr, "ERROR: the matrix is not invariant with respect to finite rotation\n");
    fflush( stderr );

    info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "finite rotation_comp\n" );
  if( MATRIX_PRINT_PLUS( stdout, rotation ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "original matrix\n" );
  if( MATRIX_PRINT_PLUS( stdout, *matrix_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "transformed matrix\n" );
  if( MATRIX_PRINT_PLUS( stdout, dummy_matrix3 ) ) info=1;
  fprintf( stdout, "\n" );

  if( info ){

    fprintf( stdout, "FAILURE!\n" );
    fprintf( stdout, "\n" );

  }

#endif /* __DEBUG_PLUS__ */
  

  /* deallocation */
  if( MATRIX_FREE( dummy_matrix3 ) ) info=1;

  if( MATRIX_FREE( dummy_matrix2 ) ) info=1;

  if( MATRIX_FREE( dummy_matrix1 ) ) info=1;

  if( MATRIX_FREE( rotation ) ) info=1;


  return info;

}

//------------------------------------------
//------------------------------------------

/* compute_finite_rotation */

int compute_finite_rotation( const constants constants , int angle, matrix_p matrix_p ){

  /* constants*/
  int    N_atoms;
  /* dummies */
  int sdim;
  int i;
  int info=0;


  N_atoms = constants.N_atoms;
  sdim    =  constants.spacial_dimension; 

  if( MATRIX_ZERO( *matrix_p ) ) info=1;


  if( sdim != 2 ){

    fprintf( stderr, "ERROR: It's implemented only for 2 dimension\n" );

    info=1;

  }

  if( matrix_p->matrix_dim_row != N_atoms *sdim || matrix_p->matrix_dim_column != N_atoms *sdim ){

    fprintf( stderr, "ERROR: the matrix has not the right dimesions\n" );
    fflush( stderr );

    info=1;

  }


  if( !info){

    for( i=0; i<N_atoms; i++ ){

      matrix_p->matrix[ COORDINATE_INDEX( 0 +i*sdim, 0 +PBC_INT( i-1, N_atoms ) *sdim) ] = CMPLX(  cos( 2 *PI/angle ) );
      matrix_p->matrix[ COORDINATE_INDEX( 0 +i*sdim, 1 +PBC_INT( i-1, N_atoms ) *sdim) ] = CMPLX( -sin( 2 *PI/angle ) );

      matrix_p->matrix[ COORDINATE_INDEX( 1 +i*sdim, 0 +PBC_INT( i-1, N_atoms ) *sdim) ] = CMPLX(  sin( 2 *PI/angle ) );
      matrix_p->matrix[ COORDINATE_INDEX( 1 +i*sdim, 1 +PBC_INT( i-1, N_atoms ) *sdim) ] = CMPLX(  cos( 2 *PI/angle ) );

    }

  }

  /*
    fprintf( stdout, "finite rotation\n" );
    if( MATRIX_PRINT_PLUS( stdout, *matrix_p ) ) info=1;
    fprintf( stdout, "\n" );
  */


  return info;

}

//------------------------------------------

/* compute_finite_rotation_comp */

int compute_finite_rotation_comp( const constants constants , int angle, matrix_p matrix_p ){

  /* constants*/
  int    N_atoms;
  /* dummies */
  int sdim;
  int i;
  int info=0;


  N_atoms = constants.N_atoms;
  sdim    =  constants.spacial_dimension; 


  if( MATRIX_ZERO( *matrix_p ) ) info=1;


  if( sdim != 2 ){

    fprintf( stderr, "ERROR: It's implemented only for 2 dimension\n" );

    info=1;

  }

  if( matrix_p->matrix_dim_row != N_atoms *sdim || matrix_p->matrix_dim_column != N_atoms *sdim ){

    fprintf( stderr, "ERROR: the matrix has not the right dimesions\n" );
    fflush( stderr );

    info=1;

  }


  if( !info){

    for( i=0; i<N_atoms; i++ ){

      matrix_p->matrix[ COORDINATE_INDEX( 0 +i*sdim, 0 +PBC_INT( i+1, N_atoms ) *sdim) ] = CMPLX(  cos( 2 *PI/angle ) );
      matrix_p->matrix[ COORDINATE_INDEX( 0 +i*sdim, 1 +PBC_INT( i+1, N_atoms ) *sdim) ] = CMPLX( -sin( 2 *PI/angle ) );

      matrix_p->matrix[ COORDINATE_INDEX( 1 +i*sdim, 0 +PBC_INT( i+1, N_atoms ) *sdim) ] = CMPLX(  sin( 2 *PI/angle ) );
      matrix_p->matrix[ COORDINATE_INDEX( 1 +i*sdim, 1 +PBC_INT( i+1, N_atoms ) *sdim) ] = CMPLX(  cos( 2 *PI/angle ) );

    }

  }

  /*
    fprintf( stdout, "finite rotation_comp\n" );
    if( MATRIX_PRINT_PLUS( stdout, *matrix_p ) ) info=1;
    fprintf( stdout, "\n" );
  */


  return info;

}


//------------------------------------------
//------------------------------------------

/* centre_of_mass_transform_vector */

int centre_of_mass_transform_vector( const constants constants, vector_p vector_p ){

  /* constants */
  int        N_atoms;
  /* dummies */
  int        sdim;
  int        i,j;
  vector     centroid;
  int        info=0;


  N_atoms = constants.N_atoms;
  sdim    =  constants.spacial_dimension; 

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: centre_of_mass_transform\n" );

#endif /* __DEBUG_PLUS__ */


  /* allocate */
  if( VECTOR_ALLOCATE( sdim, centroid ) ) info=1;

  // centroid
  for( i=0; i<N_atoms; i++ ){
    
    for( j=0; j<sdim; j++ ){
      
      centroid.vector[ j ] = CMPLX_SUM( centroid.vector[ j ], vector_p->vector[ j +i *sdim ] );
      
    }
    
  }
  
  for( j=0; j<sdim; j++ ){
    
    centroid.vector[ j ] = CMPLX_DIVISION( centroid.vector[ j ], CMPLX( (double)N_atoms ) );
    
  }
  

  for( i=0; i<N_atoms; i++ ){
    
    for( j=0; j<sdim; j++ ){
      
      vector_p->vector[ j +i *sdim ] = CMPLX_DIF( vector_p->vector[ j +i *sdim ], centroid.vector[ j ] );
      
    }
    
  }


  /* allocate */
  if( VECTOR_FREE( centroid ) ) info=1;


  return info;

}

//------------------------------------------
//------------------------------------------

/* symmetry_adapted_unitary_transform  */

int symmetry_adapted_unitary_transform( const constants constants, double angle, matrix_p matrix_p ){

  /* constants */
  int    N_coor;
  /* dummies */
  matrix rotation;
  matrix dummy_matrix1;
  int  info=0;


  N_coor = constants.N_coor;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: symmetry_adapted_unitary_transform\n" );

#endif /* __DEBUG_PLUS__ */

  /* allocation */
  if( MATRIX_ALLOCATE( N_coor, N_coor, rotation ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix1 ) ) info=1;


  /* create rotation */

  if( N_coor != 6 ){

    fprintf( stderr, "ERROR: this function works only for a trimer in dim 2\n" );
    fflush( stderr );

    info=1;

  }

  if( MATRIX_UNIT( rotation ) ) info=1;


  rotation.matrix[ COORDINATE_INDEX( 3, 3 ) ] = rotation.matrix[ COORDINATE_INDEX( 4, 4 ) ] = CMPLX( cos( angle ) ); 

  rotation.matrix[ COORDINATE_INDEX( 3, 4 ) ] = CMPLX( -sin( angle ) );
    
  rotation.matrix[ COORDINATE_INDEX( 4, 3 ) ] = CMPLX(  sin( angle ) ); 


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "symmetry_adapted_unitary_transform\n");
  if( MATRIX_PRINT_PLUS( stdout, rotation ) ) info=1;
  fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */

  /* unitary transform */

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "matrix [before]\n");
  if( MATRIX_PRINT_PLUS( stdout, *matrix_p ) ) info=1;
  fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */

  if( MATRIX_MATRIX_PRODUCT( *matrix_p, rotation, dummy_matrix1 ) ) info=1;

  if( MATRIX_COPY( *matrix_p, dummy_matrix1 ) ) info=1;

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "matrix [after]\n");
  if( MATRIX_PRINT_PLUS( stdout, *matrix_p ) ) info=1;
  fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */

  /* deallocation */
  if( MATRIX_FREE( dummy_matrix1 ) ) info=1;

  if( MATRIX_FREE( rotation ) ) info=1;


  return info;

}

//------------------------------------------
