
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
#include "PolyCEID_hamiltonian_many.h"


#define EPS0     1.0e-15

#define DELTA1   1.5e-6
#define EPS1     1.5e-7

#define DELTA2   1.5e-4
#define EPS2     1.5e-5


/*********************
  FUNCTIONS & MACROS
*********************/


/* hamiltonian_many update */

int hamiltonian_many_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  /* states */
  /* dummies */
  int                 comp;
  int                 info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: hamiltonian_many_update\n" );

#endif /* __DEBUG_PLUS__*/


#if defined( __DEBUG__ ) && defined( __HAMILTONIAN_CONSISTENCY__ )

  /* check consistency */
  if( HAMILTONIAN_SINGLE_CHECK( constants, *state_p, *config_p ) ) info=1;

#endif /* defined( __DEBUG__ ) && defined( __HAMILTONIAN_CONSISTENCY__ ) */


  /* distances update */
  // WARNING: and it is needed by the next
  if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;


  /* centre of mass update */
  // WARNING: and it is needed by the next
  if( COMPUTE_CENTRE_OF_MASS( constants, config_p->atoms ) ) info=1;


  /* H_matrix_many update */
  if( H_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;
      
    
  /* F_matrix_many update */
  if( F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p  ) ) info=1;
  

  /* K_matrix_many update */
  if( K_MATRIX_MANY_UPDATE( constants, *state_p, *config_p  ) ) info=1;
   

#if defined( __DEBUG__ ) && defined( __HAMILTONIAN_CONSISTENCY__ )
    
  /* check consistency */
  if( HAMILTONIAN_MANY_CHECK( constants, *state_p, *config_p ) ) info=1;

#endif /* defined( __DEBUG__ ) && defined( __HAMILTONIAN_CONSISTENCY__ ) */


  // WARNING: here there should be the loop over cartesian components...
  comp = 0;

  /* K_matrix_many update */
  if( CLASSICAL_DIPOLE_MANY_UPDATE( constants, *state_p, *config_p, comp ) ) info=1;


  /* forces update */
  // WARNING: it must be here because it depends on the previous matrices
  // and it is needed by the next
  if( FORCES_UPDATE( constants, *state_p, *config_p ) ) info=1;


  /* delta_F_matrix_many update */
  if( DELTA_F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;
  

  return info;

}

//------------------------------------------
//------------------------------------------

/* H_matrix_many update */

int H_matrix_many_update( const constants constants, state_p state_p, config_p config_p ){

  /* constanst */
  int        N_levels_single;
  int        N_levels_many;
  imatrix    single_level_occupation;
  imatrix    many_body_hybridisation_table;
  vector     symmetry_coefficient;
  /* state */
  matrix_p   Ehrenfest_frame_p;
  matrix_p   H_matrix_p;
  vector_p   dummy_vector_single1_p;
  vector_p   dummy_vector_single2_p;
  vector_p   dummy_vector_single3_p;
  matrix_p   dummy_matrix_single1_p;

  /* dummies */
  complex    dummy;
  double     sum;
  complex    symmetry_coefficient_loc;
  int        level1, level2;
  int        counter;
  int        i, j, k;
  int        i_single;
  int        info=0;


  N_levels_single               =  constants.N_levels_single;
  N_levels_many                 =  constants.N_levels_many;
  single_level_occupation       =  constants.single_level_occupation;
  many_body_hybridisation_table =  constants.many_body_hybridisation_table;
  symmetry_coefficient          =  constants.symmetry_coefficient;

  H_matrix_p                    = &(config_p->electrons.H_matrix);
  Ehrenfest_frame_p             = &(config_p->electrons.Ehrenfest_frame);
  dummy_vector_single1_p        = &(state_p->dummy_vector_single1);
  dummy_vector_single2_p        = &(state_p->dummy_vector_single2);
  dummy_vector_single3_p        = &(state_p->dummy_vector_single3);
  dummy_matrix_single1_p        = &(state_p->dummy_matrix_single1);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: H_matrix_many_update\n" );

#endif /* __DEBUG_PLUS__ */


  /* H_matrix_single*/
  if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, *dummy_matrix_single1_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "H_matrix_single\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* set matrix to zero */
  if( MATRIX_ZERO( *H_matrix_p ) ) info=1; // Important


  /* set counter to zero */
  counter=0;

  for( i=0; i<N_levels_many; i++ ){

    // diagonal terms

    /* set dummy to zero*/
    sum = 0.0e0;

    for( k=0; k<N_levels_single*SPIN_DEG; k++ ){

      //      fprintf( stdout, "single_level %d, single_level_occupation %d\n", k, single_level_occupation.imatrix[i][ k ] );

      if( single_level_occupation.imatrix[i][ k ] ){

	level1 = k /SPIN_DEG;

	//	fprintf( stdout, "OK, level1 = %d\n", level1 );

	/* vector copy */
	for( i_single=0; i_single<N_levels_single; i_single++ ){
	  
	  dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];

	} /* end i_single loop */


	/* compute matrix element */
	if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single2_p ) ) info=1;

	dummy = SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single2_p );

	/* dummy update*/
	sum += REAL( dummy );
	  
      } /* end conditional */

    } /* end k loop */

    /* writing */
    H_matrix_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ] = CMPLX( sum );


    for( j=(i+1); j<N_levels_many; j++ ){


      if( many_body_hybridisation_table.imatrix[counter][0] != i ||
	  many_body_hybridisation_table.imatrix[counter][1] != j ){

	fprintf( stderr, "ERROR: inconsistency found in many_body_hybridisation_table ordering [%d,%d]\n", i, j );
	fflush( stderr );

	info=1;

      }
      else{

	level1 = many_body_hybridisation_table.imatrix[counter][2];
	level2 = many_body_hybridisation_table.imatrix[counter][4];
	symmetry_coefficient_loc = symmetry_coefficient.vector[counter];

	counter++;


	if( level1 != N_levels_single && level2 != N_levels_single ){


	  /* vector copy */
	  for( i_single=0; i_single<N_levels_single; i_single++ ){

	    dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];
	    
	    dummy_vector_single2_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level2 ) ];
	    
	  } /* end i_single loop */
	  
	  
	  /* compute matrix element */
	  if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single2_p, *dummy_vector_single3_p ) ) info=1;
	  
	  dummy = CMPLX_PRODUCT( symmetry_coefficient_loc, SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single3_p ) );
	  
	  /* writing */
	  H_matrix_p->matrix[ ELECTRON_MANY_INDEX( i, j ) ] = dummy;
	  
	  H_matrix_p->matrix[ ELECTRON_MANY_INDEX( j, i ) ] = CONJ( dummy );
	  
	}

      } /* end conditional */

    } /* end j loop */

  } /* end i loop */


#ifdef __DEBUG__

  // if( CHECK_PARTICLE_HOLE_MATRIX( config_p->electrons.H_matrix ) ) info=1;

#endif /* __DEBUG__ */
  

  return info;

}

//------------------------------------------

/* F_matrix_many update */

int F_matrix_many_update( const constants constants, state_p state_p, config_p config_p ){

  /* constanst */
  int             N_coor;
  int             N_levels_single;
  int             N_levels_many;
  imatrix         single_level_occupation;
  imatrix         many_body_hybridisation_table;
  vector          symmetry_coefficient;
  /* state */
  matrix_p        Ehrenfest_frame_p;
  matrix_array_p  F_matrix_p;
  vector_p        dummy_vector_single1_p;
  vector_p        dummy_vector_single2_p;
  vector_p        dummy_vector_single3_p;
  matrix_p        dummy_matrix_single1_p;
  /* dummies */
  complex         dummy;
  double          sum;
  complex         symmetry_coefficient_loc;
  int             level1, level2;
  int             counter;
  int             i, j, k;
  int             i_single;
  int             i_coor;
  int             info=0;


  N_coor                        =  constants.N_coor;
  N_levels_single               =  constants.N_levels_single;
  N_levels_many                 =  constants.N_levels_many;
  single_level_occupation       =  constants.single_level_occupation;
  many_body_hybridisation_table =  constants.many_body_hybridisation_table;
  symmetry_coefficient          =  constants.symmetry_coefficient;

  F_matrix_p                    = &(config_p->electrons.F_matrix);
  Ehrenfest_frame_p             = &(config_p->electrons.Ehrenfest_frame);
  dummy_vector_single1_p        = &(state_p->dummy_vector_single1);
  dummy_vector_single2_p        = &(state_p->dummy_vector_single2);
  dummy_vector_single3_p        = &(state_p->dummy_vector_single3);
  dummy_matrix_single1_p        = &(state_p->dummy_matrix_single1);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: F_matrix_many_update\n" );

#endif /* __DEBUG_PLUS__ */


  /* set matrix array to zero */
  if( MATRIX_ARRAY_ZERO( *F_matrix_p ) ) info=1; // Important


  for( i_coor=0; i_coor<N_coor; i_coor++ ){

    /* F_matrix_single*/
    if( F_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, i_coor, *dummy_matrix_single1_p ) ) info=1;


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "F_matrix_single[%d]\n", i_coor );
    if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
    fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */


    /* set counter to zero */
    counter=0;

    for( i=0; i<N_levels_many; i++ ){
      
      // diagonal terms

      /* set dummy to zero*/
      sum = 0.0e0;

      for( k=0; k<N_levels_single*SPIN_DEG; k++ ){

	if( single_level_occupation.imatrix[i][ k ] ){
	  
	  level1 = k /SPIN_DEG;

	  /* vector copy */
	  for( i_single=0; i_single<N_levels_single; i_single++ ){
	  
	    dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];

	  } /* end i_single loop */


	  /* compute matrix element */
	  if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single2_p ) ) info=1;

	  dummy = SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single2_p );

	  /* dummy update*/
	  sum += REAL( dummy );
	  
	} /* end conditional */

      } /* end k loop */

      /* writing */
      F_matrix_p->array[ i_coor ].matrix[ ELECTRON_MANY_INDEX( i, i ) ] = CMPLX( sum );

      
      // off diagonal terms

      for( j=(i+1); j<N_levels_many; j++ ){


	if( many_body_hybridisation_table.imatrix[counter][0] != i ||
	    many_body_hybridisation_table.imatrix[counter][1] != j ){

	  fprintf( stderr, "ERROR: inconsistency found in many_body_hybridisation_table ordering [%d,%d]\n", i, j );
	  fflush( stderr );

	  info=1;

	}
	else{

	  level1 = many_body_hybridisation_table.imatrix[counter][2];
	  level2 = many_body_hybridisation_table.imatrix[counter][4];
	  symmetry_coefficient_loc = symmetry_coefficient.vector[counter];

	  counter++;


	  if( level1 != N_levels_single && level2 != N_levels_single ){


	    /* vector copy */
	    for( i_single=0; i_single<N_levels_single; i_single++ ){

	      dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];
	      
	      dummy_vector_single2_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level2 ) ];

	    } /* end i_single loop */


	    /* compute matrix element */
	    if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single2_p, *dummy_vector_single3_p ) ) info=1;
	  
	    dummy = CMPLX_PRODUCT( symmetry_coefficient_loc, SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single3_p ) );

	    /* writing */
	    F_matrix_p->array[ i_coor ].matrix[ ELECTRON_MANY_INDEX( i, j ) ] = dummy;
	    
	    F_matrix_p->array[ i_coor ].matrix[ ELECTRON_MANY_INDEX( j, i ) ] = CONJ( dummy );

	  }

	} /* end conditional */

      } /* end j loop */

    } /* end i loop */

#ifdef __DEBUG__

  // if( CHECK_PARTICLE_HOLE_MATRIX( state_p->electrons.F_matrix.array[ i_coor ] ) ) info=1;

#endif /* __DEBUG__ */


  } /* end i_coor loop */


  return info;

}

//------------------------------------------

/* delta_F_matrix_many update */
/* WARNING: delta_F_matrix_many update assumes that forces have been updated */

int delta_F_matrix_many_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        N_coor;
  int        N_levels_many;
  /* state */
  double*    force_p;
  matrix_p   F_matrix_many_p;
  matrix_p   delta_F_matrix_many_p;
  /* dummies */
  int        i;
  int        r;
  int        info=0;

  N_coor               =  constants.N_coor;
  N_levels_many        =  constants.N_levels_many;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: delta_F_matrix_many_update\n" );

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    for( r=0; r<N_coor; r++ ){ /* coordinate loop */

      force_p               = &config_p->atoms.forces.rvector[ r ];
      F_matrix_many_p       = &config_p->electrons.F_matrix.array[ r ];       // Warning: only one matrix updated!
      delta_F_matrix_many_p = &config_p->electrons.delta_F_matrix.array[ r ]; // Warning: only one matrix updated!

      /* delta_F_matrix_many, 1st step */
      if( MATRIX_COPY( *delta_F_matrix_many_p, *F_matrix_many_p ) ) info=1;
      
      /* delta_F_matrix_many, 2nd step */
      for( i=0;i<N_levels_many;i++){

	// real part only
	delta_F_matrix_many_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ].z[0] -= *force_p;
	
      }


#ifdef __DEBUG__

      // if( CHECK_PARTICLE_HOLE_MATRIX( state_p->electrons.delta_F_matrix.array[ r ] ) ) info=1;

#endif /* __DEBUG__ */


    }  /* coordinate loop */

  }


  return info;

}

//------------------------------------------

/* K_matrix_many update */

int K_matrix_many_update( const constants constants, state_p state_p, config_p config_p ){

  /* constanst */
  int             N_coor;
  int             N_levels_single;
  int             N_levels_many;
  imatrix         single_level_occupation;
  imatrix         many_body_hybridisation_table;
  vector          symmetry_coefficient;
  /* state */
  matrix_p        Ehrenfest_frame_p;
  matrix_array_p  K_matrix_p;
  vector_p        dummy_vector_single1_p;
  vector_p        dummy_vector_single2_p;
  vector_p        dummy_vector_single3_p;
  matrix_p        dummy_matrix_single1_p;
  /* dummies */
  complex         dummy;
  double          sum;
  complex         symmetry_coefficient_loc;
  int             level1, level2;
  int             counter;
  int             i, j, k;
  int             i_single;
  int             i_coor1;
  int             i_coor2;
  int             K_matrix_index;
  int             info=0;


  N_coor                        =  constants.N_coor;
  N_levels_single               =  constants.N_levels_single;
  N_levels_many                 =  constants.N_levels_many;
  single_level_occupation       =  constants.single_level_occupation;
  many_body_hybridisation_table =  constants.many_body_hybridisation_table;
  symmetry_coefficient          =  constants.symmetry_coefficient;

  K_matrix_p                    = &(config_p->electrons.K_matrix);
  Ehrenfest_frame_p             = &(config_p->electrons.Ehrenfest_frame);
  dummy_vector_single1_p        = &(state_p->dummy_vector_single1);
  dummy_vector_single2_p        = &(state_p->dummy_vector_single2);
  dummy_vector_single3_p        = &(state_p->dummy_vector_single3);
  dummy_matrix_single1_p        = &(state_p->dummy_matrix_single1);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: K_matrix_many_update\n" );

#endif /* __DEBUG_PLUS__ */


  /* set matrix array to zero */
  if( MATRIX_ARRAY_ZERO( *K_matrix_p ) ) info=1; // Important


  for( i_coor1=0; i_coor1<N_coor; i_coor1++ ){

    for( i_coor2=0; i_coor2<N_coor; i_coor2++ ){

      K_matrix_index = COORDINATE_INDEX( i_coor1, i_coor2 );

      /* K_matrix_single*/
      if( K_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, i_coor1, i_coor2, *dummy_matrix_single1_p ) ) info=1;


#ifdef __DEBUG_PLUS__

      fprintf( stdout, "K_matrix_single[%d,%d]\n", i_coor1, i_coor2 );
      if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
      fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */

    
      /* set counter to zero */
      counter=0;

      for( i=0; i<N_levels_many; i++ ){
      
	// diagonal terms

	/* set dummy to zero*/
	sum = 0.0e0;

	for( k=0; k<N_levels_single*SPIN_DEG; k++ ){

	  if( single_level_occupation.imatrix[ i ][ k ] ){
	  
	    level1 = k /SPIN_DEG;

	    /* vector copy */
	    for( i_single=0; i_single<N_levels_single; i_single++ ){
	  
	      dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];

	    } /* end i_single loop */


	    /* compute matrix element */
	    if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single2_p ) ) info=1;

	    dummy = SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single2_p );

	    /* dummy update*/
	    sum += REAL( dummy );
	  
	  } /* end conditional */

	} /* end k loop */

	/* writing */
	K_matrix_p->array[ K_matrix_index ].matrix[ ELECTRON_MANY_INDEX( i, i ) ] = CMPLX( sum );

      
	// off diagonal terms

	for( j=(i+1); j<N_levels_many; j++ ){

	  
	  if( many_body_hybridisation_table.imatrix[counter][0] != i ||
	      many_body_hybridisation_table.imatrix[counter][1] != j ){

	    fprintf( stderr, "ERROR: inconsistency found in many_body_hybridisation_table ordering [%d,%d]\n", i, j );
	    fflush( stderr );

	    info=1;

	  }
	  else{

	    level1 = many_body_hybridisation_table.imatrix[counter][2];
	    level2 = many_body_hybridisation_table.imatrix[counter][4];
	    symmetry_coefficient_loc = symmetry_coefficient.vector[counter];

	    counter++;


	    if( level1 != N_levels_single && level2 != N_levels_single ){


	      /* vector copy */
	      for( i_single=0; i_single<N_levels_single; i_single++ ){
	    
		dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];
		
		dummy_vector_single2_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level2 ) ];
		
	      } /* end i_single loop */

	      
	      /* compute matrix element */
	      if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single2_p, *dummy_vector_single3_p ) ) info=1;
	    
	      dummy = CMPLX_PRODUCT( symmetry_coefficient_loc, SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single3_p ) );

	      /* writing */
	      K_matrix_p->array[ K_matrix_index ].matrix[ ELECTRON_MANY_INDEX( i, j ) ] = dummy;

	      K_matrix_p->array[ K_matrix_index ].matrix[ ELECTRON_MANY_INDEX( j, i ) ] = CONJ( dummy );

	    }

	  } /* end conditional */


	} /* end j loop */

      } /* end i loop */

#ifdef __DEBUG__

      // if( CHECK_PARTICLE_HOLE_MATRIX( state_p->electrons.K_matrix.array[ K_matrix_index ] ) ) info=1;

#endif /* __DEBUG__ */


    } /* end i_coor2 loop */

  } /* end i_coor1 loop */


  return info;

}

//------------------------------------------
//
/* classical_dipole_many update */
//BUGFIX: this was originally DEF
int classical_dipole_many_update( const constants constants, state_p state_p, config_p config_p, const int comp ){

  /* constanst */
  int        N_levels_single;
  int        N_levels_many;
  imatrix    single_level_occupation;
  imatrix    many_body_hybridisation_table;
  vector     symmetry_coefficient;
  /* state */
  matrix_p   dipole_many_p;
  matrix_p   Ehrenfest_frame_p;
  vector_p   dummy_vector_single1_p;
  vector_p   dummy_vector_single2_p;
  vector_p   dummy_vector_single3_p;
  matrix_p   dummy_matrix_single1_p;
  /* dummies */
  complex    dummy;
  double     sum;
  complex    symmetry_coefficient_loc;
  int        level1, level2;
  int        counter;
  int        i, j, k;
  int        i_single;
  int        info=0;


  N_levels_single               =  constants.N_levels_single;
  N_levels_many                 =  constants.N_levels_many;
  single_level_occupation       =  constants.single_level_occupation;
  many_body_hybridisation_table =  constants.many_body_hybridisation_table;
  symmetry_coefficient          =  constants.symmetry_coefficient;

  dipole_many_p                 = &(config_p->electrons.dipole_many);
  Ehrenfest_frame_p             = &(config_p->electrons.Ehrenfest_frame);
  dummy_vector_single1_p        = &(state_p->dummy_vector_single1);
  dummy_vector_single2_p        = &(state_p->dummy_vector_single2);
  dummy_vector_single3_p        = &(state_p->dummy_vector_single3);
  dummy_matrix_single1_p        = &(state_p->dummy_matrix_single1);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: classical_dipole_matrix_update\n" );

#endif /* __DEBUG_PLUS__ */

  if( strcmp( constants.hamiltonian.class, "CHAIN_MOD2" ) && 
      strcmp( constants.hamiltonian.class, "SSH" ) ){ 

    fprintf( stderr, "# WARNING: dipole perturbation for this model_hamiltonian is not supported yet [%s]\n", constants.hamiltonian.class );
    fflush( stderr );

  }
  else{
	  
    /* classical_dipole_updating */
    if( CLASSICAL_DIPOLE_SINGLE_UPDATE( constants, *state_p, *config_p, comp, *dummy_matrix_single1_p ) ) info=1;


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "dipole_matrix_single\n" );
    if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
    fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


    /* set matrix to zero */
    if( MATRIX_ZERO( *dipole_many_p ) ) info=1; // Important


    /* set counter to zero */
    counter=0;

    for( i=0; i<N_levels_many; i++ ){

      // diagonal terms

      /* set dummy to zero*/
      sum = 0.0e0;

      for( k=0; k<N_levels_single*SPIN_DEG; k++ ){

        //      fprintf( stdout, "single_level %d, single_level_occupation %d\n", k, single_level_occupation.imatrix[i][ k ] );

        if( single_level_occupation.imatrix[i][ k ] ){

	  level1 = k /SPIN_DEG;

	  //	fprintf( stdout, "OK, level1 = %d\n", level1 );

	  /* vector copy */
	  for( i_single=0; i_single<N_levels_single; i_single++ ){
	  
	    dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];

	  } /* end i_single loop */


	  /* compute matrix element */
	  if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single2_p ) ) info=1;

	  dummy = SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single2_p );

	  /* dummy update*/
	  sum += REAL( dummy );
	  
        } /* end conditional */

      } /* end k loop */

      /* writing */
      dipole_many_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ] = CMPLX( sum );


      for( j=(i+1); j<N_levels_many; j++ ){


        if( many_body_hybridisation_table.imatrix[counter][0] != i ||
	    many_body_hybridisation_table.imatrix[counter][1] != j ){

	  fprintf( stderr, "ERROR: inconsistency found in many_body_hybridisation_table ordering [%d,%d]\n", i, j );
	  fflush( stderr );

	  info=1;

        }
        else{

	  level1 = many_body_hybridisation_table.imatrix[counter][2];
	  level2 = many_body_hybridisation_table.imatrix[counter][4];
	  symmetry_coefficient_loc = symmetry_coefficient.vector[counter];

	  counter++;


	  if( level1 != N_levels_single && level2 != N_levels_single ){


	    /* vector copy */
	    for( i_single=0; i_single<N_levels_single; i_single++ ){

	      dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];
	    
	      dummy_vector_single2_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level2 ) ];
	    
	    } /* end i_single loop */
	  
	  
	    /* compute matrix element */
	    if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single2_p, *dummy_vector_single3_p ) ) info=1;
	  
	    dummy = CMPLX_PRODUCT( symmetry_coefficient_loc, SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single3_p ) );
	  
	    /* writing */
	    dipole_many_p->matrix[ ELECTRON_MANY_INDEX( i, j ) ] = dummy;
	  
	    dipole_many_p->matrix[ ELECTRON_MANY_INDEX( j, i ) ] = CONJ( dummy );
	  
	  }

        } /* end conditional */

      } /* end j loop */

    } /* end i loop */

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "dipole_many\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dipole_many_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

#ifdef __DEBUG__

//------------------------------------------


/* H_matrix_many_check */

int H_matrix_many_check( const constants constants, state_p state_p, config_p config_p, rvector_p positions_p, matrix_p matrix_check_p ){

  /* constanst */
  int        N_levels_single;
  int        N_levels_many;
  imatrix    single_level_occupation;
  imatrix    many_body_hybridisation_table;
  vector     symmetry_coefficient;
  /* state */
  matrix_p   Ehrenfest_frame_p;
  vector_p   dummy_vector_single1_p;
  vector_p   dummy_vector_single2_p;
  vector_p   dummy_vector_single3_p;
  matrix_p   dummy_matrix_single1_p;

  /* dummies */
  complex    dummy;
  double     sum;
  complex    symmetry_coefficient_loc;
  int        level1, level2;
  int        counter;
  int        i, j, k;
  int        i_single;
  int        info=0;


  N_levels_single               =  constants.N_levels_single;
  N_levels_many                 =  constants.N_levels_many;
  single_level_occupation       =  constants.single_level_occupation;
  many_body_hybridisation_table =  constants.many_body_hybridisation_table;
  symmetry_coefficient          =  constants.symmetry_coefficient;

  Ehrenfest_frame_p             = &(config_p->electrons.Ehrenfest_frame);
  dummy_vector_single1_p        = &(state_p->dummy_vector_single1);
  dummy_vector_single2_p        = &(state_p->dummy_vector_single2);
  dummy_vector_single3_p        = &(state_p->dummy_vector_single3);
  dummy_matrix_single1_p        = &(state_p->dummy_matrix_single1);



  //  if( H_MATRIX_SINGLE_CHECK( constants, *positions_p, *matrix_check_p ) ) info=1;

  /* H_matrix_single_check */
  if( H_MATRIX_SINGLE_CHECK( constants, *positions_p, *dummy_matrix_single1_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "H_matrix_single\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
  fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */


  /* set matrix to zero */
  if( MATRIX_ZERO( *matrix_check_p ) ) info=1; // Important


  /* set counter to zero */
  counter=0;

  for( i=0; i<N_levels_many; i++ ){

    // diagonal terms

    /* set dummy to zero*/
    sum = 0.0e0;

    for( k=0; k<N_levels_single*SPIN_DEG; k++ ){

      if( single_level_occupation.imatrix[i][ k ] ){

	level1 = k /SPIN_DEG;

	/* vector copy */
	for( i_single=0; i_single<N_levels_single; i_single++ ){
	  
	  dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];

	} /* end i_single loop */


	/* compute matrix element */
	if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single2_p ) ) info=1;

	dummy = SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single2_p );

	/* dummy update*/
	sum += REAL( dummy );
	  
      } /* end conditional */

    } /* end k loop */

    /* writing */
    matrix_check_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ] = CMPLX( sum );

	
    // off diagonal terms

    for( j=(i+1); j<N_levels_many; j++ ){


      if( many_body_hybridisation_table.imatrix[counter][0] != i ||
	  many_body_hybridisation_table.imatrix[counter][1] != j ){

	fprintf( stderr, "ERROR: inconsistency found in many_body_hybridisation_table ordering [%d,%d]\n", i, j );
	fflush( stderr );

	info=1;

      }
      else{

	level1 = many_body_hybridisation_table.imatrix[counter][2];
	level2 = many_body_hybridisation_table.imatrix[counter][4];
        symmetry_coefficient_loc = symmetry_coefficient.vector[counter];

	counter++;


	if( level1 != N_levels_single && level2 != N_levels_single ){

	  /* vector copy */
	  for( i_single=0; i_single<N_levels_single; i_single++ ){

	    dummy_vector_single1_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level1 ) ];

	    dummy_vector_single2_p->vector[ i_single ] = Ehrenfest_frame_p->matrix[ ELECTRON_SINGLE_INDEX( i_single, level2 ) ];
	    
	  } /* end i_single loop */


	  /* compute matrix element */
	  if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single2_p, *dummy_vector_single3_p ) ) info=1;
	  
          dummy = CMPLX_PRODUCT( symmetry_coefficient_loc, SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single3_p ) );

	  /* writing */
	  matrix_check_p->matrix[ ELECTRON_MANY_INDEX( i, j ) ] = dummy;
	  
	  matrix_check_p->matrix[ ELECTRON_MANY_INDEX( j, i ) ] = CONJ( dummy );

	}

      } /* end conditional */

    } /* end j loop */

  } /* end i loop */


  return info;

}

//------------------------------------------

/* Hamiltonian_many_check */

int Hamiltonian_many_check( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  /* state */
  rvector_p       positions_p;
  ivector_p       mask_p;
  matrix_p        H_matrix_many_p;
  matrix_array_p  F_matrix_many_p;
  matrix_array_p  K_matrix_many_p;
  rvector_p       dummy_positions_p;
  matrix_p        dummy_matrix1_p;
  matrix_p        dummy_matrix2_p;
  matrix_p        dummy_matrix3_p;
  matrix_p        dummy_matrix4_p;
  /*dummies*/
  int i_coor1, i_coor2;
  int index;
  int info=0;


  N_coor                  =  constants.N_coor;

  positions_p             = &(config_p->atoms.positions);
  mask_p                  = &(config_p->atoms.mask);
  H_matrix_many_p         = &(config_p->electrons.H_matrix);
  F_matrix_many_p         = &(config_p->electrons.F_matrix);
  K_matrix_many_p         = &(config_p->electrons.K_matrix);
  dummy_positions_p       = &(state_p->dummy_positions);
  dummy_matrix1_p         = &(state_p->dummy_matrix1);
  dummy_matrix2_p         = &(state_p->dummy_matrix2);
  dummy_matrix3_p         = &(state_p->dummy_matrix3);
  dummy_matrix4_p         = &(state_p->dummy_matrix4);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: Hamiltonian_many_check\n" );

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    /* H_matrix_many check */
#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: H_matrix_many_check\n" );

#endif /* __DEBUG_PLUS__ */


    if( RVECTOR_COPY( *dummy_positions_p, *positions_p ) ) info=1;

    if( H_MATRIX_MANY_CHECK( constants, *state_p, *config_p, *dummy_positions_p, *dummy_matrix1_p ) ) info=1;

    if( MATRIX_COMPARE_ENTRIES( *dummy_matrix1_p, *H_matrix_many_p, EPS0  ) ) info=1;

    if( info ){

      fprintf( stderr, "----------\n" );
      fprintf( stderr, "ERROR in H_matrix_many check [The first one is the numerical estimate]\n" );
      fprintf( stderr, "----------\n" );

    }

  }

  if( !info ){

    /* F_matrix_many check */
#ifdef __DEBUG_PLUS__

    fprintf( stdout, "DOING: F_matrix_many_check\n" );

#endif /* __DEBUG_PLUS__ */


    for( i_coor1=0; i_coor1<N_coor; i_coor1++ ){

      if( info ) break;
      
      if( mask_p->ivector[ i_coor1 ] ){ 
      
        index = i_coor1;
      
        //    if( RVECTOR_COPY( *dummy_positions_p, *positions_p ) ) info=1;
      
        // first term;
      
        dummy_positions_p->rvector[ i_coor1 ] += DELTA1;
      
        if( H_MATRIX_MANY_CHECK( constants, *state_p, *config_p, *dummy_positions_p, *dummy_matrix1_p ) ) info=1;
      
        dummy_positions_p->rvector[ i_coor1 ] -= DELTA1;

        // second term;

        dummy_positions_p->rvector[ i_coor1 ] -= DELTA1;

        if( H_MATRIX_MANY_CHECK( constants, *state_p, *config_p, *dummy_positions_p, *dummy_matrix2_p ) ) info=1;
      
        dummy_positions_p->rvector[ i_coor1 ] += DELTA1;

        // final

        if( MATRIX_DIF( *dummy_matrix1_p, *dummy_matrix2_p, *dummy_matrix1_p ) ) info=1; // WARNING: overwriting
      
        /* WARNING: the F_matrix_many is MINUS the derivative of H_matrix_many */
        if( MATRIX_SCALAR_DIVISION( *dummy_matrix1_p, CMPLX( -2.0e0 *DELTA1 ), *dummy_matrix1_p ) ) info=1; // WARNING: overwriting

        if( MATRIX_COMPARE_ENTRIES( *dummy_matrix1_p, F_matrix_many_p->array[ index ], EPS1  ) ) info=1;

        if( info ){

	  fprintf( stderr, "----------\n" );
	  fprintf( stderr, "ERROR in F_matrix_many check [The first one is the numerical estimate]\n" );
	  fprintf( stderr, "---> coor: (%d)\n", i_coor1 );
	  fprintf( stderr, "----------\n" );

        }

      }  

    } /* i_coor1 loop */
    
  }

  if( !info ){

    /* K_matrix_many check */
#ifdef __DEBUG_PLUS__

    fprintf( stdout, "DOING: K_matrix_many_check\n" );

#endif /* __DEBUG_PLUS__ */


    for( i_coor1=0; i_coor1<N_coor; i_coor1++ ){
      
      for( i_coor2=0; i_coor2<N_coor; i_coor2++ ){

	if( info ) break;
	
        if( mask_p->ivector[ i_coor1 ] && mask_p->ivector[ i_coor2 ] ){ 
      
	  index = COORDINATE_INDEX( i_coor1, i_coor2 );

	  //      if( RVECTOR_COPY( *dummy_positions_p, *positions_p ) ) info=1;

	  // first term

	  dummy_positions_p->rvector[ i_coor1 ] += DELTA2;

	  dummy_positions_p->rvector[ i_coor2 ] += DELTA2;

	  if( H_MATRIX_MANY_CHECK( constants, *state_p, *config_p, *dummy_positions_p, *dummy_matrix1_p ) ) info=1;

	  dummy_positions_p->rvector[ i_coor2 ] -= DELTA2;
	
	  dummy_positions_p->rvector[ i_coor1 ] -= DELTA2;

	  // second term

	  dummy_positions_p->rvector[ i_coor1 ] += DELTA2;

	  dummy_positions_p->rvector[ i_coor2 ] -= DELTA2;

	  if( H_MATRIX_MANY_CHECK( constants, *state_p, *config_p, *dummy_positions_p, *dummy_matrix2_p ) ) info=1;

	  dummy_positions_p->rvector[ i_coor2 ] += DELTA2;

	  dummy_positions_p->rvector[ i_coor1 ] -= DELTA2;

	  // third term

	  dummy_positions_p->rvector[ i_coor1 ] -= DELTA2;

	  dummy_positions_p->rvector[ i_coor2 ] += DELTA2;

	  if( H_MATRIX_MANY_CHECK( constants, *state_p, *config_p, *dummy_positions_p, *dummy_matrix3_p ) ) info=1;

	  dummy_positions_p->rvector[ i_coor2 ] -= DELTA2;

	  dummy_positions_p->rvector[ i_coor1 ] += DELTA2;

	  // fourth term

	  dummy_positions_p->rvector[ i_coor1 ] -= DELTA2;

	  dummy_positions_p->rvector[ i_coor2 ] -= DELTA2;

	  if( H_MATRIX_MANY_CHECK( constants, *state_p, *config_p, *dummy_positions_p, *dummy_matrix4_p ) ) info=1;

	  dummy_positions_p->rvector[ i_coor2 ] += DELTA2;

	  dummy_positions_p->rvector[ i_coor1 ] += DELTA2;

	  // final

	  if( MATRIX_DIF( *dummy_matrix1_p, *dummy_matrix2_p, *dummy_matrix1_p ) ) info=1; // WARNING: overwriting

	  if( MATRIX_DIF( *dummy_matrix1_p, *dummy_matrix3_p, *dummy_matrix1_p ) ) info=1; // WARNING: overwriting

	  if( MATRIX_SUM( *dummy_matrix1_p, *dummy_matrix4_p, *dummy_matrix1_p ) ) info=1; // WARNING: overwriting

	  if( MATRIX_SCALAR_DIVISION( *dummy_matrix1_p, CMPLX( 4.0e0 *DELTA2 *DELTA2 ), *dummy_matrix1_p ) ) info=1; // WARNING: overwriting

	  if( MATRIX_COMPARE_ENTRIES( *dummy_matrix1_p, K_matrix_many_p->array[ index ], EPS2  ) ) info=1;

	  if( info ){

	    fprintf( stderr, "----------\n" );
	    fprintf( stderr, "ERROR in K_matrix_many check [The first one is the numerical estimate]\n" );
	    fprintf( stderr, "---> coor: (%d, %d)\n", i_coor1, i_coor2 );
	    fprintf( stderr, "----------\n" );
	  
	  }

	}  

      } /* i_coor2 loop */

    } /* i_coor1 loop */

  }


  if( info ){

    fprintf( stdout, "WARNING: Hamiltonian_many after a consistency error occurred\n\n" );

    fprintf( stdout, "H_matrix_many\n" );
    if( MATRIX_PRINT_PLUS( stdout, *H_matrix_many_p) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "F_matrix_many\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *F_matrix_many_p) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "K_matrix_many\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *K_matrix_many_p) ) info=1;
    fprintf( stdout, "\n" );

    fflush( stdout );

  }


  return info;

}

//------------------------------------------

#endif /* __DEBUG__ */

//------------------------------------------

/* check_particle_hole_matrix */

int check_particle_hole_matrix( matrix_p matrix_p, char* name ){

  /* dummies */
  matrix particle_hole;
  matrix dummy_matrix1;
  matrix dummy_matrix2;
  int dim=21;
  int info=0;


  /* allocation */
  if( MATRIX_ALLOCATE( dim, dim, particle_hole ) ) info=1;

  if( MATRIX_ALLOCATE( dim, dim, dummy_matrix1 ) ) info=1;

  if( MATRIX_ALLOCATE( dim, dim, dummy_matrix2 ) ) info=1;


  if( MATRIX_ZERO( particle_hole ) ) info=1;

  particle_hole.matrix[ MATRIX_ARG( 0, 0, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 1, 1, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 2, 3, dim, dim )].z[0]   =  1.0e0;
  particle_hole.matrix[ MATRIX_ARG( 3, 2, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 4, 4, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 5, 5, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 6, 6, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 7, 8, dim, dim )].z[0]   =  1.0e0;
  particle_hole.matrix[ MATRIX_ARG( 8, 7, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 9, 9, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 10, 10, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 11, 15, dim, dim )].z[0]   =  1.0e0;
  particle_hole.matrix[ MATRIX_ARG( 15, 11, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 12, 12, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 13, 16, dim, dim )].z[0]   =  1.0e0;
  particle_hole.matrix[ MATRIX_ARG( 16, 13, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 14, 17, dim, dim )].z[0]  =  1.0e0;
  particle_hole.matrix[ MATRIX_ARG( 17, 14, dim, dim )].z[0]  =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 18, 18, dim, dim )].z[0]   =  1.0e0;
 
  particle_hole.matrix[ MATRIX_ARG( 19, 19, dim, dim )].z[0]   =  1.0e0;

  particle_hole.matrix[ MATRIX_ARG( 20, 20, dim, dim )].z[0]   =  1.0e0;


  /*
    if( MATRIX_COPY( dummy_matrix2, particle_hole ) ) info=1;
    
    if( MATRIX_MATRIX_PRODUCT( particle_hole, dummy_matrix2, dummy_matrix1 ) ) info=1;
    
    if( MATRIX_UNIT( dummy_matrix2 ) ) info=1;
    
    
    if( MATRIX_COMPARE( dummy_matrix1, dummy_matrix2 ) ){
    
    fprintf( stderr, "ERROR: particle_hole trasform is not well defined\n");
    fflush( stderr );
    
    info=1;
    
    }
  */


  if( MATRIX_MATRIX_PRODUCT( particle_hole, *matrix_p, dummy_matrix1 ) )  info=1;

  if( MATRIX_MATRIX_PRODUCT( dummy_matrix1, particle_hole, dummy_matrix2  ) ) info=1;


  if( MATRIX_COMPARE( *matrix_p, dummy_matrix2 ) ){

    fprintf( stderr, "ERROR: matrix %s is not invariant with respect to particle_hole\n", name );
    fflush( stderr );

    fprintf( stdout, "ERROR: matrix %s is not invariant with respect to particle_hole\n", name );
    fflush( stdout );


    info=1;

  }


  if( info ){

    fprintf( stdout, "particle_hole\n" );
    if( MATRIX_PRINT_PLUS( stdout, particle_hole ) ) info=1;
    fprintf( stdout, "\n" );
    
    fprintf( stdout, "original matrix [no diagonal entries]\n" );
    if( MATRIX_PRINT_PLUS( stdout, *matrix_p ) ) info=1;
    fprintf( stdout, "\n" );
    
    fprintf( stdout, "transformed matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, dummy_matrix2 ) ) info=1;
    fprintf( stdout, "\n" );
    
  }


  /* deallocation */
  if( MATRIX_FREE( dummy_matrix2 ) ) info=1;

  if( MATRIX_FREE( dummy_matrix1 ) ) info=1;

  if( MATRIX_FREE( particle_hole ) ) info=1;


  return info;

}

//------------------------------------------
