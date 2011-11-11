
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
#include "PolyCEID_one_body_electronic_transition_matrix.h"


/*********************
  FUNCTIONS & MACROS
*********************/

//BUGFIX: this was originally DEF
int compute_one_body_electronic_transition_matrix( const constants constants, state_p state_p, config_p config_p ){

  /* constanst */
  int        N_levels_single;
  /* state */
  matrix_p   one_body_electronic_hole_matrix_p;
  matrix_p   one_body_electronic_hole_matrix_Ehrenfest_p;
  matrix_p   one_body_electronic_particle_matrix_p;
  matrix_p   one_body_electronic_particle_matrix_Ehrenfest_p;
  rvector_p  hole_energies_p;
  rvector_p  particle_energies_p;
  matrix_p   hole_orbitals_p;
  matrix_p   particle_orbitals_p;
  rvector_p  hole_populations_p;
  rvector_p  particle_populations_p;
  matrix_p   Ehrenfest_frame_p;
  matrix_p   dummy_matrix_single1_p;
  matrix_p   dummy_matrix_single2_p;
  // matrix_p   dummy_matrix_single3_p;
  vector_p   dummy_vector_single1_p;
  vector_p   dummy_vector_single2_p;
  /* dummies */
  complex    dummy;
  int        i, j;
  int        info=0;


  N_levels_single                                 =  constants.N_levels_single;

  one_body_electronic_hole_matrix_p               = &(state_p->one_body_electronic_hole_matrix);
  one_body_electronic_hole_matrix_Ehrenfest_p     = &(state_p->one_body_electronic_hole_matrix_Ehrenfest);
  one_body_electronic_particle_matrix_p           = &(state_p->one_body_electronic_particle_matrix);
  one_body_electronic_particle_matrix_Ehrenfest_p = &(state_p->one_body_electronic_particle_matrix_Ehrenfest);
  hole_orbitals_p                                 = &(state_p->hole_orbitals);
  particle_orbitals_p                             = &(state_p->particle_orbitals);
  hole_populations_p                              = &(state_p->hole_populations);
  particle_populations_p                          = &(state_p->particle_populations);
  hole_energies_p                                 = &(state_p->hole_energies);
  particle_energies_p                             = &(state_p->particle_energies);
  Ehrenfest_frame_p                               = &(config_p->electrons.Ehrenfest_frame);
  dummy_matrix_single1_p                          = &(state_p->dummy_matrix_single1);
  dummy_matrix_single2_p                          = &(state_p->dummy_matrix_single2);
  // dummy_matrix_single3_p                          = &(state_p->dummy_matrix_single3);
  dummy_vector_single1_p                          = &(state_p->dummy_vector_single1);
  dummy_vector_single2_p                          = &(state_p->dummy_vector_single2);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_one_body_electronic_hole_matrix\n" );

#endif /* __DEBUG_PLUS__ */


  if( COMPUTE_ONE_BODY_ELECTRONIC_TRANSITION_MATRIX_EHRENFEST( constants, *state_p ) ) info=1;

  /*
    fprintf( stdout, "one_body_electronic_hole_matrix_Ehrenfest\n" );
    if( MATRIX_PRINT_PLUS( stdout, *one_body_electronic_hole_matrix_Ehrenfest_p ) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "one_body_electronic_particle_matrix_Ehrenfest\n" );
    if( MATRIX_PRINT_PLUS( stdout, *one_body_electronic_particle_matrix_Ehrenfest_p ) ) info=1;
    fprintf( stdout, "\n" );
    
  */

  /*
    fprintf( stdout, "Ehrenfest_frame\n" );
    if( MATRIX_PRINT_PLUS( stdout, *Ehrenfest_frame_p ) ) info=1;
    fprintf( stdout, "\n" );
  */


  // hole
  if( MATRIX_TRANS( *Ehrenfest_frame_p, *dummy_matrix_single1_p ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( *one_body_electronic_hole_matrix_Ehrenfest_p, *dummy_matrix_single1_p, *dummy_matrix_single2_p ) ) info=1;

  if( MATRIX_CONJ( *Ehrenfest_frame_p, *dummy_matrix_single1_p ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_matrix_single2_p, *one_body_electronic_hole_matrix_p ) ) info=1;

  // diagonalisation (hole)
  if( DIAGONALISATION( *one_body_electronic_hole_matrix_p, *hole_orbitals_p, *hole_populations_p ) ) info=1;


  // particle
  if( MATRIX_TRANS( *Ehrenfest_frame_p, *dummy_matrix_single1_p ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( *one_body_electronic_particle_matrix_Ehrenfest_p, *dummy_matrix_single1_p, *dummy_matrix_single2_p ) ) info=1;

  if( MATRIX_CONJ( *Ehrenfest_frame_p, *dummy_matrix_single1_p ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_matrix_single2_p, *one_body_electronic_particle_matrix_p ) ) info=1;

  // diagonalisation (particle)
  if( DIAGONALISATION( *one_body_electronic_particle_matrix_p, *particle_orbitals_p, *particle_populations_p ) ) info=1;


  /*
    fprintf( stdout, "one_body_electronic_hole_matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, *one_body_electronic_hole_matrix_p ) ) info=1;
    fprintf( stdout, "\n" );
    
    fprintf( stdout, "one_body_electronic_particle_matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, *one_body_electronic_particle_matrix_p ) ) info=1;
    fprintf( stdout, "\n" );
  */


  // compute energies

  // single-body Hamiltonian
  if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, *dummy_matrix_single1_p ) ) info=1;
  
  
  /* temporary part */
  /*
    if( MATRIX_ADJOINT( *particle_orbitals_p, *dummy_matrix_single2_p ) ) info=1;
  
    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single2_p, *dummy_matrix_single1_p, *dummy_matrix_single3_p ) ) info=1;

    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single3_p, *hole_orbitals_p, *dummy_matrix_single2_p ) ) info=1;

    fprintf( stdout, "Transition_Hamiltonian\n" );
    if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
    fprintf( stdout, "\n" );
  */        
  /* end temporary part */

  
  for( i=0; i<N_levels_single; i++ ){

    // hole

    // dummy vector
    for( j=0; j< N_levels_single; j++ ){

      dummy_vector_single1_p->vector[ j ] = hole_orbitals_p->matrix[ ELECTRON_SINGLE_INDEX( j, i ) ];
	    
    } /* end j loop */

    if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single2_p ) ) info=1;

    dummy = SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single2_p );

    // fprintf( stderr, "[H] dummy = ( %le, %le) \n", dummy.z[0], dummy.z[1] );

    hole_energies_p->rvector[ i ] = REAL( dummy);


    // particle

    // dummy vector
    for( j=0; j< N_levels_single; j++ ){

      dummy_vector_single1_p->vector[ j ] = particle_orbitals_p->matrix[ ELECTRON_SINGLE_INDEX( j, i ) ];
	    
    } /* end j loop */

    if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single2_p ) ) info=1;

    dummy = SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single2_p );

    // fprintf( stderr, "[P] dummy = ( %le, %le) \n", dummy.z[0], dummy.z[1] );

    particle_energies_p->rvector[ i ] = REAL( dummy);


  } /* end i loop */
  

  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int compute_one_body_electronic_transition_matrix_Ehrenfest( const constants constants, state_p state_p ){

  /* constanst */
  unsigned short int  flag_no_spin_flip;
  int                 N_levels_single;
  int                 N_levels_many;
  imatrix             many_body_hybridisation_table;
  vector              symmetry_coefficient;
  /* state */
  matrix_p            adiabatic_states_many_p;
  matrix_p            electronic_density_eigenvectors_p;
  rvector_p           electronic_density_eigenvalues_p;
  matrix_p            one_body_electronic_hole_matrix_Ehrenfest_p;
  matrix_p            one_body_electronic_particle_matrix_Ehrenfest_p;
  vector_p            dummy_vector1_p;
  vector_p            dummy_vector2_p;
  matrix_p            dummy_matrix_single1_p;
  matrix_p            dummy_matrix_single2_p;
  matrix_p            dummy_matrix_single3_p;
  /* dummies */
  complex             dummy;
  complex             symmetry_coefficient_loc;
  int                 alpha;
  int                 level1, level2;
  int                 counter;
  int                 i, j;
  // int                 k;
  int                 index;
  int                 info=0;


  flag_no_spin_flip                                 =  constants.flag_no_spin_flip;
  N_levels_single                                   =  constants.N_levels_single;
  N_levels_many                                     =  constants.N_levels_many;
  many_body_hybridisation_table                     =  constants.many_body_hybridisation_table;
  symmetry_coefficient                              =  constants.symmetry_coefficient;

  adiabatic_states_many_p                           = &(state_p->adiabatic_states_many);
  electronic_density_eigenvalues_p                  = &(state_p->electronic_density_eigenvalues);
  electronic_density_eigenvectors_p                 = &(state_p->electronic_density_eigenvectors);
  one_body_electronic_hole_matrix_Ehrenfest_p       = &(state_p->one_body_electronic_hole_matrix_Ehrenfest);
  one_body_electronic_particle_matrix_Ehrenfest_p   = &(state_p->one_body_electronic_particle_matrix_Ehrenfest);
  dummy_vector1_p                                   = &(state_p->dummy_vector1);
  dummy_vector2_p                                   = &(state_p->dummy_vector2);
  dummy_matrix_single1_p                            = &(state_p->dummy_matrix_single1);
  dummy_matrix_single2_p                            = &(state_p->dummy_matrix_single2);
  dummy_matrix_single3_p                            = &(state_p->dummy_matrix_single3);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_one_body_electronic_hole_matrix_Ehrenfest\n" );

#endif /* __DEBUG_PLUS__ */


  /* set matrix to zero */
  if( MATRIX_ZERO( *one_body_electronic_hole_matrix_Ehrenfest_p ) )      info=1; // Important
  if( MATRIX_ZERO( *one_body_electronic_particle_matrix_Ehrenfest_p ) )  info=1; // Important


  /* initialise inistantaneous GS */
  for( i=0; i<N_levels_many; i++ ){

    dummy_vector1_p->vector[ i ] = adiabatic_states_many_p->matrix[ ELECTRON_MANY_INDEX( i, 0 ) ];

  }

  /*
    fprintf( stdout, "Adiabatic GS\n" );
    if( VECTOR_PRINT_PLUS( stdout, *dummy_vector1_p ) ) info=1;
    fprintf( stdout, "\n" );
  */

  // natural many-body state index
  for( alpha=0; alpha<N_levels_many; alpha++ ){


    /* initialise inistantaneous natural many-body state */
    for( i=0; i<N_levels_many; i++ ){
      
      dummy_vector2_p->vector[ i ] = electronic_density_eigenvectors_p->matrix[ ELECTRON_MANY_INDEX( i, alpha ) ];
      
    }

    /*
      fprintf( stdout, "Electroni_density_eigenvector [%d]\n", alpha );
      if( VECTOR_PRINT_PLUS( stdout, *dummy_vector2_p ) ) info=1;
      fprintf( stdout, "\n" );
    */

    /* set dummy matrix to zero */
    if( MATRIX_ZERO( *dummy_matrix_single1_p ) ) info=1; // Important


    /* set counter to zero */
    counter=0;
    
    // first many-body index
    for( i=0; i<N_levels_many; i++ ){


      // diagonal terms
      /* Diagonal terms no longer included according to R.L. Martin JCP 118, 4775 (2003) */

      // for( k=0; k<N_levels_single*SPIN_DEG; k++ ){
	
	// fprintf( stdout, "single_level %d, single_level_occupation %d\n", k, single_level_occupation.imatrix[i][ k ] );

	// if( single_level_occupation.imatrix[i][ k ] ){
	  
	  // level1 = k /SPIN_DEG;

	  // fprintf( stdout, "OK, level1 = %d\n", level1 );

	  // index = ELECTRON_SINGLE_INDEX( level1, level1 ); 


	  // dummy = CMPLX_PRODUCT( CONJ( dummy_vector2_p->vector[ i ] ), dummy_vector1_p->vector[ i ] );


	  /* updating */
	  // dummy_matrix_single1_p->matrix[ index ] = CMPLX_SUM( dummy_matrix_single1_p->matrix[ index ], dummy ); // WARNING: overwriting


	// } /* end conditional */

      // } /* end k loop */

      
      // second many-body index, off-diagonal entries only
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


	    /* updating */
	    index = ELECTRON_SINGLE_INDEX( level1, level2 );

	    dummy = CMPLX_PRODUCT( CONJ( dummy_vector2_p->vector[ i ] ), dummy_vector1_p->vector[ j ] );

	    dummy = CMPLX_PRODUCT( dummy, symmetry_coefficient_loc );                                              // WARNING: overwriting

	    dummy_matrix_single1_p->matrix[ index ] = CMPLX_SUM( dummy_matrix_single1_p->matrix[ index ], dummy ); // WARNING: overwriting


	    //
	    index = ELECTRON_SINGLE_INDEX( level2, level1 );                                                       // WARNING: reverse ordering

	    dummy = CMPLX_PRODUCT( CONJ( dummy_vector2_p->vector[ j ] ), dummy_vector1_p->vector[ i ] );           // WARNING: reverse ordering

	    dummy = CMPLX_PRODUCT( dummy, CONJ( symmetry_coefficient_loc ) );                                      // WARNING: overwriting 
	    // WARNING: in principle this CONJ is redondant

	    dummy_matrix_single1_p->matrix[ index ] = CMPLX_SUM( dummy_matrix_single1_p->matrix[ index ], dummy ); // WARNING: overwriting
	    
	  }

	} /* end conditional */

      } /* end j loop */

    } /* end i loop */


    /*
      fprintf( stdout, "Transition matrix (Martin) [%d]\n", alpha );
      if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
      fprintf( stdout, "\n" );
    */


    /* one_body_electronic_hole_matrix_Ehrenfest updating */

    // forming the "hole matrix" for a single natural many-body
    // state alpha.
    // This is the only covariant way to do it, equivalent to
    // SVD in the subsector given by the state alpha
    if( MATRIX_ADJOINT( *dummy_matrix_single1_p, *dummy_matrix_single2_p ) ) info=1;

    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single2_p, *dummy_matrix_single1_p, *dummy_matrix_single3_p ) ) info=1;
    

    // here multiple each entries of the dummy matrix
    // by the populalation of the natural orbital alpha
    // and add to "one_body_electronic_hole_matrix_Ehrenfest"
    for( i=0; i<N_levels_single; i++ ){

      for( j=0; j<N_levels_single; j++ ){

	index = ELECTRON_SINGLE_INDEX( i, j );


        dummy = CMPLX_PRODUCT( CMPLX( electronic_density_eigenvalues_p->rvector[ alpha ] ), dummy_matrix_single3_p->matrix[ index ] );

	one_body_electronic_hole_matrix_Ehrenfest_p->matrix[ index ] = CMPLX_SUM( one_body_electronic_hole_matrix_Ehrenfest_p->matrix[ index ], dummy ); // WARNING: overwriting

      } /* end j loop */

    } /* end i loop */
    


    /* one_body_electronic_particle_matrix_Ehrenfest updating */

    // forming the "particle matrix" for a single natural many-body
    // state alpha.
    // This is the only covariant way to do it, equivalent to
    // SVD in the subsector given by the state alpha
    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, *dummy_matrix_single2_p, *dummy_matrix_single3_p ) ) info=1; // WARNING:reversed order for particles!
    

    // here multiples each entries of the dummy matrix
    // by the populalation of the natural orbital alpha
    // and add to "one_body_electronic_particle_matrix_Ehrenfest"
    for( i=0; i<N_levels_single; i++ ){

      for( j=0; j<N_levels_single; j++ ){

	index = ELECTRON_SINGLE_INDEX( i, j );


        dummy = CMPLX_PRODUCT( CMPLX( electronic_density_eigenvalues_p->rvector[ alpha ] ), dummy_matrix_single3_p->matrix[ index ] );

	one_body_electronic_particle_matrix_Ehrenfest_p->matrix[ index ] = CMPLX_SUM( one_body_electronic_particle_matrix_Ehrenfest_p->matrix[ index ], dummy ); // WARNING: overwriting
	
      } /* end j loop */

    } /* end i loop */
    
    /*
      fprintf( stdout, "one_body_electronic_hole_matrix_Ehrenfest [%d]\n", alpha );
      if( MATRIX_PRINT_PLUS( stdout, *one_body_electronic_hole_matrix_Ehrenfest_p ) ) info=1;
      fprintf( stdout, "\n" );

      fprintf( stdout, "one_body_electronic_particle_matrix_Ehrenfest [%d]\n", alpha );
      if( MATRIX_PRINT_PLUS( stdout, *one_body_electronic_particle_matrix_Ehrenfest_p ) ) info=1;
      fprintf( stdout, "\n" );
    */

  } /* end loop alpha */


  if( flag_no_spin_flip ){

    if( MATRIX_SCALAR_DIVISION( *one_body_electronic_hole_matrix_Ehrenfest_p, CMPLX( (double) SPIN_DEG ), 
                                *one_body_electronic_hole_matrix_Ehrenfest_p) ) info=1;                         // WARNING: overwriting

    if( MATRIX_SCALAR_DIVISION( *one_body_electronic_particle_matrix_Ehrenfest_p, CMPLX( (double) SPIN_DEG ), 
                                *one_body_electronic_particle_matrix_Ehrenfest_p) ) info=1;                     // WARNING: overwriting

  }


  return info;

}

//------------------------------------------
