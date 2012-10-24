
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
#include "PolyCEID_nonadiabatic_coupling.h"


/*********************
  FUNCTIONS & MACROS
*********************/

//BUGFIX: this was originally TMP
int compute_nonadiabaticity( constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  rvector_p       nonadiabaticity_p;
  rvector_p       nonadiabatic_coupling_p;
  rvector_p       ap_p;
  matrix_p        H_matrix_p;
  rvector_p       dummy_rvector_p;
  matrix_p        dummy_matrix1_p;
  rvector_p       masses_aux_p;
  /* dummies */
  int             i_coor;
  int             index;
  int             i, j;
  double          delta_E;
  rvector         adiabatic_populations;
  matrix          nonadiabatic_forces;
  rvector         momenta_new;
  int             info=0;


  N_coor                   =  constants.N_coor;
  N_levels_many            =  constants.N_levels_many;
  masses_aux_p             = &config_p->atoms.masses_aux;

  nonadiabatic_coupling_p  = &(state_p->nonadiabatic_coupling);
  nonadiabaticity_p        = &(state_p->nonadiabaticity);
  ap_p                     = &(state_p->phonons.ap);
  H_matrix_p               = &(config_p->electrons.H_matrix);
  dummy_rvector_p          = &(state_p->dummy_rvector);
  dummy_matrix1_p          = &(state_p->dummy_matrix1);


  // initial allocations
  if( RVECTOR_ALLOCATE( N_levels_many, adiabatic_populations ) )              info=0;

  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, nonadiabatic_forces ) )  info=0;

  if( RVECTOR_ALLOCATE( N_coor, momenta_new ) )                               info=0;


  // transform momenta
  if( TRANSFORM_MOMENTA( constants, *state_p, *config_p, momenta_new ) )      info=0; 

    
  // if( RVECTOR_PRINT( stdout, momenta_new ) ) info=0;


  /* Hamiltonian diagonalisation */
  if( DIAGONALISATION( *H_matrix_p, *dummy_matrix1_p, *dummy_rvector_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "eigenvectors\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "eigenvalues\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *dummy_rvector_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  // compute adiabatic populations
  if( COMPUTE_ADIABATIC_POPULATIONS( constants, *state_p, *config_p, *dummy_matrix1_p, adiabatic_populations ) ) info=0;


  // loop on the (generalised) coordinates
  for( i_coor=0; i_coor<N_coor; i_coor++ ){

    // fprintf( stdout, "---> i_coor = %d\n", i_coor );

    // setting to zero
    nonadiabaticity_p->rvector[ i_coor ] = ZERO;

    if( masses_aux_p->rvector[ i_coor ] > EPS ){
    
      // compute forces and potentials
      if( COMPUTE_NONADIABATIC_FORCES( constants, *state_p, *config_p, *dummy_matrix1_p, i_coor, nonadiabatic_forces ) ) info=0;

      // Strategy here is: "sum on the final states and average over the initial ones"
      // loop on the possible many-body electronic inital states
      for( i=0; i<N_levels_many; i++ ){  // loop over the initial states

        // loop on the possible final many-body electronic states
        for( j=0; j<N_levels_many; j++ ){ // loop over the final sattes

          // The vectorised index
	  index = i_coor *N_levels_many *N_levels_many +i *N_levels_many +j; 

          // energy difference [adiabatic]
          delta_E = ( dummy_rvector_p->rvector[ j ] ) -( dummy_rvector_p->rvector[ i ] );

          // BUGFIX: frustrated hops excluded. Is that what we want?
	  if( fabs( delta_E ) > EPS ){

            // BUGFIX: hbar missing?
            nonadiabatic_coupling_p->rvector[ index ] = HBAR *CMPLX_NORM( nonadiabatic_forces.matrix[ ELECTRON_MANY_INDEX( i, j ) ] ) /fabs( delta_E ); // WARNING: overwriting

	  }
	  else{

            nonadiabatic_coupling_p->rvector[ index ] = ZERO;

	  }        

          // BUGFIX: frustrated hops not included
	  if( adiabatic_populations.rvector[ i ] > EPS && delta_E < ONEO2 *( momenta_new.rvector[ i_coor ] ) *( momenta_new.rvector[ i_coor ] ) /( masses_aux_p->rvector[ i_coor ] ) ){

            // updating the sum (coupling)
            nonadiabaticity_p->rvector[ i_coor ] += nonadiabatic_coupling_p->rvector[ index ] *adiabatic_populations.rvector[ i ]; // WARNING: overwriting
              
	  }        

        } /* end j loop */

      } /* end i loop */

      // giving the right dimensions [eV]
      // OSQRT2 *( ap_p->rvector[ i_coor ] ) is the zero-point momentum for mode i_coor
      nonadiabaticity_p->rvector[ i_coor ] *= TWO /HBAR *sqrt( ONEO2 *( ( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] ) ) +( momenta_new.rvector[ i_coor ] ) *( momenta_new.rvector[ i_coor ] ) ) /( masses_aux_p->rvector[ i_coor ] );

    }

  } /* end i_coor loop */


  // final deallocations
  if( RVECTOR_FREE( momenta_new ) )            info=0;

  if( MATRIX_FREE( nonadiabatic_forces ) )     info=0;

  if( RVECTOR_FREE( adiabatic_populations ) )  info=0;


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int compute_adiabatic_populations( const constants constants, state_p state_p, config_p config_p, matrix_p adiabatic_states_p, rvector_p adiabatic_populations_p ){

  /* constants */
  int             N_levels_many;
  /* state */
  matrix_p        mu00_p;
  vector_p        dummy_vector1_p;
  vector_p        dummy_vector2_p;
  /* dummies */
  int             i, k;
  int             info=0;


  N_levels_many            =  constants.N_levels_many;

  mu00_p                   = &(config_p->electrons.mu00);
  dummy_vector1_p          = &(state_p->dummy_vector1);
  dummy_vector2_p          = &(state_p->dummy_vector2);


  for( i=0; i<N_levels_many; i++ ){

    // copy the electronic many-body states to tmp vectors
    for( k=0; k<N_levels_many; k++ ){

      dummy_vector1_p->vector[ k ] = adiabatic_states_p->matrix[ ELECTRON_MANY_INDEX( k, i ) ];

    }

  
    //  computing the population, < i | mu00 | i >, in two steps
    if( MATRIX_VECTOR_PRODUCT( *mu00_p, *dummy_vector1_p, *dummy_vector2_p ) ) info=1;

    adiabatic_populations_p->rvector[ i ] = CMPLX_NORM( SCALAR_PRODUCT( *dummy_vector1_p, *dummy_vector2_p ) );

  }


  return info;

}

//------------------------------------------

int compute_nonadiabatic_forces( const constants constants, state_p state_p, config_p config_p, matrix_p adiabatic_states_p, int mode, matrix_p nonadiabatic_forces_p ){


  /* constants */
  int             N_levels_many;
  /* state */
  matrix_array_p  delta_F_matrix_p;
  vector_p        dummy_vector1_p;
  vector_p        dummy_vector2_p;
  vector_p        dummy_vector3_p;
  /* dummies */
  int             i, j, k;
  /* dummies */
  int info=0;


  N_levels_many            =  constants.N_levels_many;

  delta_F_matrix_p         = &(config_p->electrons.delta_F_matrix);
  dummy_vector1_p          = &(state_p->dummy_vector1);
  dummy_vector2_p          = &(state_p->dummy_vector2);
  dummy_vector3_p          = &(state_p->dummy_vector3);


  for( i=0; i<N_levels_many; i++ ){

    for( j=0; j<N_levels_many; j++ ){

      // copy the electronic many-body states to tmp vectors
      for( k=0; k<N_levels_many; k++ ){

        dummy_vector1_p->vector[ k ] = adiabatic_states_p->matrix[ ELECTRON_MANY_INDEX( k, i ) ];

      }

      for( k=0; k<N_levels_many; k++ ){

        dummy_vector2_p->vector[ k ] = adiabatic_states_p->matrix[ ELECTRON_MANY_INDEX( k, j ) ];

      }

  
      //  computing the population, < i | F | j >, in three steps
      //  WARNING: use delta_F because is already rotated!
      if( MATRIX_VECTOR_PRODUCT( delta_F_matrix_p->array[ mode ], *dummy_vector2_p, *dummy_vector3_p ) ) info=1;

      nonadiabatic_forces_p->matrix[ ELECTRON_MANY_INDEX( i, j ) ] = SCALAR_PRODUCT( *dummy_vector1_p, *dummy_vector3_p );

    }

  }


  return info;

}

//------------------------------------------

