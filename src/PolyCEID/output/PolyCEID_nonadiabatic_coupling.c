
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
  rvector_p       nonadiabatic_coupling_p;
  rvector_p       ap_p;
  matrix_p        H_matrix_p;
  rvector_p       dummy_rvector_p;
  matrix_p        dummy_matrix1_p;
  rvector_p       masses_aux_p;
  /* dummies */
  int             i_coor;
  int             i, j;
  double          delta_E;
  double          matrix_element_loc;
  rvector         adiabatic_populations;
  matrix          nonadiabatic_forces;
  rvector         momenta_new;
  int             info=0;


  N_coor                   =  constants.N_coor;
  N_levels_many            =  constants.N_levels_many;
  masses_aux_p             = &config_p->atoms.masses_aux;

  nonadiabatic_coupling_p  = &(state_p->nonadiabatic_coupling);
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
    nonadiabatic_coupling_p->rvector[ i_coor ] = ZERO;

    if( masses_aux_p->rvector[ i_coor ] > EPS ){
    
      // compute forces and potentials
      if( COMPUTE_NONADIABATIC_FORCES( constants, *state_p, *config_p, *dummy_matrix1_p, i_coor, nonadiabatic_forces ) )    info=0;

    
      // Strategy here is: "sum on the final states and average over the initial ones"
      // loop on the possible many-body electronic inital states
      for( i=0; i<N_levels_many; i++ ){  // loop over the initial states


        matrix_element_loc = ZERO;


        if( adiabatic_populations.rvector[ i ] > EPS ){

          // loop on the possible final many-body electronic states
          for( j=0; j<N_levels_many; j++ ){ // loop over the final sattes

            // energy difference [adiabatic]
            delta_E = ( dummy_rvector_p->rvector[ j ] ) -( dummy_rvector_p->rvector[ i ] );

            // BUGFIX: conical intersections: delta_E = 0 although i != j
	    // frustrated hopping condition added
            if( j !=i && masses_aux_p->rvector[ i_coor ] > EPS && delta_E < ONEO2 *( momenta_new.rvector[ i_coor ] ) *( momenta_new.rvector[ i_coor ] ) /( masses_aux_p->rvector[ i_coor ] ) ){

              matrix_element_loc  += CMPLX_NORM( nonadiabatic_forces.matrix[ ELECTRON_MANY_INDEX( i, j ) ] ) /fabs( delta_E ); // WARNING: overwriting

	    }     

          } /* end j loop */


          // updating the sum (coupling)
          nonadiabatic_coupling_p->rvector[ i_coor ] += matrix_element_loc *adiabatic_populations.rvector[ i ]; // WARNING: overwriting

        }

      } /* end i loop */

      // giving the right dimensions [s^{-1}]
      // OSQRT2 *( ap_p->rvector[ i_coor ] ) is the zero-point momentum for mode i_coor
      nonadiabatic_coupling_p->rvector[ i_coor ] *= TWO *sqrt( ONEO2 *( ( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] ) ) +( momenta_new.rvector[ i_coor ] ) *( momenta_new.rvector[ i_coor ] ) ) /( masses_aux_p->rvector[ i_coor ] );

    }

  } /* end i_coor loop */


  // final deallocations
  if( RVECTOR_FREE( momenta_new ) )            info=0;

  if( MATRIX_FREE( nonadiabatic_forces ) )     info=0;

  if( RVECTOR_FREE( adiabatic_populations ) )  info=0;


  return info;

}


//BUGFIX: this was originally TMP
int compute_nonadiabatic_rate( constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  rvector_p       omega_p;
  rvector_p       masses_aux_p;
  rvector_p       nonadiabatic_rate_p;
  matrix_p        H_matrix_p;
  rvector_p       dummy_rvector_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  int             i_coor;
  int             i, j;
  double          delta_E, delta_G;
  double          lambda_i, lambda_j;
  double          matrix_element_loc;
  double          temperature;
  rvector         adiabatic_populations;
  matrix          nonadiabatic_forces;
  rvector         adiabatic_potentials;
  rvector         momenta_new;
  rvector         forces_new;
  int             info=0;


  N_coor                   =  constants.N_coor;
  N_levels_many            =  constants.N_levels_many;
  masses_aux_p             = &config_p->atoms.masses_aux;

  omega_p                  = &(state_p->phonons.omega);
  nonadiabatic_rate_p      = &(state_p->nonadiabatic_rate);
  H_matrix_p               = &(config_p->electrons.H_matrix);
  dummy_rvector_p          = &(state_p->dummy_rvector);
  dummy_matrix1_p          = &(state_p->dummy_matrix1);


  // initial allocations
  if( RVECTOR_ALLOCATE( N_levels_many, adiabatic_populations ) )              info=0;

  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, nonadiabatic_forces ) )  info=0;

  if( RVECTOR_ALLOCATE( N_levels_many, adiabatic_potentials ) )               info=0;

  if( RVECTOR_ALLOCATE( N_coor, momenta_new ) )                               info=0;

  if( RVECTOR_ALLOCATE( N_coor, forces_new ) )                                info=0;


  // transform momenta
  if( TRANSFORM_MOMENTA( constants, *state_p, *config_p, momenta_new ) )      info=0; 

  // transform forces
  if( TRANSFORM_FORCES( constants, *state_p, *config_p, forces_new ) )        info=0; 

  /*  
    if( RVECTOR_PRINT( stdout, momenta_new ) ) info=0;
   
    if( RVECTOR_PRINT( stdout, forces_new ) ) info=0;
  */ 

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

    /*
      fprintf( stdout, "------------------------------------\n" );
      fprintf( stdout, "---> i_coor = %d\n", i_coor );
    */

    // setting to zero
    nonadiabatic_rate_p->rvector[ i_coor ] = ZERO;

    if( masses_aux_p->rvector[ i_coor ] > EPS ){
    
      // compute forces and potentials
      if( COMPUTE_NONADIABATIC_FORCES( constants, *state_p, *config_p, *dummy_matrix1_p, i_coor, nonadiabatic_forces ) )    info=0;

      if( COMPUTE_ADIABATIC_POTENTIALS( constants, *state_p, *config_p, *dummy_matrix1_p, i_coor, adiabatic_potentials ) )  info=0;
    

      // Strategy here is: "sum on the final states and average over the initial ones"
      // loop on the possible many-body electronic inital states
      for( i=0; i<N_levels_many; i++ ){  // loop over the initial states


        // reorganisation energy [i]
        lambda_i  = REAL( nonadiabatic_forces.matrix[ ELECTRON_MANY_INDEX( i, i ) ] ) +forces_new.rvector[ i_coor ]; 
        lambda_i *= lambda_i *ONEO2 /( adiabatic_potentials.rvector[ i ] ); //WARNING: overwriting


        // effective temperature
        // WARNING: masses_aux is for the normal modes
        temperature  = ONEO2 *( momenta_new.rvector[ i_coor ] ) *( momenta_new.rvector[ i_coor ] ) /( masses_aux_p->rvector[ i_coor ] );

        temperature += lambda_i; // adding reorganisation energy


        // loop on the possible final many-body electronic states
        for( j=0; j<N_levels_many; j++ ){ // loop over the final sattes


          // energy difference [adiabatic]
          delta_E = ( dummy_rvector_p->rvector[ j ] ) -( dummy_rvector_p->rvector[ i ] );

          // reorganisation energy [j]
          lambda_j  = REAL( nonadiabatic_forces.matrix[ ELECTRON_MANY_INDEX( j, j ) ] ) +forces_new.rvector[ i_coor ]; 
          lambda_j *= lambda_j *ONEO2 /( adiabatic_potentials.rvector[ j ] ); //WARNING: overwriting


          // energy difference [minima]
	  delta_G = delta_E -( lambda_j - lambda_i );


          // we want emission, only
	  if( delta_G < -EPS && adiabatic_populations.rvector[ i ] > EPS ){

          
            /*
              fprintf( stdout, "%le\n", REAL( nonadiabatic_forces.matrix[ ELECTRON_MANY_INDEX( i, i ) ] ) );
              fprintf( stdout, "%le\n", REAL( nonadiabatic_forces.matrix[ ELECTRON_MANY_INDEX( j, j ) ] ) );
              fprintf( stdout, "\n" );
              fprintf( stdout, "%le\n", forces_new.rvector[ i_coor ] );
              fprintf( stdout, "\n" );
              fprintf( stdout, "%le\n", adiabatic_potentials.rvector[ i ] );
              fprintf( stdout, "%le\n", adiabatic_potentials.rvector[ j ] );
              fprintf( stdout, "\n" );
      
              fprintf( stdout, "i = %d, j = %d\n", i, j );
              fprintf( stdout, "\n" );
              fprintf( stdout, "lambda[%d] = %le, lambda[%d] = %le\n", i, lambda_i, j, lambda_j );
              fprintf( stdout, "delta_E    = %le, delta_G    = %le\n", delta_E, delta_G );
              fprintf( stdout, "temperature[%d]              = %le\n", i, temperature );
              fprintf( stdout, "adiabatic_populations[%d]    = %le\n", i, adiabatic_populations.rvector[ i ] );
              fprintf( stdout, "\n" );
            */

            // norm of the matrix element
            matrix_element_loc  = CMPLX_NORM( nonadiabatic_forces.matrix[ ELECTRON_MANY_INDEX( i, j ) ] );

            // squaring
            matrix_element_loc *= matrix_element_loc;

	    // averaging
	    matrix_element_loc *= adiabatic_populations.rvector[ i ]; // WARNING: overwriting


            // fprintf( stdout, "matrix_element_loc    [%d,%d]  = %le\n", i, j, matrix_element_loc );


            // including Frank-Condon factor
            matrix_element_loc *= FC_FACTOR( delta_G, HBAR *( omega_p->rvector[ i_coor ] ), lambda_i, temperature );

            /*
              fprintf( stdout, "hw_ph  = %le\n", HBAR *( omega_p->rvector[ i_coor ] ) );
              fprintf( stdout, "matrix_element_loc+FC [%d,%d]  = %le\n", i, j, matrix_element_loc );
              fprintf( stdout, "\n" );
            */

	    // updating the sum (rate)
	    nonadiabatic_rate_p->rvector[ i_coor ] += matrix_element_loc; // WARNING: overwriting


	  } /* end delta_E conditional */

        } /* end j loop */

      } /* end i loop */

      // giving the right dimensions [s^{-1}]
      nonadiabatic_rate_p->rvector[ i_coor ]  *= TWO *PI /HBAR; // WARNING: overwriting

    }

  } /* end i_coor loop */


  // final deallocations
  if( RVECTOR_FREE( forces_new ) )             info=0;

  if( RVECTOR_FREE( momenta_new ) )            info=0;

  if( RVECTOR_FREE( adiabatic_potentials ) )   info=0;

  if( MATRIX_FREE( nonadiabatic_forces ) )     info=0;

  if( RVECTOR_FREE( adiabatic_populations ) )  info=0;


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

double line_shape( double energy, double denergy ){

  /* dummies */
  double  dummy;
  double  value;


  // fprintf( stdout, "energy = %le, denergy = %le\n", energy, denergy );


  dummy = energy /denergy;  

  // Lorentzian
  // value = 1.0e0 /( PI *denergy *( dummy*dummy + 1.0e0 ) );
  
  // Gaussian
  value = exp( -ONEO2 *dummy *dummy ) /( sqrt( TWO *PI ) *denergy );


  return value;

}	

//------------------------------------------

double FC_factor( double energy, double ph_energy, double lambda, double temperature ){

  /* dummies */
  double  nph;
  double  dummy;
  double  HR_factor;
  double  value;


  dummy = -( energy +lambda ) /ph_energy ;

  // nearest integer
  nph   = floor( dummy );
  if( dummy -nph > ONEO2 ){

     nph += ONE;

  }        

  // non-negative number of phonons
  if( nph < ZERO ) nph=ZERO;

  
  // fprintf( stdout, "nph = %le\n", nph );


  // Huang-Rhys factor
  HR_factor = lambda /ph_energy;


  // fprintf( stdout, "HR_factor = %le\n", HR_factor );


  // max of the Poisson distribution
  if( nph < 10 ){

    value = pow( HR_factor, nph ) *exp( -HR_factor ) /FACTORIAL( (int) nph );

  }
  else{

    value = ONE /sqrt(TWO *PI *HR_factor); // Stirling's approximation

  }


  // fprintf( stdout, "value    = %le\n", value );


  // adding line shape
  value *= LINE_SHAPE( energy +lambda +nph *ph_energy, sqrt( TWO *lambda *temperature ) );

  /*
    fprintf( stdout, "value+LS = %le\n", value );
    fprintf( stdout, "\n", value );
  */


  return value;

}	

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

int compute_adiabatic_density_matrix( const constants constants, state_p state_p, config_p config_p, matrix_p adiabatic_states_p, matrix_p adiabatic_density_matrix_p ){

  /* constants */
  int             N_levels_many;
  /* state */
  matrix_p        mu00_p;
  vector_p        dummy_vector1_p;
  vector_p        dummy_vector2_p;
  vector_p        dummy_vector3_p;
  /* dummies */
  int             i, j, k;
  int             info=0;


  N_levels_many            =  constants.N_levels_many;

  mu00_p                   = &(config_p->electrons.mu00);
  dummy_vector1_p          = &(state_p->dummy_vector1);
  dummy_vector2_p          = &(state_p->dummy_vector2);
  dummy_vector3_p          = &(state_p->dummy_vector3);


  for( i=0; i<N_levels_many; i++ ){

    // copy the electronic many-body states to tmp vectors
    for( k=0; k<N_levels_many; k++ ){

      dummy_vector1_p->vector[ k ] = adiabatic_states_p->matrix[ ELECTRON_MANY_INDEX( k, i ) ];

    }

    for( j=0; j<N_levels_many; j++ ){

      // copy the electronic many-body states to tmp vectors
      for( k=0; k<N_levels_many; k++ ){

        dummy_vector2_p->vector[ k ] = adiabatic_states_p->matrix[ ELECTRON_MANY_INDEX( k, j ) ];

      }

  
      //  computing the population, < i | mu00 | j >, in two steps
      if( MATRIX_VECTOR_PRODUCT( *mu00_p, *dummy_vector2_p, *dummy_vector3_p ) ) info=1;

      adiabatic_density_matrix_p->matrix[ ELECTRON_MANY_INDEX( i, j) ] = SCALAR_PRODUCT( *dummy_vector1_p, *dummy_vector3_p );

    }

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

// WARNING: missing fan term!
int compute_adiabatic_potentials( const constants constants, state_p state_p, config_p config_p, matrix_p adiabatic_states_p, int mode, rvector_p adiabatic_potentials_p ){


  /* constants */
  int             N_levels_many;
  /* state */
  matrix_array_p  K_matrix_p;
  vector_p        dummy_vector1_p;
  vector_p        dummy_vector2_p;
  /* dummies */
  int             i, k;
  /* dummies */
  int info=0;


  N_levels_many            =  constants.N_levels_many;

  K_matrix_p               = &(config_p->electrons.K_matrix);
  dummy_vector1_p          = &(state_p->dummy_vector1);
  dummy_vector2_p          = &(state_p->dummy_vector2);


  for( i=0; i<N_levels_many; i++ ){

    // copy the electronic many-body states to tmp vectors
    for( k=0; k<N_levels_many; k++ ){

      dummy_vector1_p->vector[ k ] = adiabatic_states_p->matrix[ ELECTRON_MANY_INDEX( k, i ) ];

    }

  
    //  computing the population, < i | K | i >, in two steps
    if( MATRIX_VECTOR_PRODUCT( K_matrix_p->array[ COORDINATE_INDEX(mode, mode) ], *dummy_vector1_p, *dummy_vector2_p ) ) info=1;

    adiabatic_potentials_p->rvector[ i ] = CMPLX_NORM( SCALAR_PRODUCT( *dummy_vector1_p, *dummy_vector2_p ) );

  }

  return info;

}

//------------------------------------------

