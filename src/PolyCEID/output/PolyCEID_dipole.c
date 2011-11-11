
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
// Sia lode alla Vergine Immacolata!
#include "PolyCEID_dipole.h"


/*********************
  FUNCTIONS & MACROS
*********************/


/* compute superpositions */

int compute_superpositions( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             max_rho_index;
  /* state */
  matrix_array_p  rho_p;
  matrix_array_p  initial_state_recorded_p;
  matrix_array_p  excited_state_recorded_p;
  double*         superposition_instantaneous_and_initial_state_p;
  double*         superposition_instantaneous_and_excited_state_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  complex         dummy;
  int             index, index_aux;
  int             info=0;
  

  max_rho_index                                   =  constants.max_rho_index;
 
  rho_p                                            = &(config_p->electrons.rho);
  initial_state_recorded_p                        = &(state_p->initial_state_recorded);
  excited_state_recorded_p                        = &(state_p->excited_state_recorded);
  superposition_instantaneous_and_initial_state_p = &(state_p->observables.superposition_instantaneous_and_initial_state);
  superposition_instantaneous_and_excited_state_p = &(state_p->observables.superposition_instantaneous_and_excited_state);
  dummy_matrix1_p                                 = &(state_p->dummy_matrix1);


  /* set to zero */
  *superposition_instantaneous_and_initial_state_p = 0.0e0;
  *superposition_instantaneous_and_excited_state_p = 0.0e0;

  for( index=0; index<max_rho_index; index++ ){

    index_aux = rho_index_conjugate.ivector[ index ];


    /*  superposition_instantaneous_and_initial_state */
    if( MATRIX_MATRIX_PRODUCT( rho_p->array[ index_aux ], initial_state_recorded_p->array[ index ], *dummy_matrix1_p ) ) info=1;

    dummy = MATRIX_TRACE( *dummy_matrix1_p );


    *superposition_instantaneous_and_initial_state_p += REAL( dummy );

    /*
    fprintf( stdout, "rho\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *rho_p ) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "excited_state_recorded\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *excited_state_recorded_p ) ) info=1;
    fprintf( stdout, "\n" );
    */

    /*  superposition_instantaneous_and_excited_state */
    if( MATRIX_MATRIX_PRODUCT( rho_p->array[ index_aux ], excited_state_recorded_p->array[ index ], *dummy_matrix1_p ) ) info=1;

    dummy = MATRIX_TRACE( *dummy_matrix1_p );


    *superposition_instantaneous_and_excited_state_p += REAL( dummy );


  } /* end index loop */


  return info;

}

//------------------------------------------

/* compute classical_dipole */

int compute_classical_dipole( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        N_levels_many;
  int        sdim;
  /* state */
  matrix_p   dipole_many_p;
  rvector_p  oscillator_strength_many_p;
  rvector_p  oscillator_frequency_many_p;
  rvector_p  adiabatic_PES_many_p;
  matrix_p   adiabatic_states_many_p;
  matrix_p   dummy_matrix1_p;
  matrix_p   dummy_matrix2_p;
  /* dummies */
  int        i,j;
  int        index;
  int        i_comp;
  double     norm2;
  int        info=0;
  

  N_levels_many               =  constants.N_levels_many;
  sdim                        =  constants.spacial_dimension; 

  dipole_many_p               = &config_p->electrons.dipole_many;
  oscillator_strength_many_p  = &state_p->observables.oscillator_strength_many;
  oscillator_frequency_many_p = &state_p->observables.oscillator_frequency_many;
  adiabatic_PES_many_p        = &state_p->adiabatic_PES_many;
  adiabatic_states_many_p     = &state_p->adiabatic_states_many;
  dummy_matrix1_p             = &state_p->dummy_matrix1;
  dummy_matrix2_p             = &state_p->dummy_matrix2;


  for( i_comp=0; i_comp<sdim; i_comp++ ){

    // Change to adiabatic reference frame
    if( MATRIX_ADJOINT( *adiabatic_states_many_p, *dummy_matrix1_p ) );

    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, *dipole_many_p, *dummy_matrix2_p ) )         info=1;

    if( MATRIX_MATRIX_PRODUCT(  *dummy_matrix2_p, *adiabatic_states_many_p, *dummy_matrix1_p ) )  info=1;

#ifdef __DEBUG_PLUS__

    fprintf( stdout, "dipole_many\n" );
    if( MATRIX_PRINT_PLUS( stdout, *dipole_many_p ) )  info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "dipole_many [adiabatic frame]\n" );
    if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p ) ) info=1;
    fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */

    // oscillator strength and frequencies

    index = i_comp *N_levels_many *(N_levels_many -1) /2;

    for( i=0; i<N_levels_many; i++ ){

      for( j=(i+1); j<N_levels_many; j++ ){
            
         // oscillator frequency
         oscillator_frequency_many_p->rvector[ index ] = adiabatic_PES_many_p->rvector[ j ] - 
                                                         adiabatic_PES_many_p->rvector[ i ];

         // oscillator strength
         norm2 = CMPLX_NORM( dummy_matrix1_p->matrix[ ELECTRON_MANY_INDEX( i, j ) ] );

         norm2 = norm2 *norm2;

         oscillator_strength_many_p->rvector[ index ]  = TWO /sdim *E_MASS_OVER_AMU *AMU_TO_INTERNAL /( HBAR *HBAR ) 
	                                               *( oscillator_frequency_many_p->rvector[ index ] ) *norm2;
         // WARNING: not multiplied by SPIN_DEG because I want emission (singlet ground-state assumed here)
         // I don't think it must be divided by the numeber of electrons, 
	 // I always used single particle operators [not sure]
	 // f-sum rule not valid, beause of the discretised Hamiltonian

         index++;

      } /* end j loop */        

    } /* end i loop */        

  } /* end i_comp loop */

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "oscillator_strength_many\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *oscillator_strength_many_p ) )  info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "oscillator_frequency_many\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *oscillator_frequency_many_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------
