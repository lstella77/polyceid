
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
#include "PolyCEID_dipole_single.h"


/*********************
  FUNCTIONS & MACROS
*********************/


/* compute_dipole_single */
//BUGFIX: this was originally DEF
int compute_dipole_single( const constants constants, state_p state_p, config_p config_p ){

  /* constants*/
  int        N_levels_single;
  int        sdim;
  /* state */
  matrix_p   hole_orbitals_p;
  matrix_p   particle_orbitals_p;
  rvector_p  hole_populations_p;
  rvector_p  particle_populations_p;
  rvector_p  hole_energies_p;
  rvector_p  particle_energies_p;
  rvector_p  oscillator_strength_single_p;
  rvector_p  oscillator_frequency_single_p;
  matrix_p   dummy_matrix_single1_p;
  vector_p   dummy_vector_single1_p;
  vector_p   dummy_vector_single2_p;
  /* dummies */
  complex    dummy;
  double     norm2;
  int        i;
  int        j;
  int        comp;
  int        info=0;


  N_levels_single                = constants.N_levels_single;
  sdim                           =  constants.spacial_dimension; 

  hole_orbitals_p                = &(state_p->hole_orbitals);
  particle_orbitals_p            = &(state_p->particle_orbitals);
  hole_populations_p             = &(state_p->hole_populations);
  particle_populations_p         = &(state_p->particle_populations);
  hole_energies_p                = &(state_p->hole_energies);
  particle_energies_p            = &(state_p->particle_energies);
  oscillator_strength_single_p   = &(state_p->observables.oscillator_strength_single);
  oscillator_frequency_single_p  = &(state_p->observables.oscillator_frequency_single);
  dummy_matrix_single1_p         = &(state_p->dummy_matrix_single1);  
  dummy_vector_single1_p         = &(state_p->dummy_vector_single1); 
  dummy_vector_single2_p         = &(state_p->dummy_vector_single2); 
  

   if( strcmp( constants.hamiltonian.class, "CHAIN_MOD2" ) &&
       strcmp( constants.hamiltonian.class, "SSH" ) ){
     
     fprintf( stderr, "# WARNING: dipole perturbation for this model_hamiltonian is not supported yet [%s]\n", constants.hamiltonian.class );
     fflush( stderr );
 
   }
   else{
	   
     // compute single-body dipole operator
	     
    /* compute centre_of_mass */
    if( COMPUTE_CENTRE_OF_MASS( constants, config_p->atoms ) ) info=1;
		     
    /* classical_dipole_updating */
    comp=0; //BUGFIX: here there should be the loop over cartesian components...

    if( CLASSICAL_DIPOLE_SINGLE_UPDATE( constants, *state_p, *config_p, comp, *dummy_matrix_single1_p ) ) info=1;

    /*
      fprintf( stdout, "Dipole operator (single)\n" );
      if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
      fprintf( stdout, "\n" );
    */

    for( i=0; i<N_levels_single; i++ ){
	  
      // check
      if( fabs( hole_populations_p->rvector[ i ] - particle_populations_p->rvector[ i ] ) > EPS ){

        fprintf( stderr, "ERROR: hole and particle populations differ. [%le, %le]\n",  
            hole_populations_p->rvector[ i ], particle_populations_p->rvector[ i ] );
        fflush( stderr );

        info=1;

        break;

      }
    
      // check if the population whether the population is negligeable
      if( hole_populations_p->rvector[ i ] >  10 *EPS ){

        // compute oscillator frequency
        oscillator_frequency_single_p->rvector[ i ] = particle_energies_p->rvector[ i ] - hole_energies_p->rvector[ i ];
    
        // check
        //if( oscillator_frequency_single_p->rvector[ i ] < -EPS ){

          //fprintf( stderr, "ERROR: transition energy is negative. [%le]\n", oscillator_frequency_single_p->rvector[ i ] );      
          //fflush( stderr );

          //info=1;

          //break;
    
        //}
      
        // compute weighted oscillator strength
        
        // copy hole
        for( j=0; j<N_levels_single; j++ ){

          dummy_vector_single1_p->vector[ j ] = hole_orbitals_p->matrix[ ELECTRON_SINGLE_INDEX( j, i ) ]; 

        } /* end j loop */        


        /*
          fprintf( stdout, "hole_orbital %d\n", i );
          if( VECTOR_PRINT_PLUS( stdout, *dummy_vector_single1_p ) ) info=1;
          fprintf( stdout, "\n" );
        */


       if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single2_p ) ) info=1;


        // copy particle
        for( j=0; j<N_levels_single; j++ ){

          dummy_vector_single1_p->vector[ j ] = particle_orbitals_p->matrix[ ELECTRON_SINGLE_INDEX( j, i ) ];

        } /* end j loop */


        /*
          fprintf( stdout, "particle_orbital %d\n", i );
          if( VECTOR_PRINT_PLUS( stdout, *dummy_vector_single1_p ) ) info=1;
          fprintf( stdout, "\n" );
        */


        // matrix element
        dummy = SCALAR_PRODUCT( *dummy_vector_single1_p, *dummy_vector_single2_p );


        /*
          fprintf( stdout, "dummy [1]\n" );
          if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
          fprintf( stdout, "\n" );
        */


        // expectation value
        norm2 = CMPLX_NORM( dummy );

        norm2 = norm2 *norm2;
        norm2 *= 2.0;
	// The factor 2 is due to the degeneracy (only single-particle states)

      
        /*
          fprintf( stdout, "dummy [2]\n" );
          if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
          fprintf( stdout, "\n" );
      
          fprintf( stdout, "hole_populations [%d] = %le\n", i, hole_populations_p->rvector[ i ] );
          fprintf( stdout, "oscillator_frequency [%d] = %le\n", i, oscillator_frequency_single_p->rvector[ i ] );
        */


        // waighted oscillator strength
        oscillator_strength_single_p->rvector[ i ] = TWO /sdim 
          *E_MASS_OVER_AMU *AMU_TO_INTERNAL /( HBAR *HBAR ) *( oscillator_frequency_single_p->rvector[ i ] ) *norm2;
	// WARNING: not multiplied by SPIN_DEG because I want emission (singlet ground-state assumed here)
	// I don't think it must be divided by the numeber of electrons, I always used single particle operators [not sure]
	// f-sum rule not valid, beause of the discretised Hamiltonian
     
      }
      else{

         oscillator_frequency_single_p->rvector[ i ] = ZERO;

         oscillator_strength_single_p->rvector[ i ]  = ZERO;

      }


    } /* end i loop */



//#ifdef __DEBUG__

//    if( COMPUTE_ADIABATIC_TRANSITIONS_SINGLE( constants, *state_p ) ) info=1;

//#endif /* __DEBUG__ */

  }


  return info;

}

//------------------------------------------

/* compute_adiabatic_transitions_single */
//BUGFIX: this was originally DEF
int compute_adiabatic_transitions_single( const constants constants, state_p state_p, config_p config_p ){

  /* constants*/
  int        N_levels_single;
  int        sdim;
  /* state */
  rvector_p  adiabatic_PES_single_p;
  matrix_p   adiabatic_states_single_p;
  matrix_p   dummy_matrix_single1_p;
  matrix_p   dummy_matrix_single2_p;
  matrix_p   dummy_matrix_single3_p;
  vector_p   dummy_vector_single1_p;
  vector_p   dummy_vector_single2_p;
  vector_p   dummy_vector_single3_p;
  /* dummies */
  complex    dummy;
  double     norm2;
  double     energy;
  double     strength;
  double     summ;
  int        i, j, k;
  int        comp;
  int        info=0;


  N_levels_single                = constants.N_levels_single;
  sdim                           =  constants.spacial_dimension; 

  adiabatic_PES_single_p         = &(state_p->adiabatic_PES_single);
  adiabatic_states_single_p      = &(state_p->adiabatic_states_single);
  dummy_matrix_single1_p         = &(state_p->dummy_matrix_single1);  
  dummy_matrix_single2_p         = &(state_p->dummy_matrix_single2);  
  dummy_matrix_single3_p         = &(state_p->dummy_matrix_single3);  
  dummy_vector_single1_p         = &(state_p->dummy_vector_single1); 
  dummy_vector_single2_p         = &(state_p->dummy_vector_single2); 
  dummy_vector_single3_p         = &(state_p->dummy_vector_single3); 
  
	    
  if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, *dummy_matrix_single2_p ) ) info=1;
  
  fprintf( stdout, "H_matrix_single\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
  fprintf( stdout, "\n" );
	
			
  /* compute centre_of_mass */
  if( COMPUTE_CENTRE_OF_MASS( constants, config_p->atoms ) ) info=1;
		     
  /* classical_dipole_updating */
  comp=0; //BUGFIX: here there should be the loop over cartesian components...

  if( CLASSICAL_DIPOLE_SINGLE_UPDATE( constants, *state_p, *config_p, comp, *dummy_matrix_single1_p ) ) info=1;

  fprintf( stdout, "Dipole operator (single)\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
  fprintf( stdout, "\n" );


  // sum rule condition
  if( MATRIX_COMMUTATOR( *dummy_matrix_single1_p, *dummy_matrix_single2_p, *dummy_matrix_single3_p ) ) info=1;
  if( MATRIX_COMMUTATOR( *dummy_matrix_single1_p, *dummy_matrix_single3_p, *dummy_matrix_single2_p ) ) info=1;

  fprintf( stdout, "Sum rule identity\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single2_p ) ) info=1;
  fprintf( stdout, "\n" );


  for( i=0; i<N_levels_single; i++ ){
	  
    // copy 
    for( k=0; k<N_levels_single; k++ ){

      dummy_vector_single1_p->vector[ k ] = adiabatic_states_single_p->matrix[ ELECTRON_SINGLE_INDEX( k, i ) ]; 

    } /* end k loop */        
    
    
    /*
      fprintf( stdout, "hole_orbital %d\n", i );
      if( VECTOR_PRINT_PLUS( stdout, *dummy_vector_single1_p ) ) info=1;
      fprintf( stdout, "\n" );
    */


    // init to zero
    summ=0;

    for( j=0; j<N_levels_single; j++ ){

      // copy
      for( k=0; k<N_levels_single; k++ ){

        dummy_vector_single2_p->vector[ k ] = adiabatic_states_single_p->matrix[ ELECTRON_SINGLE_INDEX( k, j ) ];

      } /* end k loop */


      /*
        fprintf( stdout, "particle_orbital %d\n", i );
        if( VECTOR_PRINT_PLUS( stdout, *dummy_vector_single2_p ) ) info=1;
        fprintf( stdout, "\n" );
      */


      // compute oscillator frequency
      energy = adiabatic_PES_single_p->rvector[ i ] - adiabatic_PES_single_p->rvector[ j ];
    
    
      // compute weighted oscillator strength
      if( MATRIX_VECTOR_PRODUCT( *dummy_matrix_single1_p, *dummy_vector_single1_p, *dummy_vector_single3_p ) ) info=1;

      // matrix element
      dummy = SCALAR_PRODUCT( *dummy_vector_single2_p, *dummy_vector_single3_p );


      /*
        fprintf( stdout, "dummy [1]\n" );
        if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
       fprintf( stdout, "\n" );
      */


      // expectation value
      norm2 = CMPLX_NORM( dummy );

      norm2 = norm2 *norm2;

      fprintf( stdout, "# --> dipole[%d,%d] = %le\n", i, j, norm2 );

      
      /*
        fprintf( stdout, "dummy [2]\n" );
        if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
        fprintf( stdout, "\n" );
      */


      // waighted oscillator strength
      strength = TWO /sdim *E_MASS_OVER_AMU *AMU_TO_INTERNAL /( HBAR *HBAR ) *energy *norm2;
      // WARNING: not multiplied by SPIN_DEG because I want emission (singlet ground-state assumed here)
      // I don't think it must be divided by the numeber of electrons, I always used single particle operators [not sure]
      // f-sum rule not valid, beause of the discretised Hamiltonian
	 
      summ += strength;
	 
    } /* end j loop */

    fprintf( stdout, "# Verify sum rule[%d] = %le\n", i, summ );

  } /* end i loop */

 fprintf( stdout, "\n" );

  
  return info;

}

