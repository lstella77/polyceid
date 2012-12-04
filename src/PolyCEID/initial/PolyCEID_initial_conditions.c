
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
#include "PolyCEID_initial_conditions.h"


/* private variables */
matrix          effective_hamiltonian;
matrix          eigenvectors;
rvector         eigenvalues;


/*********************
  FUNCTIONS & MACROS
*********************/

/* initial_condition_rho */

int initial_condition_rho( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             max_rho_index;
  int             sqrt_max_rho_index;
  /* state */
  matrix_p        initial_rho_electron_p;
  matrix_array_p  rho_p;
  matrix_array_p  initial_state_recorded_p;
  matrix_array_p  excited_state_recorded_p;
  matrix_p        initial_condition_atoms_p;
  /* dummies */
  int             i; 
  int             j;
  int             index;
  complex         dummy;
  int             info=0;


  max_rho_index              =  constants.max_rho_index;
  sqrt_max_rho_index         =  constants.sqrt_max_rho_index;

  initial_rho_electron_p     = &(state_p->initial_rho_electron);
  rho_p                      = &(config_p->electrons.rho); 
  initial_state_recorded_p   = &(state_p->initial_state_recorded);
  excited_state_recorded_p   = &(state_p->excited_state_recorded);
  initial_condition_atoms_p  = &(state_p->phonons.initial_condition_atoms);


#ifdef __DEBUG__

  fprintf( stdout, "DOING: initial_condition_rho\n" );

#endif /* __DEBUG__ */


  /* allocations */
  if( ALLOCATE_EFFECTIVE_HAMILTONIAN_UPDATE_WORKSPACE( constants ) ) info=1;

  if( RVECTOR_ALLOCATE( sqrt_max_rho_index, eigenvalues ) ) info=1;

  if( MATRIX_ALLOCATE( sqrt_max_rho_index, sqrt_max_rho_index, eigenvectors ) ) info=1;

  if( MATRIX_ALLOCATE( sqrt_max_rho_index, sqrt_max_rho_index, effective_hamiltonian ) ) info=1;


  //BUGFIX: why here?
  /* initial copy */
  if( RVECTOR_COPY( state_p->positions_saved, config_p->atoms.positions ) ) info=1;


  /* initial guess */
  if( MATRIX_COPY( rho_p->array[ 0 ], *initial_rho_electron_p ) ) info=1;
  
  for( i=1; i<max_rho_index; i++ ){ // WARNING: it starts form 1!
    
    if( MATRIX_ZERO( rho_p->array[ i ] ) ) info=1;

  }


  if( !info ){
	

    /*--------------------*/
    /* Ground state first */
    /*--------------------*/

    /* compute initial_condition_atoms */
    if( COMPUTE_INITIAL_CONDITION_ATOMS( constants, *state_p, *config_p ) ) info=1;


    /* compute the new rho */	
    for( i=0; i<sqrt_max_rho_index; i++ ){
      
      for( j=0; j<sqrt_max_rho_index; j++ ){
	
	//	index = i +sqrt_max_rho_index *j;
	
	index = RHO_INDEX_AUX( i, j );

	/* 
	   WARNING: Here I make use of the excited many-body electronic state, 
	   while the atomic coefficients have been calculated by means of the initial one 
	*/

        dummy = initial_condition_atoms_p->matrix[ index ];

	if( MATRIX_SCALAR_PRODUCT( *initial_rho_electron_p, dummy, rho_p->array[ index ] ) ) info=1;
	
	if( info ) break;
	
      } /* end j loop */
      
    } /* end i loop */



#ifdef __DEBUG_PLUS__
  
    fprintf( stdout, "rho [ground]\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *rho_p ) ) info=1;
    fprintf( stdout, "\n" );
    
#endif /* __DEBUG_PLUS__ */

    //BUGFIX: DEF or TMP?
    /* initial state recorded */
    if( MATRIX_ARRAY_COPY( *initial_state_recorded_p, *rho_p ) ) info=1;

#ifdef __DEBUG_PLUS__
  
    fprintf( stdout, "initial_state_recorded\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *initial_state_recorded_p ) ) info=1;
    fprintf( stdout, "\n" );
    
#endif /* __DEBUG_PLUS__ */


    //BUGFIX: why here?
    if( MU00_UPDATE( constants, *state_p, *config_p ) ) info=1;
  
    /*
      fprintf( stdout, "mu00 [ground]\n" );
      if( MATRIX_PRINT_PLUS( stdout, state_p->electrons.mu00 ) ) info=1;
      fprintf( stdout, "\n" );
    */

    //BUGFIX: why here?
    //if( MATRIX_COPY( config_def_p->electrons.Ehrenfest_frame, state_tmp_p->electrons.Ehrenfest_frame ) ) info=1;

    /*
      fprintf( stdout, "Ehrenfest_frame [ground]\n" );
      if( MATRIX_PRINT_PLUS( stdout, state_p->electrons.Ehrenfest_frame ) ) info=1;
      fprintf( stdout, "\n" );
    */


    /*-----------------------------------------*/
    /* WARNING: from here is the excited state */
    /*-----------------------------------------*/

    /* compute_excited_rho_electron*/
    if( COMPUTE_EXCITED_RHO_ELECTRON( constants, *state_p ) ) info=1;  //WARNING: That's the excited state!
    
    
    /* compute the new rho */	
    for( i=0; i<sqrt_max_rho_index; i++ ){
      
      for( j=0; j<sqrt_max_rho_index; j++ ){
	
	//	index = i +sqrt_max_rho_index *j;
	
	index = RHO_INDEX_AUX( i, j );

	dummy = initial_condition_atoms_p->matrix[ index ];

	/* 
	   WARNING: Here I make use of the excited many-body electronic state, 
	   while the atomic coefficients have been calculated by means of the initial one 
	*/
	
	if( MATRIX_SCALAR_PRODUCT( *initial_rho_electron_p, dummy, rho_p->array[ index ] ) ) info=1;  //WARNING: That's the excited state! 
	
	if( info ) break;
	
      } /* end j loop */
      
    } /* end i loop */

  }


  //BUGFIX: why here?
  /* final copy */
  if( RVECTOR_COPY( config_p->atoms.positions, state_p->positions_saved ) ) info=1;
  
  
#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "rho [excited]\n" );
  if( MATRIX_ARRAY_PRINT_PLUS( stdout, *rho_p ) ) info=1;
  fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */
    
  /* excited state recorded */
  if( MATRIX_ARRAY_COPY( *excited_state_recorded_p, *rho_p ) ) info=1;


#ifdef __DEBUG_PLUS__
  
    fprintf( stdout, "excited_state_recorded\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *excited_state_recorded_p ) ) info=1;
    fprintf( stdout, "\n" );
    
#endif /* __DEBUG_PLUS__ */


  /* deallocations */
  if( MATRIX_FREE( effective_hamiltonian ) ) info=1;

  if( MATRIX_FREE( eigenvectors ) ) info=1;

  if( RVECTOR_FREE( eigenvalues ) ) info=1;

  if( FREE_EFFECTIVE_HAMILTONIAN_UPDATE_WORKSPACE() ) info=1;


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

/*  compute_initial_condition_atoms */

int compute_initial_condition_atoms( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             sqrt_max_rho_index;
  int             initial_ionic_state;
  double          temperature;
  /* state */
  matrix_p        initial_condition_atoms_p;
  /* dummies */
  rvector         bfactor;
  double          egs;
  double          pfunction;
  complex         dummy;
  int             i, j, k;
  int             index;
  int             info=0;


  sqrt_max_rho_index         =  constants.sqrt_max_rho_index;
  initial_ionic_state        =  constants.initial_ionic_state;
  temperature                =  constants.temperature;

  initial_condition_atoms_p  = &(state_p->phonons.initial_condition_atoms);

  // allocate
  if( RVECTOR_ALLOCATE( sqrt_max_rho_index, bfactor ) ) info=1;

  // set to zero
  MATRIX_ZERO( *initial_condition_atoms_p );


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_initial_condition_atoms\n" );

#endif /* __DEBUG_PLUS__ */


  /* create the effective hamiltonian for atoms only */
  if( EFFECTIVE_HAMILTONIAN_UPDATE( constants, *state_p, *config_p, effective_hamiltonian ) ) info=1;
  

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "effective_hamiltonian\n" );
  if( MATRIX_PRINT_PLUS( stdout, effective_hamiltonian ) ) info=1; 
  fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */


  /* solve the variational equation */
  if( DIAGONALISATION( effective_hamiltonian, eigenvectors, eigenvalues ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "eigenvectors\n" );
  if( MATRIX_PRINT_PLUS( stdout, eigenvectors ) ) info=1; 
  fprintf( stdout, "\n");

  fprintf( stdout, "eigenvalues\n" );
  if( RVECTOR_PRINT_PLUS( stdout, eigenvalues ) ) info=1; 
  fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */

  // to seta scale of the energies and avoid overflows
  egs = eigenvalues.rvector[ 0 ];

  pfunction=0.0e0;

  for( k=0; k<sqrt_max_rho_index; k++ ){

     bfactor.rvector[ k ] = exp( -( eigenvalues.rvector[ k ] -egs )/( BOLTZMANN_K *temperature ) );

     pfunction += bfactor.rvector[ k ];

  } /* k loop*/

  // normalisation
  for( k=0; k<sqrt_max_rho_index; k++ ){

    bfactor.rvector[ k ] /= pfunction;

  } /* k loop*/

#ifdef __DEBUG__

  fprintf( stdout, "Boltzmann's factors\n" );
  if( RVECTOR_PRINT_PLUS( stdout, bfactor ) ) info=1; 
  fprintf( stdout, "\n");

#endif /* __DEBUG__ */

  /* copy the coefficients */
  for( i=0; i<sqrt_max_rho_index; i++ ){
      
    for( j=0; j<sqrt_max_rho_index; j++ ){
	
      index = RHO_INDEX_AUX( i, j );

      if( temperature < EPS ){

        dummy = CMPLX_PRODUCT( ( eigenvectors.matrix[ RHO_INDEX_AUX( i, initial_ionic_state ) ] ), CONJ( eigenvectors.matrix[ RHO_INDEX_AUX( j, initial_ionic_state ) ] ) );

        initial_condition_atoms_p->matrix[ index ] = CMPLX_SUM( initial_condition_atoms_p->matrix[ index ], dummy ); // WARNING: overwriting

      }
      else{

        for( k=0; k<sqrt_max_rho_index; k++ ){

	  if( bfactor.rvector[ k ] > 0.0e0 ){
	
	    dummy = CMPLX_PRODUCT( ( eigenvectors.matrix[ RHO_INDEX_AUX( i, k ) ] ), CONJ( eigenvectors.matrix[ RHO_INDEX_AUX( j, k ) ] ) );

            dummy =  CMPLX_PRODUCT( dummy, CMPLX( bfactor.rvector[ k ] ) ); // WARNING: overwriting 

            initial_condition_atoms_p->matrix[ index ] = CMPLX_SUM( initial_condition_atoms_p->matrix[ index ], dummy ); // WARNING: overwriting

	  }  

        } /* k loop*/

      } /* conditional */

    } /* j loop */

  } /* i loop */


  // allocate
  if( RVECTOR_FREE( bfactor ) ) info=1;

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "initial_condition_atoms\n" );
  if( VECTOR_PRINT_PLUS( stdout, *initial_condition_atoms_p ) ) info=1; 
  fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------
