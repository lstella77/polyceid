
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
// Sia benedetto S.Leandro!
#include "PolyCEID_energy_corrections.h"


double          kinetic_energy_corrections=0.0e0;
double          kinetic_energy_corrections_old=0.0e0;

double          potential_energy_corrections=0.0e0;
double          potential_energy_corrections_old=0.0e0;


/*********************
  FUNCTIONS & MACROS
*********************/

int compute_kinetic_energy_corrections( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  int             N_chain;
  int             N_coor_red;
  int             rho_index_border_length;
  int             rho_index_next_to_border_length;
  ivector         relevant_modes;
  rvector         masses_aux;
  double          dt;
  /* state */
  rvector_p       momenta_p;
  matrix_array_p  rho_p;
  matrix_array_p  F_matrix_p;
  matrix_array_p  delta_F_matrix_p;
  matrix_array_p  K_matrix_p;
  matrix_p        dummy_matrix1_p;
  rvector_p       ar_p;
  rvector_p       ap_p;
  double*         kinetic_energy_system_correction_p;
  /* dummies */
  double          kinetic_energy_corrections_partial=0.0e0;
  complex         dummy_n, dummy_m;
  complex         dummy_tot;
  complex         matrix_contribution;
  int             i_rho;
  int             i_mode, j_mode;
  int             i_coor, j_coor;
  int             h_mode;
  int             h_coor;
  int             K_matrix_index1;
  int             n_border, m_border;
  int             n_border_extra, m_border_extra;
  int             n_border_index, m_border_index;
  int             n_border_changed, m_border_changed;
  int             n_next_to_border, m_next_to_border;
  int             n_next_to_border_extra, m_next_to_border_extra;
  int             n_next_to_border_index, m_next_to_border_index;
  int             n_next_to_border_changed, m_next_to_border_changed;
  int             info=0;


  flag_normal_mode_expansion         =  constants.flag_normal_mode_expansion;
  N_chain                            =  constants.N_chain;
  N_coor_red                         =  constants.N_coor_red;
  rho_index_border_length            =  constants.rho_index_border_length;
  rho_index_next_to_border_length    =  constants.rho_index_next_to_border_length;
  relevant_modes                     =  constants.relevant_modes;
  masses_aux                         =  config_p->atoms.masses_aux;
  dt                                 =  constants.dt;

  momenta_p                          = &(config_p->atoms.momenta);
  rho_p                              = &(config_p->electrons.rho);
  F_matrix_p                         = &(config_p->electrons.F_matrix);
  delta_F_matrix_p                   = &(config_p->electrons.delta_F_matrix);
  K_matrix_p                         = &(config_p->electrons.K_matrix);
  dummy_matrix1_p                    = &(state_p->dummy_matrix1);
  ar_p                               = &(state_p->phonons.ar);
  ap_p                               = &(state_p->phonons.ap);
  kinetic_energy_system_correction_p = &(config_p->kinetic_energy_system_correction);


  if( N_chain ){

#ifdef __TIME_SCALE__

    dt *= exp( config_p->thermostat.positions.rvector[ 0 ] );

#endif /* __TIME_SCALE__ */

  }
		  

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_kinetic_energy_corrections\n" );
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  // copy
  kinetic_energy_corrections_old = kinetic_energy_corrections;


  if( flag_normal_mode_expansion ){

    if( RVECTOR_COPY( state_p->momenta_saved, *momenta_p ) ) info=1;

    if( MATRIX_ARRAY_COPY( state_p->F_matrix_saved, *F_matrix_p ) ) info=1;

    if( MATRIX_ARRAY_COPY( state_p->delta_F_matrix_saved, *delta_F_matrix_p ) ) info=1;

    if( MATRIX_ARRAY_COPY( state_p->K_matrix_saved, *K_matrix_p ) ) info=1;

    /* transform hamiltonian */
    if( TRANSFORM_HAMILTONIAN( constants, *state_p, *config_p, config_p->electrons.delta_F_matrix, config_p->electrons.K_matrix ) ) info=1;

    if( TRANSFORM_HAMILTONIAN_AUX( constants, *state_p, *config_p, *F_matrix_p, *momenta_p ) ) info=1;

  }


  // set to zero
  kinetic_energy_corrections = 0.0e0;


  if( !info ){

    /*
      1ST TERM: 1B
    */

    /*
      
    fprintf( stdout, "--> 1ST TERM: [1B]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){

	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion

	      matrix_contribution = MATRIX_TRACE( rho_p->array[ i_rho ] );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ i_coor ] *( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			       );

	      /* dummy_m */
	      dummy_m = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor] ) *( ap_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			       );

	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;

	      
	      // partial updating
	      kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!
	      

	    } /* end conditional */
	    
	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [1B] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      2ND TERM: 2A
    */

    /*
    
    fprintf( stdout, "--> 2ND TERM: [2A]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_next_to_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_next_to_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_next_to_border_changed = %d, m_border_changed = %d\n", n_next_to_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_next_to_border ); // WARNING: note the index inversion

	      matrix_contribution = MATRIX_TRACE( rho_p->array[ i_rho ] );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ i_coor ] *( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			      );

	      /* dummy_m */
	      dummy_m = CMPLX( 
			      OSQRT2 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor ] ) *momenta_p->rvector[ j_coor ]
			      *sqrt( dummy_rho_index2.ivector[ j_mode ] +1.0e0 )
			      );

	      dummy_m = CMPLX_TIMES_I( dummy_m ); //WARNING: overwriting


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting


	      // breaking instruction
	      if( info ) break;

	    
	      // partial updating
	      kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [2A] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      3RD TERM: 2B
    */

    /*
    
    fprintf( stdout, "--> 3RD TERM: [2B]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){

	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	    
	    m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_next_to_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_next_to_border_changed == n_next_to_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d\n", n_next_to_border, i_coor, m_next_to_border, j_coor );
	      fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );
	      
	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion

	      matrix_contribution = MATRIX_TRACE( rho_p->array[ i_rho ] );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ i_coor ] *( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			      );

	      /* dummy_m */
	      dummy_m = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor] ) *( ap_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			      );

	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	    
	      // breaking instruction
	      if( info ) break;

	      
	      // partial updating
	      kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!
	      
	      
	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [2B] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      4TH TERM: 3A
    */

    /*
    
    fprintf( stdout, "--> 4TH TERM: [3A]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( n_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion

	      if( MATRIX_MATRIX_PRODUCT( delta_F_matrix_p->array[ i_coor ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -OSQRT2 *( ar_p->rvector[ i_coor ] )
			      *sqrt( dummy_rho_index1.ivector[ i_mode ] +1.0e0 )
			       );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      OSQRT2 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor ] ) *momenta_p->rvector[ j_coor ]
			      *sqrt( dummy_rho_index2.ivector[ j_mode ] +1.0e0 )
			      );

	      dummy_m = CMPLX_TIMES_I( dummy_m ); //WARNING: overwriting


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	    
	      // breaking instruction
	      if( info ) break;

	    
	      // partial updating
	      kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [3A] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      5TH TERM: 3B
    */

    /*
    
    fprintf( stdout, "--> 5TH TERM: [3B]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( n_border_extra, i_mode );
	    
	    m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_next_to_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_next_to_border_changed == n_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d\n", n_border, i_coor, m_next_to_border, j_coor );
	      fprintf( stdout, "---> n_border_changed = %d, m_next_to_border_changed = %d\n", n_border_changed, m_next_to_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_next_to_border, n_border ); // WARNING: note the index inversion

	      if( MATRIX_MATRIX_PRODUCT( delta_F_matrix_p->array[ i_coor ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -OSQRT2 *( ar_p->rvector[ i_coor ] )
			      *sqrt( dummy_rho_index1.ivector[ i_mode ] +1.0e0 )
			       );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor] ) *( ap_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			      );

	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting
	  
	    
	      // breaking instruction
	      if( info ) break;

	      
	      // partial updating
	      kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [3B] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      6TH TERM: 4B
    */

    /*
    
    fprintf( stdout, "--> 6TH TERM: [4B]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion

	      K_matrix_index1 = COORDINATE_INDEX( i_coor, i_coor );

	      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			       );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor] ) *( ap_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			      );

	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;
	    
	    
	      // partial updating
	      kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	    
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [4B] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      7TH TERM: 5A
    */

   /*
    
    fprintf( stdout, "--> 7TH TERM: [5A]\n" );
    fprintf( stdout, "\n" );
    
   */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_next_to_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_next_to_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_next_to_border_changed = %d, m_border_changed = %d\n", n_next_to_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */
	      

	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_next_to_border ); // WARNING: note the index inversion

	      K_matrix_index1 = COORDINATE_INDEX( i_coor, i_coor );

	      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			       );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      OSQRT2 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor ] ) *momenta_p->rvector[ j_coor ]
			      *sqrt( dummy_rho_index2.ivector[ j_mode ] +1.0e0 )
			      );

	      dummy_m = CMPLX_TIMES_I( dummy_m ); //WARNING: overwriting


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;

	    
	      // partial updating
	      kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [5A] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      8TH TERM: 5B
    */

    /*
    
    fprintf( stdout, "--> 8TH TERM: [5B]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];

	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	    
	    m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_next_to_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_next_to_border_changed == n_next_to_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d\n", n_next_to_border, i_coor, m_next_to_border, j_coor );
	      fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion

	      K_matrix_index1 = COORDINATE_INDEX( i_coor, i_coor );

	      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			       );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor] ) *( ap_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			      );

	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;

	    
	      // partial updating
	      kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [5B] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      9TH TERM: 6B
    */

   /*
    
    fprintf( stdout, "--> 9TH TERM: [6B]\n" );
    fprintf( stdout, "\n" );
    
   */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( h_mode=0; h_mode<N_coor_red; h_mode++ ){
	  
	    h_coor = relevant_modes.ivector[ h_mode ]; // i_mode to cartesian coordinate
	  
	    for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	      j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	      //
	      n_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( n_border_extra, i_mode, h_mode );
	    
	      m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_border_extra, j_mode );
	
	    
	      /* conditional */
	      if( m_border_changed == n_border_changed && h_coor != i_coor){
	    
		/*

		fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d, h_coor = %d\n", n_border, i_coor, m_border, j_coor, h_coor );
		fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );

		*/


		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion
		
		K_matrix_index1 = COORDINATE_INDEX( i_coor, h_coor );

		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;
		
		matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );
		

		/* dummy_n */
		dummy_n = CMPLX( 
				0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ h_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ h_mode ] +1.0e0 ) )
				 );
		

		/* dummy_m */
		dummy_m = CMPLX( 
				-0.25e0 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor] ) *( ap_p->rvector[ j_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
				);

		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		
		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
		// breaking instruction
		if( info ) break;
	      
	    
		// partial updating
		kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!

	      
	      } /* end conditional */

	      
	    } /* j_mode loop */

	    // breaking instruction
	    if( info ) break;
	    
	  } /* h_mode loop */

	  // breaking instruction
	  if( info ) break;
	    
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [6B] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      10TH TERM: 7A
    */

   /*
    
    fprintf( stdout, "--> 10TH TERM: [7A]\n" );
    fprintf( stdout, "\n" );
    
   */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( h_mode=0; h_mode<N_coor_red; h_mode++ ){
	  
	    h_coor = relevant_modes.ivector[ h_mode ]; // i_mode to cartesian coordinate
	  

	    for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	      j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	      
	      //
	      n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( n_next_to_border_extra, i_mode, h_mode );
	      
	      m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( m_border_extra, j_mode );
	      
	    
	      /* conditional */
	      if( m_border_changed == n_next_to_border_changed && h_coor != i_coor ){
		
		/*
		
		fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_border = %d, j_coor = %d, h_coor = %d\n", n_next_to_border, i_coor, m_border, j_coor, h_coor );
		fprintf( stdout, "---> n_next_to_border_changed = %d, m_border_changed = %d\n", n_next_to_border_changed, m_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );

		*/


		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_border, n_next_to_border ); // WARNING: note the index inversion
		
		K_matrix_index1 = COORDINATE_INDEX( i_coor, h_coor );
		
		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

		matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );
		
		
		/* dummy_n */
		dummy_n = CMPLX( 
				0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ h_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ h_mode ] +1.0e0 ) )
				 );

		
		/* dummy_m */
		dummy_m = CMPLX( 
				OSQRT2 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor ] ) *momenta_p->rvector[ j_coor ]
				*sqrt( dummy_rho_index2.ivector[ j_mode ] +1.0e0 )
				);

		dummy_m = CMPLX_TIMES_I( dummy_m ); //WARNING: overwriting


		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
		// breaking instruction
		if( info ) break;
		
		
		// partial updating
		kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */

	  
	    } /* j_mode loop */

	    // breaking instruction
	    if( info ) break;
	  
	  } /* h_mode loop */

	  // breaking instruction
	  if( info ) break;

	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [7A] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      11TH TERM: 7B
    */

    /*
    
    fprintf( stdout, "--> 11TH TERM: [7B]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    kinetic_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( h_mode=0; h_mode<N_coor_red; h_mode++ ){
	  
	    h_coor = relevant_modes.ivector[ h_mode ]; // i_mode to cartesian coordinate
	    	  
	    for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	      
	      j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	      
	      //
	      n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( n_next_to_border_extra, i_mode, h_mode );
	      
	      m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_next_to_border_extra, j_mode );
	
	    
	      /* conditional */
	      if( m_next_to_border_changed == n_next_to_border_changed && h_coor != i_coor ){
	    
		/*

		fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d, h_coor = %d\n", 
			 n_next_to_border, i_coor, m_next_to_border, j_coor, h_coor );
		fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );

		*/


		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion

		K_matrix_index1 = COORDINATE_INDEX( i_coor, h_coor );
		
		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

		matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


		/* dummy_n */
		dummy_n = CMPLX( 
				0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ h_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ h_mode ] +1.0e0 ) )
				 );


		/* dummy_m */
		dummy_m = CMPLX( 
				-0.25e0 /masses_aux.rvector[ j_coor ] *( ap_p->rvector[ j_coor] ) *( ap_p->rvector[ j_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
				);

		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		
		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
		// breaking instruction
		if( info ) break;
		
	    
		// partial updating
		kinetic_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */

	  
	    } /* j_mode loop */
	    
	    // breaking instruction
	    if( info ) break;
	    
	  } /* j_mode loop */
	    
	    // breaking instruction
	    if( info ) break;

	} /* i_mode loop */
	  
	// breaking instruction
	if( info ) break;
	  
      } /* m_next_to_border_index loop */
      
      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [7B] partial correction is %le\n", kinetic_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    kinetic_energy_corrections += kinetic_energy_corrections_partial;

  }

  //----------

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "--------------\n" );
  fprintf( stdout, "total kinetic correction is %le\n", kinetic_energy_corrections );
  fprintf( stdout, "--------------\n" );
  
#endif /* __DEBUG_PLUS__ */


  // update
  *kinetic_energy_system_correction_p += 0.5e0 *( kinetic_energy_corrections_old +kinetic_energy_corrections ) *dt;


  /*

  fprintf( stdout, "# kinetic_energy_corrections_old is %le\n", kinetic_energy_corrections_old );
  fprintf( stdout, "# kinetic_energy_corrections     is %le\n", kinetic_energy_corrections     );
    
  */


  /* copy back */
  if( flag_normal_mode_expansion ){

    if( RVECTOR_COPY( *momenta_p, state_p->momenta_saved ) ) info=1;

    if( MATRIX_ARRAY_COPY( *F_matrix_p, state_p->F_matrix_saved ) ) info=1;

    if( MATRIX_ARRAY_COPY( *delta_F_matrix_p, state_p->delta_F_matrix_saved ) ) info=1;

    if( MATRIX_ARRAY_COPY( *K_matrix_p, state_p->K_matrix_saved ) ) info=1;

  }


  return info;

}

//------------------------------------------

int compute_potential_energy_corrections( const constants constants, state_p state_p, config_p config_p ){


  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  unsigned short int  flag_no_Ehrenfest_frame;
  int             N_chain;
  int             N_coor;
  int             N_coor_red;
  int             rho_index_border_length;
  int             rho_index_next_to_border_length;
  ivector         relevant_modes;
  rvector         masses_aux;
  double          dt;
  /* state */
  rvector_p       momenta_p;
  matrix_array_p  rho_p;
  matrix_array_p  mu10_p;
  matrix_array_p  mu20_p;
  matrix_p        H_matrix_p;
  matrix_array_p  F_matrix_p;
  matrix_array_p  delta_F_matrix_p;
  matrix_array_p  K_matrix_p;
  matrix_p        dummy_matrix1_p;
  matrix_p        dummy_matrix2_p;
  rvector_p       ar_p;
  rvector_p       ap_p;
  double*         potential_energy_system_correction_p;
  /* dummies */
  double          potential_energy_corrections_partial=0.0e0;
  complex         dummy_n, dummy_m;
  complex         dummy_tot;
  complex         matrix_contribution;
  int             i_rho;
  int             i_mode, j_mode;
  int             i_coor, j_coor;
  int             h_mode, k_mode;
  int             h_coor, k_coor;
  int             K_matrix_index1, K_matrix_index2;
  int             n_border, m_border;
  int             n_border_extra, m_border_extra;
  int             n_border_index, m_border_index;
  int             n_border_changed, m_border_changed;
  int             n_next_to_border, m_next_to_border;
  int             n_next_to_border_extra, m_next_to_border_extra;
  int             n_next_to_border_index, m_next_to_border_index;
  int             n_next_to_border_changed, m_next_to_border_changed;
  int             info=0;


  flag_normal_mode_expansion            =  constants.flag_normal_mode_expansion;
  flag_no_Ehrenfest_frame               =  constants.flag_no_Ehrenfest_frame;
  N_chain                               =  constants.N_chain;
  N_coor                                =  constants.N_coor;
  N_coor_red                            =  constants.N_coor_red;
  rho_index_border_length               =  constants.rho_index_border_length;
  rho_index_next_to_border_length       =  constants.rho_index_next_to_border_length;
  relevant_modes                        =  constants.relevant_modes;
  masses_aux                            =  config_p->atoms.masses_aux;
  dt                                    =  constants.dt;

  momenta_p                             = &(config_p->atoms.momenta);
  rho_p                                 = &(config_p->electrons.rho);
  mu10_p                                = &(config_p->electrons.mu10);
  mu20_p                                = &(config_p->electrons.mu20);
  H_matrix_p                            = &(config_p->electrons.H_matrix);
  F_matrix_p                            = &(config_p->electrons.F_matrix);
  delta_F_matrix_p                      = &(config_p->electrons.delta_F_matrix);
  K_matrix_p                            = &(config_p->electrons.K_matrix);
  dummy_matrix1_p                       = &(state_p->dummy_matrix1);
  dummy_matrix2_p                       = &(state_p->dummy_matrix2);
  ar_p                                  = &(state_p->phonons.ar);
  ap_p                                  = &(state_p->phonons.ap);
  potential_energy_system_correction_p  = &(config_p->potential_energy_system_correction);


  if( N_chain ){

#ifdef __TIME_SCALE__

    dt *= exp( config_p->thermostat.positions.rvector[ 0 ] );

#endif /* __TIME_SCALE__ */    

  }
		  

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_potential_energy_corrections\n" );
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  // copy
  potential_energy_corrections_old = potential_energy_corrections;


  if( flag_normal_mode_expansion ){

    if( RVECTOR_COPY( state_p->momenta_saved, *momenta_p ) ) info=1;

    if( MATRIX_ARRAY_COPY( state_p->F_matrix_saved, *F_matrix_p ) ) info=1;

    if( MATRIX_ARRAY_COPY( state_p->delta_F_matrix_saved, *delta_F_matrix_p ) ) info=1;

    if( MATRIX_ARRAY_COPY( state_p->K_matrix_saved, *K_matrix_p ) ) info=1;

    /* transform hamiltonian */
    if( TRANSFORM_HAMILTONIAN( constants, *state_p, *config_p, config_p->electrons.delta_F_matrix, config_p->electrons.K_matrix ) ) info=1;

    if( TRANSFORM_HAMILTONIAN_AUX( constants, *state_p, *config_p, *F_matrix_p, *momenta_p ) ) info=1;

  }


  // set to zero
  potential_energy_corrections = 0.0e0;


  if( !info ){

    /*
      1ST TERM: 1D
    */

    /*
    
    fprintf( stdout, "--> 1ST TERM: [1D]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion

	      K_matrix_index2 = COORDINATE_INDEX( j_coor, j_coor );

	      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index2 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ i_coor ] *( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			      );

	      /* dummy_m */
	      dummy_m = CMPLX( 
			      0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			       );


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;

	      
	      // partial updating
	      potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [1D] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      2ND TERM: 1E
    */

    /*
    
    fprintf( stdout, "--> 2ND TERM: [1E]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate

	    for( k_mode=0; k_mode<N_coor_red; k_mode++ ){
	    
	      k_coor = relevant_modes.ivector[ k_mode ]; // i_mode to cartesian coordinate

	    
	      //
	      n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_border_extra, i_mode );
	      
	      m_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( m_border_extra, j_mode, k_mode );
	      
	      
	      /* conditional */
	      if( m_border_changed == n_border_changed && k_coor != j_coor ){
		
		/*
		
		fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d, k_coor = %d\n", n_border, i_coor, m_border, j_coor, k_coor );
		fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );
		
		*/
		
		
		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion
		
		K_matrix_index2 = COORDINATE_INDEX( j_coor, k_coor );
		
		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index2 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;
		
		matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );
		
		
		/* dummy_n */
		dummy_n = CMPLX( 
				-0.25e0 /masses_aux.rvector[ i_coor ] *( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
				);
		
		/* dummy_m */
		dummy_m = CMPLX( 
				0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ k_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ k_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		
		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting
		
		
		// breaking instruction
		if( info ) break;

	    
		// partial updating
		potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */

	  
	    } /* k_mode loop */

	    // breaking instruction
	    if( info ) break;
	    
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;

	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [1E] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      3RD TERM: 2C
    */

    /*
    
    fprintf( stdout, "--> 3RD TERM: [2C]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_next_to_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_next_to_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_next_to_border_changed = %d, m_border_changed = %d\n", n_next_to_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );
	      
	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_next_to_border ); // WARNING: note the index inversion

	      if( MATRIX_MATRIX_PRODUCT( F_matrix_p->array[ j_coor ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ i_coor ] *( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			      );

	      /* dummy_m */
	      dummy_m = CMPLX( 
			      -OSQRT2 *( ar_p->rvector[ j_coor] )
			      *sqrt( dummy_rho_index2.ivector[ j_mode ] +1.0e0 )
			       );


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;
	      
	    
	      // partial updating
	      potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [2C] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      4TH TERM: 2D
    */

    /*
    
    fprintf( stdout, "--> 4TH TERM: [2D]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	    
	    m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_next_to_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_next_to_border_changed == n_next_to_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d\n", n_next_to_border, i_coor, m_next_to_border, j_coor );
	      fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion

	      K_matrix_index2 = COORDINATE_INDEX( j_coor, j_coor );

	      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index2 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -0.25e0 /masses_aux.rvector[ i_coor ] *( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			      );

	      /* dummy_m */
	      dummy_m = CMPLX( 
			      0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			       );


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;
	      
	      
	      // partial updating
	      potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!
	      
	      
	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [2D] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      5TH TERM: 2E
    */

    /*
    
    fprintf( stdout, "--> 5TH TERM: [2E]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate

	    for( k_mode=0; k_mode<N_coor_red; k_mode++ ){
	    
	      k_coor = relevant_modes.ivector[ k_mode ]; // i_mode to cartesian coordinate

	    
	      //
	      n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	      
	      m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( m_next_to_border_extra, j_mode, k_mode );
	      
	      
	      /* conditional */
	      if( m_next_to_border_changed == n_next_to_border_changed && k_coor != j_coor){
		
		/*
		
		fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d, k_coor = %d\n", 
			 n_next_to_border, i_coor, m_next_to_border, j_coor, k_coor );
		fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );
		
		*/
		
		
		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion
		
		K_matrix_index2 = COORDINATE_INDEX( j_coor, k_coor );
		
		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index2 ], rho_p->array[ i_rho ], *dummy_matrix1_p ) ) info=1;
		
		matrix_contribution = MATRIX_TRACE( *dummy_matrix1_p );
		
		
		/* dummy_n */
		dummy_n = CMPLX( 
				-0.25e0 /masses_aux.rvector[ i_coor ] *( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
				);
		
		/* dummy_m */
		dummy_m = CMPLX( 
				0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ k_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ k_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		
		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting
		
		
		// breaking instruction
		if( info ) break;

	    
		// partial updating
		potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */

	  
	    } /* k_mode loop */

	    // breaking instruction
	    if( info ) break;
	    
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;

	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [2E] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      6TH TERM: 3C
    */

   /*
    
    fprintf( stdout, "--> 3TH TERM: [3C]\n" );
    fprintf( stdout, "\n" );
    
   */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( n_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion

	      if( MATRIX_MATRIX_PRODUCT( delta_F_matrix_p->array[ i_coor ], F_matrix_p->array[ j_coor ], *dummy_matrix1_p ) ) info=1;

	      if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -OSQRT2 *( ar_p->rvector[ i_coor ] )
			      *sqrt( dummy_rho_index1.ivector[ i_mode ] +1.0e0 )
			       );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      -OSQRT2 *( ar_p->rvector[ j_coor] )
			      *sqrt( dummy_rho_index2.ivector[ j_mode ] +1.0e0 )
			       );


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;
	      
	      
	      // partial updating
	      potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [3C] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      7TH TERM: 3D
    */

    /*
    
    fprintf( stdout, "--> 4TH TERM: [3D]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( n_border_extra, i_mode );
	    
	    m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_next_to_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_next_to_border_changed == n_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d\n", n_border, i_coor, m_next_to_border, j_coor );
	      fprintf( stdout, "---> n_border_changed = %d, m_next_to_border_changed = %d\n", n_border_changed, m_next_to_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_next_to_border, n_border ); // WARNING: note the index inversion

	      K_matrix_index2 = COORDINATE_INDEX( j_coor, j_coor );

	      if( MATRIX_MATRIX_PRODUCT( delta_F_matrix_p->array[ i_coor ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;

	      if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      -OSQRT2 *( ar_p->rvector[ i_coor ] )
			      *sqrt( dummy_rho_index1.ivector[ i_mode ] +1.0e0 )
			       );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			       );


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;
	      
	    
	      // partial updating
	      potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [3D] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      8TH TERM: 3E
    */

    /*
    
    fprintf( stdout, "--> 8TH TERM: [3E]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate

	    for( k_mode=0; k_mode<N_coor_red; k_mode++ ){
	    
	      k_coor = relevant_modes.ivector[ k_mode ]; // i_mode to cartesian coordinate

	    
	      //
	      n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( n_border_extra, i_mode );
	      
	      m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( m_next_to_border_extra, j_mode, k_mode );

	      
	      /* conditional */
	      if( m_next_to_border_changed == n_border_changed && k_coor != j_coor ){
		
		/*
		
		fprintf( stdout, "---> n_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d, k_coor = %d\n", n_border, i_coor, m_next_to_border, j_coor, k_coor );
		fprintf( stdout, "---> n_border_changed = %d, m_next_to_border_changed = %d\n", n_border_changed, m_next_to_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );
		
		*/
		
		
		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_next_to_border, n_border ); // WARNING: note the index inversion
		
		K_matrix_index2 = COORDINATE_INDEX( j_coor, k_coor );
		
		if( MATRIX_MATRIX_PRODUCT( delta_F_matrix_p->array[ i_coor ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;

		if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
		
		matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );
		
		
		/* dummy_n */
		dummy_n = CMPLX( 
				-OSQRT2 *( ar_p->rvector[ i_coor ] )
				*sqrt( dummy_rho_index1.ivector[ i_mode ] +1.0e0 )
				);
		
		
		/* dummy_m */
		dummy_m = CMPLX( 
				0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ k_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ k_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		
		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting
		
		
		// breaking instruction
		if( info ) break;

	    
		// partial updating
		potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */

	  
	    } /* k_mode loop */

	    // breaking instruction
	    if( info ) break;
	    
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;

	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [3E] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      9TH TERM: 4D
    */

    /*
    
    fprintf( stdout, "--> 9TH TERM: [4D]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion

	      K_matrix_index1 = COORDINATE_INDEX( i_coor, i_coor );

	      K_matrix_index2 = COORDINATE_INDEX( j_coor, j_coor );

	      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;

	      if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			      );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			      );


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;

	    
	      // partial updating
	      potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	    } /* end conditional */

	  
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [4D] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      10TH TERM: 4E
    */

    /*
    
    fprintf( stdout, "--> 10TH TERM: [4E]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate

	    for( k_mode=0; k_mode<N_coor_red; k_mode++ ){
	    
	      k_coor = relevant_modes.ivector[ k_mode ]; // i_mode to cartesian coordinate

	    
	      //
	      n_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_border_extra, i_mode );
	      
	      m_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( m_border_extra, j_mode, k_mode );
	      
	      
	      /* conditional */
	      if( m_border_changed == n_border_changed && k_coor != j_coor ){
		
		/*
		
		fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d, k_coor = %d\n", n_border, i_coor, m_border, j_coor, k_coor );
		fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );
		
		*/
		
		
		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion

		K_matrix_index1 = COORDINATE_INDEX( i_coor, i_coor );

		K_matrix_index2 = COORDINATE_INDEX( j_coor, k_coor );
		
		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;

		if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
		
		matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );
		
		
		/* dummy_n */
		dummy_n = CMPLX( 
				0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ i_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_m */
		dummy_m = CMPLX( 
				0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ k_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ k_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		
		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting
		
		
		// breaking instruction
		if( info ) break;

	    
		// partial updating
		potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */
	  

	    } /* k_mode loop */

	    // breaking instruction
	    if( info ) break;
	    
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;

	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [4E] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      11TH TERM: 5C
    */

    /*
    
    fprintf( stdout, "--> 11TH TERM: [5C]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	    
	    m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( m_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_border_changed == n_next_to_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_border = %d, j_coor = %d\n", n_next_to_border, i_coor, m_border, j_coor );
	      fprintf( stdout, "---> n_next_to_border_changed = %d, m_border_changed = %d\n", n_next_to_border_changed, m_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_border, n_next_to_border ); // WARNING: note the index inversion

	      K_matrix_index1 = COORDINATE_INDEX( i_coor, i_coor );
	      
	      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], F_matrix_p->array[ j_coor ], *dummy_matrix1_p ) ) info=1;
	      
	      if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
	      
	      matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			      );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      -OSQRT2 *( ar_p->rvector[ j_coor] )
			      *sqrt( dummy_rho_index2.ivector[ j_mode ] +1.0e0 )
			       );


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;

	    
	      // partial updating
	      potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!

	  
	    } /* end conditional */
	  

	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [5C] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      12TH TERM: 5D
    */

   /*
    
    fprintf( stdout, "--> 9TH TERM: [5D]\n" );
    fprintf( stdout, "\n" );
    
   */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	    //
	    n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	    
	    m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_next_to_border_extra, j_mode );
	
	    
	    /* conditional */
	    if( m_next_to_border_changed == n_next_to_border_changed ){
	    
	      /*

	      fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d\n", n_next_to_border, i_coor, m_next_to_border, j_coor );
	      fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
	      fprintf( stdout, "\n" );
	      fprintf( stdout, "dummy_rho_indices\n" );
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
	      fprintf( stdout, "\n" );

	      */


	      /* matrix contribution */
	      i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion

	      K_matrix_index1 = COORDINATE_INDEX( i_coor, i_coor );

	      K_matrix_index2 = COORDINATE_INDEX( j_coor, j_coor );

	      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;

	      if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;

	      matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );


	      /* dummy_n */
	      dummy_n = CMPLX( 
			      0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ i_coor ] )
			      *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
			      );


	      /* dummy_m */
	      dummy_m = CMPLX( 
			      0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ j_coor] )
			      *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
			      );


	      /* dummy_tot */
	      dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );

	      dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

	      dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
	      // breaking instruction
	      if( info ) break;

	    
	      // partial updating
	      potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!
	      

	    } /* end conditional */
	  

	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [5D] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      13TH TERM: 5E
    */

    /*
    
    fprintf( stdout, "--> 13TH TERM: [5E]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];

	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	    j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate

	    for( k_mode=0; k_mode<N_coor_red; k_mode++ ){
	    
	      k_coor = relevant_modes.ivector[ k_mode ]; // i_mode to cartesian coordinate

	    
	      //
	      n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( n_next_to_border_extra, i_mode );
	      
	      m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( m_next_to_border_extra, j_mode, k_mode );
	      
	      
	      /* conditional */
	      if( m_next_to_border_changed == n_next_to_border_changed && k_coor != j_coor ){
		
		/*
		
		fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, m_next_to_border = %d, j_coor = %d, k_coor = %d\n", 
		n_next_to_border, i_coor, m_next_to_border, j_coor, k_coor );
		fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );
		
		*/
		
		
		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion

		K_matrix_index1 = COORDINATE_INDEX( i_coor, i_coor );

		K_matrix_index2 = COORDINATE_INDEX( j_coor, k_coor );
		
		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;

		if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
		
		matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );
		
		
		/* dummy_n */
		dummy_n = CMPLX( 
				0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ i_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_m */
		dummy_m = CMPLX( 
				0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ k_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ k_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		
		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting
		
		
		// breaking instruction
		if( info ) break;
		
	    
		// partial updating
		potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */

	  
	    } /* k_mode loop */

	    // breaking instruction
	    if( info ) break;
	    
	  } /* j_mode loop */

	  // breaking instruction
	  if( info ) break;

	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [5E] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      14TH TERM: 6D
    */

   /*
    
    fprintf( stdout, "--> 14TH TERM: [6D]\n" );
    fprintf( stdout, "\n" );
    
   */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( h_mode=0; h_mode<N_coor_red; h_mode++ ){
	  
	    h_coor = relevant_modes.ivector[ h_mode ]; // i_mode to cartesian coordinate
	  
	    for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	      j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	      
	      //
	      n_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( n_border_extra, i_mode, h_mode );
	      
	      m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_border_extra, j_mode );
	      
	      
	      /* conditional */
	      if( m_border_changed == n_border_changed && h_coor != i_coor ){
		
		/*
		
		fprintf( stdout, "---> n_border = %d, i_coor = %d, h_coor = %d, m_border = %d, j_coor = %d\n", n_border, i_coor, h_coor, m_border, j_coor );
		fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );
		
		*/
		

		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion
		
		K_matrix_index1 = COORDINATE_INDEX( i_coor, h_coor );
		
		K_matrix_index2 = COORDINATE_INDEX( j_coor, j_coor );

		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;
		
		if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
		
		matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );
	      
		
		/* dummy_n */
		dummy_n = CMPLX( 
				0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ h_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ h_mode ] +1.0e0 ) )
				);


		/* dummy_m */
		dummy_m = CMPLX( 
				0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ j_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
		// breaking instruction
		if( info ) break;

	    
		// partial updating
		potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */

	      
	    } /* j_mode loop */
	    
	    // breaking instruction
	    if( info ) break;
	  
	  } /* h_mode loop */
	    
	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [6D] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      15TH TERM: 6E
    */

    /*
    
    fprintf( stdout, "--> 15TH TERM: [6E]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_border_index=0; n_border_index<rho_index_border_length; n_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_border = rho_index_border.ivector[ n_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_border_extra = rho_index_aux_translate.ivector[ n_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( h_mode=0; h_mode<N_coor_red; h_mode++ ){
	  
	    h_coor = relevant_modes.ivector[ h_mode ]; // i_mode to cartesian coordinate
	  
	    for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	      
	      j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	      
	      for( k_mode=0; k_mode<N_coor_red; k_mode++ ){
		
		k_coor = relevant_modes.ivector[ k_mode ]; // i_mode to cartesian coordinate
		
	    
		//
		n_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( n_border_extra, i_mode, h_mode );
		
		m_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( m_border_extra, j_mode, k_mode );
		
	      
		/* conditional */
		if( m_border_changed == n_border_changed && h_coor != i_coor && k_coor != j_coor ){
		  
		  /*
		
		  fprintf( stdout, "---> n_border = %d, i_coor = %d, m_border = %d, j_coor = %d, k_coor = %d\n", n_border, i_coor, m_border, j_coor, k_coor );
		  fprintf( stdout, "---> n_border_changed = %d, m_border_changed = %d\n", n_border_changed, m_border_changed );
		  fprintf( stdout, "\n" );
		  fprintf( stdout, "dummy_rho_indices\n" );
		  if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		  if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		  fprintf( stdout, "\n" );
		
		  */
		
		  
		  /* matrix contribution */
		  i_rho = RHO_INDEX_AUX( m_border, n_border ); // WARNING: note the index inversion
		  
		  K_matrix_index1 = COORDINATE_INDEX( i_coor, h_coor );
		  
		  K_matrix_index2 = COORDINATE_INDEX( j_coor, k_coor );
		
		  if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;
		  
		  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
		
		  matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );
		
		
		  /* dummy_n */
		  dummy_n = CMPLX( 
				  0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ h_coor ] )
				  *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ h_mode ] +1.0e0 ) )
				  );
		
		
		  /* dummy_m */
		  dummy_m = CMPLX( 
				  0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ k_coor] )
				  *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ k_mode ] +1.0e0 ) )
				  );
		
		
		  /* dummy_tot */
		  dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		  dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		
		  dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting
		
		
		  // breaking instruction
		  if( info ) break;
		
		
		  // partial updating
		  potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


		} /* end conditional */

	  
	      } /* k_mode loop */

	      // breaking instruction
	    if( info ) break;
	    
	    } /* j_mode loop */

	    // breaking instruction
	    if( info ) break;

	  } /* h_mode loop */

	  // breaking instruction
	  if( info ) break;

	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [6E] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      16TH TERM: 7C
    */

    /*
    
    fprintf( stdout, "--> 16TH TERM: [7C]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_border_index=0; m_border_index<rho_index_border_length; m_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_border = rho_index_border.ivector[ m_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_border_extra = rho_index_aux_translate.ivector[ m_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( h_mode=0; h_mode<N_coor_red; h_mode++ ){
	  
	    h_coor = relevant_modes.ivector[ h_mode ]; // i_mode to cartesian coordinate
	  
	    for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	    
	      j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	    
	      //
	      n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( n_next_to_border_extra, i_mode, h_mode );
	      
	      m_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( m_border_extra, j_mode );
	
	    
	      /* conditional */
	      if( m_border_changed == n_next_to_border_changed && h_coor != i_coor ){
		
		/*
		
		fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, h_coor = %d, m_border = %d, j_coor = %d\n", n_next_to_border, i_coor, h_coor, m_border, j_coor );
		fprintf( stdout, "---> n_next_to_border_changed = %d, m_border_changed = %d\n", n_next_to_border_changed, m_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );
		
		*/


		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_border, n_next_to_border ); // WARNING: note the index inversion

		K_matrix_index1 = COORDINATE_INDEX( i_coor, h_coor );
	      
		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], F_matrix_p->array[ j_coor ], *dummy_matrix1_p ) ) info=1;
	      
		if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
	      
		matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );


		/* dummy_n */
		dummy_n = CMPLX( 
				0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ h_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ h_mode ] +1.0e0 ) )
				);

		
		/* dummy_m */
		dummy_m = CMPLX( 
				-OSQRT2 *( ar_p->rvector[ j_coor] )
				*sqrt( dummy_rho_index2.ivector[ j_mode ] +1.0e0 )
				);


		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

		
		// breaking instruction
		if( info ) break;

	    
		// partial updating
		potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!


	      } /* end conditional */

	  
	    } /* j_mode loop */
	    
	    // breaking instruction
	    if( info ) break;
	  
	  } /* h_mode loop */

	  // breaking instruction
	  if( info ) break;

	  } /* i_mode loop */

	  // breaking instruction
	  if( info ) break;
      
      } /* m_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [7C] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      17TH TERM: 7D
    */

   /*
    
    fprintf( stdout, "--> 17TH TERM: [7D]\n" );
    fprintf( stdout, "\n" );
    
   */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( h_mode=0; h_mode<N_coor_red; h_mode++ ){
	    
	    h_coor = relevant_modes.ivector[ h_mode ]; // i_mode to cartesian coordinate
	  
	    for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	      
	      j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	      
	      //
	      n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( n_next_to_border_extra, i_mode, h_mode );
	      
	      m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( m_next_to_border_extra, j_mode );
	      
	    
	      /* conditional */
	      if( m_next_to_border_changed == n_next_to_border_changed && h_coor != i_coor ){
		
		/*

		fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, h_coor = %d, m_next_to_border = %d, j_coor = %d\n", 
			 n_next_to_border, i_coor, h_coor, m_next_to_border, j_coor );
		fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
		fprintf( stdout, "\n" );
		fprintf( stdout, "dummy_rho_indices\n" );
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		fprintf( stdout, "\n" );
		
		*/


		/* matrix contribution */
		i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion
		
		K_matrix_index1 = COORDINATE_INDEX( i_coor, h_coor );
		
		K_matrix_index2 = COORDINATE_INDEX( j_coor, j_coor );

		if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;
		
		if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
		
		matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );


		/* dummy_n */
		dummy_n = CMPLX( 
				0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ h_coor ] )
				*sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ h_mode ] +1.0e0 ) )
				);
		
		
		/* dummy_m */
		dummy_m = CMPLX( 
				0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ j_coor] )
				*sqrt( ( dummy_rho_index2.ivector[ j_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) )
				);


		/* dummy_tot */
		dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
	      
		dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting

		dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting

	      
		// breaking instruction
		if( info ) break;

	    
		// partial updating
		potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!
	      
	  
	      } /* end conditional */


	    } /* j_mode loop */

	    // breaking instruction
	    if( info ) break;
	  
	  } /* h_mode loop */

	  // breaking instruction
	  if( info ) break;
	  
	} /* i_mode loop */
	
	  // breaking instruction
	if( info ) break;
	
      } /* m_next_to_border_index loop */
      
      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [7D] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------

  if( !info ){

    /*
      18TH TERM: 7E
    */

    /*
    
    fprintf( stdout, "--> 18TH TERM: [7E]\n" );
    fprintf( stdout, "\n" );
    
    */
    

    // set to zero
    potential_energy_corrections_partial = 0.0e0;


    for( n_next_to_border_index=0; n_next_to_border_index<rho_index_next_to_border_length; n_next_to_border_index++ ){
      
      for( m_next_to_border_index=0; m_next_to_border_index<rho_index_next_to_border_length; m_next_to_border_index++ ){
	
	n_next_to_border = rho_index_next_to_border.ivector[ n_next_to_border_index ];

	m_next_to_border = rho_index_next_to_border.ivector[ m_next_to_border_index ];

	// rho_extra
	n_next_to_border_extra = rho_index_aux_translate.ivector[ n_next_to_border ];

	m_next_to_border_extra = rho_index_aux_translate.ivector[ m_next_to_border ];


	// rho_index_inverse_aux
	if( RHO_INDEX_INVERSE_AUX( n_next_to_border, dummy_rho_index1 ) ) info=1;

	if( RHO_INDEX_INVERSE_AUX( m_next_to_border, dummy_rho_index2 ) ) info=1;


	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  i_coor = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate
	  
	  for( h_mode=0; h_mode<N_coor_red; h_mode++ ){
	  
	    h_coor = relevant_modes.ivector[ h_mode ]; // i_mode to cartesian coordinate
	  
	    for( j_mode=0; j_mode<N_coor_red; j_mode++ ){
	      
	      j_coor = relevant_modes.ivector[ j_mode ]; // i_mode to cartesian coordinate
	      
	      for( k_mode=0; k_mode<N_coor_red; k_mode++ ){
		
		k_coor = relevant_modes.ivector[ k_mode ]; // i_mode to cartesian coordinate
		
	    
		//
		n_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( n_next_to_border_extra, i_mode, h_mode );
		
		m_next_to_border_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( m_next_to_border_extra, j_mode, k_mode );
		
	      
		/* conditional */
		if( m_next_to_border_changed == n_next_to_border_changed && h_coor != i_coor && k_coor != j_coor ){
		  
		  /*
		
		  fprintf( stdout, "---> n_next_to_border = %d, i_coor = %d, h_coor = %d, m_next_to_border = %d, j_coor = %d, k_coor = %d\n", 
			   n_next_to_border, i_coor, h_coor, m_next_to_border, j_coor, k_coor );
		  fprintf( stdout, "---> n_next_to_border_changed = %d, m_next_to_border_changed = %d\n", n_next_to_border_changed, m_next_to_border_changed );
		  fprintf( stdout, "\n" );
		  fprintf( stdout, "dummy_rho_indices\n" );
		  if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
		  if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
		  fprintf( stdout, "\n" );
		  
		  */
		
		
		  /* matrix contribution */
		  i_rho = RHO_INDEX_AUX( m_next_to_border, n_next_to_border ); // WARNING: note the index inversion
		  
		  K_matrix_index1 = COORDINATE_INDEX( i_coor, h_coor );
		  
		  K_matrix_index2 = COORDINATE_INDEX( j_coor, k_coor );
		  
		  if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index1 ], K_matrix_p->array[ K_matrix_index2 ], *dummy_matrix1_p ) ) info=1;
		  
		  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, rho_p->array[ i_rho ], *dummy_matrix2_p ) ) info=1;
		  
		  matrix_contribution = MATRIX_TRACE( *dummy_matrix2_p );
		
		  
		  /* dummy_n */
		  dummy_n = CMPLX( 
				  0.25e0 *( ar_p->rvector[ i_coor ] ) *( ar_p->rvector[ h_coor ] )
				  *sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ h_mode ] +1.0e0 ) )
				  );
		  
		
		  /* dummy_m */
		  dummy_m = CMPLX( 
				  0.25e0 *( ar_p->rvector[ j_coor] ) *( ar_p->rvector[ k_coor] )
				  *sqrt( ( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ k_mode ] +1.0e0 ) )
				  );
		  
		  
		  /* dummy_tot */
		  dummy_tot = CMPLX_PRODUCT( dummy_n, dummy_m );
		  
		  dummy_tot = CMPLX_PRODUCT( dummy_tot, matrix_contribution ); //WARNING: overwriting
		  
		  dummy_tot = CMPLX_DIVISION( dummy_tot, IHBAR ); //WARNING: overwriting
		
		
		  // breaking instruction
		  if( info ) break;
	      
	    
		  // partial updating
		  potential_energy_corrections_partial -= 2.0e0 *REAL( dummy_tot ); //WARNING: the sign is MINUS!

	      
		} /* end conditional */

	  
	      } /* k_mode loop */

	      // breaking instruction
	      if( info ) break;
	    
	    } /* j_mode loop */

	    // breaking instruction
	    if( info ) break;

	  } /* h_mode loop */

	  // breaking instruction
	  if( info ) break;

	} /* i_mode loop */

	// breaking instruction
	if( info ) break;
      
      } /* m_next_to_border_index loop */

      // breaking instruction
      if( info ) break;
    
    } /* n_next_to_border_index loop */
  
  } /* end info conditional */

  if( !info ){


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "--------------\n" );
    fprintf( stdout, "Term [7E] partial correction is %le\n", potential_energy_corrections_partial );
    fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


    // updating
    potential_energy_corrections += potential_energy_corrections_partial;

  }

  //----------


  /* corrections due to the Ehrenfest frame transform */
  if( !flag_no_Ehrenfest_frame ){

    if( !info ){

      /*
        19TH TERM: Ehrenfest corrections I
      */

      /*
    
      fprintf( stdout, "--> 19TH TERM: [Ehrenfest I]\n" );
      fprintf( stdout, "\n" );
    
      */
    
    
      // set to zero
      potential_energy_corrections_partial = 0.0e0;


      for( i_coor=0; i_coor<N_coor; i_coor++ ){

        if( MATRIX_COMMUTATOR( F_matrix_p->array[ i_coor ], *H_matrix_p, *dummy_matrix1_p ) ) info=1;

        if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, mu10_p->array[ i_coor ], *dummy_matrix2_p ) ) info=1;

        potential_energy_corrections_partial -= IMAG( MATRIX_TRACE( *dummy_matrix2_p ) );

      }

      potential_energy_corrections_partial /= HBAR; 

    } /* end info conditional */


    if( !info ){

#ifdef __DEBUG_PLUS__

      fprintf( stdout, "--------------\n" );
      fprintf( stdout, "Term [Ehrenfest I] partial correction is %le\n", potential_energy_corrections_partial );
      fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


      // updating
      //    potential_energy_corrections += potential_energy_corrections_partial;

    }

    //----------

    if( !info ){

      /*
        20TH TERM: Ehrenfest corrections II
      */

     /*
      
      fprintf( stdout, "--> 20TH TERM: [Ehrenfest II]\n" );
      fprintf( stdout, "\n" );
    
      */
    
    
      // set to zero
      potential_energy_corrections_partial = 0.0e0;

      for( i_coor=0; i_coor<N_coor; i_coor++ ){

        for( j_coor=0; j_coor<N_coor; j_coor++ ){

	  K_matrix_index1 = COORDINATE_INDEX( i_coor, j_coor );

	  if( MATRIX_COMMUTATOR( K_matrix_p->array[ K_matrix_index1 ], *H_matrix_p, *dummy_matrix1_p ) ) info=1;
	
	  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, mu20_p->array[ K_matrix_index1 ], *dummy_matrix2_p ) ) info=1;

          potential_energy_corrections_partial += IMAG( MATRIX_TRACE( *dummy_matrix2_p ) );

        } /* end j_coor loop */

      } /* end i_coor loop */

      potential_energy_corrections_partial /= 2.0e0 *HBAR; 

    } /* end info conditional */


    if( !info ){

#ifdef __DEBUG_PLUS__

      fprintf( stdout, "--------------\n" );
      fprintf( stdout, "Term [Ehrenfest II] partial correction is %le\n", potential_energy_corrections_partial );
      fprintf( stdout, "--------------\n" );
    
#endif /* __DEBUG_PLUS__ */


      // updating
      //    potential_energy_corrections += potential_energy_corrections_partial;

    }

   //----------

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "--------------\n" );
  fprintf( stdout, "total potential correction is %le\n", potential_energy_corrections );
  fprintf( stdout, "--------------\n" );
  
#endif /* __DEBUG_PLUS__ */


  // update
  *potential_energy_system_correction_p += 0.5e0 *( potential_energy_corrections_old +potential_energy_corrections ) *dt;
  

  /*

    fprintf( stdout, "# potential_energy_corrections_old is %le\n", potential_energy_corrections_old );
    fprintf( stdout, "# potential_energy_corrections     is %le\n", potential_energy_corrections     );
    
  */
  

  /* copy back */
  if( flag_normal_mode_expansion ){

    if( RVECTOR_COPY( *momenta_p, state_p->momenta_saved ) ) info=1;

    if( MATRIX_ARRAY_COPY( *F_matrix_p, state_p->F_matrix_saved ) ) info=1;

    if( MATRIX_ARRAY_COPY( *delta_F_matrix_p, state_p->delta_F_matrix_saved ) ) info=1;

    if( MATRIX_ARRAY_COPY( *K_matrix_p, state_p->K_matrix_saved ) ) info=1;

  }


  return info;

}

//------------------------------------------
