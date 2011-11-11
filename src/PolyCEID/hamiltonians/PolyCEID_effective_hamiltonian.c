
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
#include "PolyCEID_effective_hamiltonian.h"


/* privates */
vector    delta_F_matrix_reduced;
matrix    K_matrix_reduced;
ivector   dummy_index_tmp1;
ivector   dummy_index_tmp2;


/*********************
  FUNCTIONS & MACROS
*********************/

/* effective_hamiltonian_update */

//BUGFIX: this was originally TMP
int effective_hamiltonian_update( const constants constants, state_p state_p, config_p config_p, matrix_p effective_hamiltonian_p ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  int             N_coor;
  int             N_coor_red;
  int             max_rho_index;
  ivector         relevant_modes;
  rvector         masses_aux;
  /* state */
  matrix_p        initial_rho_electron_p;
  matrix_p        H_matrix_p;
  matrix_array_p  delta_F_matrix_p;
  matrix_array_p  K_matrix_p;
  matrix_p        dummy_matrix1_p;
  rvector_p       ar_p;
  rvector_p       ap_p;
  /* dummies */
  int             i, j;
  int             index;
  int             i_mode;
  int             j_mode;
  int             i_coor;
  int             j_coor;
  int             K_matrix_index;
  int             info=0;
  complex         H_matrix_reduced;
  complex         dummy;
  complex         dummy_sum;
  complex         dummy_mode;
#ifdef __DEBUG_PLUS__
  double          norm;
#endif /* __DEBUG_PLUS__ */


  flag_normal_mode_expansion =  constants.flag_normal_mode_expansion;
  N_coor                     =  constants.N_coor;
  N_coor_red                 =  constants.N_coor_red;
  max_rho_index              =  constants.max_rho_index;
  relevant_modes             =  constants.relevant_modes;
  masses_aux                 =  config_p->atoms.masses_aux;

  initial_rho_electron_p     = &(state_p->initial_rho_electron);
  H_matrix_p                 = &(config_p->electrons.H_matrix);
  delta_F_matrix_p           = &(config_p->electrons.delta_F_matrix);
  K_matrix_p                 = &(config_p->electrons.K_matrix);
  dummy_matrix1_p            = &(state_p->dummy_matrix1);
  ar_p                       = &(state_p->phonons.ar);
  ap_p                       = &(state_p->phonons.ap);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: effective_hamiltonian_update\n" );

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    // BUGFIX:config?
    /* compute initial_rho_electron */
    if( COMPUTE_INITIAL_RHO_ELECTRON( constants, *state_p) ) info=1;
      
  }


  // WARNING: mu are needed by force calculation

  if( !info ){

    /* update mu */
    if( MU_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }


  if( !info && flag_normal_mode_expansion ){

    /* transform mu */
    if( TRANSFORM_MU( constants, *state_p, *config_p ) ) info=1;    //WARNING: is it really needed?

  }

  
  if( !info ){

    /* update distance */
    if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;
    
  }


  if( !info ){
    
    /* initialise Ehrenfest frame */
    if( INITIALISE_EHRENFEST_FRAME( constants, *state_p, *config_p ) ) info=1;
    
  }


  if( !info ){

    /* update hamiltonian */
    if( HAMILTONIAN_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;
    
  }


  if( !info && flag_normal_mode_expansion ){

    /* transform hamiltonian */
    if( TRANSFORM_HAMILTONIAN( constants, *state_p, *config_p, config_p->electrons.delta_F_matrix, config_p->electrons.K_matrix ) ) info=1;    //WARNING: is it really needed?

  }


  // H_matrix_reduced
  if( MATRIX_MATRIX_PRODUCT( *H_matrix_p, *initial_rho_electron_p, *dummy_matrix1_p ) ) info=1;
    
  H_matrix_reduced = MATRIX_TRACE( *dummy_matrix1_p );
  
  /*
    fprintf( stdout, "H_matrix_reduced\n" );
    if( CMPLX_PRINT_PLUS( stdout, H_matrix_reduced ) ) info=1;
    fprintf( stdout, "\n\n" );
  */

  /*
    fprintf( stdout, "initial_rho_electron\n" );
    if( MATRIX_PRINT_PLUS( stdout, *initial_rho_electron_p ) ) info=1;
    fprintf( stdout, "\n\n" );

    fprintf( stdout, "delta_F_matrix\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *delta_F_matrix_p ) ) info=1;
    fprintf( stdout, "\n\n" );
  */


  // delta_F_matrix_reduced
  for( i=0; i<N_coor; i++ ){

    if( MATRIX_MATRIX_PRODUCT( delta_F_matrix_p->array[ i ], *initial_rho_electron_p, *dummy_matrix1_p ) ) info=1;
      
    delta_F_matrix_reduced.vector[ i ] = MATRIX_TRACE( *dummy_matrix1_p );

  } /* i loop */

  /*
    fprintf( stdout, "delta_F_matrix_reduced\n" );
    if( VECTOR_PRINT_PLUS( stdout, delta_F_matrix_reduced ) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "forces\n" );
    if( RVECTOR_PRINT_PLUS( stdout, state_p->atoms.forces ) ) info=1;
    fprintf( stdout, "\n" );
  */


  // K_matrix_reduced
  for( i=0; i<N_coor; i++ ){

    for( j=0; j<N_coor; j++ ){

      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ COORDINATE_INDEX( i, j ) ], *initial_rho_electron_p, *dummy_matrix1_p ) ) info=1;
      
      K_matrix_reduced.matrix[ COORDINATE_INDEX( i, j ) ] = MATRIX_TRACE( *dummy_matrix1_p );

    } /* j loop */

  } /* i loop */


#ifndef __NO_HESSIAN_CORRECTIONS__

  /*
    fprintf( stdout, "K_matrix_reduced [before corrections]\n" );
    if( MATRIX_PRINT_PLUS( stdout, K_matrix_reduced ) ) info=1;
    fprintf( stdout, "\n" );
  */


  if( flag_normal_mode_expansion ){

    /* transform F_matrix as well*/
    if( TRANSFORM_F_MATRIX( constants, *state_p, *config_p, config_p->electrons.F_matrix ) ) info=1;

  }


  /* correction to the bare K_matrix_reduced */
  if( HESSIAN_CORRECTIONS( constants, *state_p, *config_p, K_matrix_reduced ) ) info=1;


#endif /* __NO_HESSIAN_CORRECTIONS__ */

  /*
    fprintf( stdout, "K_matrix_reduced\n" );
    if( MATRIX_PRINT_PLUS( stdout, K_matrix_reduced ) ) info=1;
    fprintf( stdout, "\n" );
  */

  //----------
  /* updating */

  if( !info ){

    /* main loop */
    for( index=0; index<max_rho_index; index++ ){

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "->index = %d\n", index );
#endif /* __DEBUG_PLUS__ */


      /* indeces initilisation */
      if( RHO_INDEX_INVERSE( index, dummy_rho_index1, dummy_rho_index2 ) ) info=1;

      /*
	fprintf( stdout, "  dummy_rho_indices\n" );
	if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
      */

      if( IVECTOR_COPY( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

      if( IVECTOR_COPY( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

      /*
	fprintf( stdout, "  dummy_rho_indices\n" );
	if( IVECTOR_PRINT_PLUS( stdout, dummy_index_tmp1 ) ) info=1;
	if( IVECTOR_PRINT_PLUS( stdout, dummy_index_tmp2 ) ) info=1;
      */

      /* set to zero */
      effective_hamiltonian_p->matrix[ index ] = CMPLX_ZERO;

      /*----------------------
	First part
	----------------------*/

      //      fprintf( stdout, "  ---------- First part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;
      
      for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

	i_coor = relevant_modes.ivector[ i_mode ];
	  
	//	fprintf( stdout, "   i_coor = %d, i_atom = %d \n", i_coor, i_atom );

	/* set to zero */
	dummy_sum = CMPLX_ZERO;

	// 1st term

	dummy_index_tmp1.ivector[ i_mode ] -= 2;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] -1 ) *( dummy_rho_index1.ivector[ i_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp1.ivector[ i_mode ] += 2;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// 2nd term

	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( 2 *dummy_rho_index1.ivector[ i_mode ] +1 );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_DIF( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */


	// 3rd term

	dummy_index_tmp2.ivector[ i_mode ] -= 2;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] -1 ) *( dummy_rho_index2.ivector[ i_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp2.ivector[ i_mode ] += 2;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// final

	dummy = CMPLX( ( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] ) /( masses_aux.rvector[ i_coor ] ) );

	dummy_sum = CMPLX_PRODUCT( dummy, dummy_sum ); //WARNING: overwtiting

	//
	dummy_mode = CMPLX_SUM( dummy_mode, dummy_sum ); //WARNING: overwtiting

    } /* i_mode loop */

      // global

      dummy = CMPLX( 0.25e0 );

      dummy_mode = CMPLX_PRODUCT( dummy, dummy_mode ); //WARNING: overwtiting

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [1] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

      effective_hamiltonian_p->matrix[ index ] = CMPLX_DIF( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting

      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /*----------------------
	Second part
	----------------------*/

      //      fprintf( stdout, "  ---------- Second part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;

      if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  //	  fprintf( stdout, "-->i_mode = %d\n", i_mode );
	  	  
	  dummy = H_matrix_reduced;
	  
	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  dummy_mode = CMPLX_SUM( dummy_mode, dummy ); //WARNING: overwtiting

	} /* i_mode loop */

	// final

	dummy_mode = CMPLX_PRODUCT( CMPLX( 0.5e0 ), dummy_mode ); //WARNING: overwtiting

      }

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [2] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */
      
      effective_hamiltonian_p->matrix[ index ] = CMPLX_SUM( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting
	
      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */


      /*----------------------
	Third part
	----------------------*/

      //      fprintf( stdout, "  ---------- Third part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;

      for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

	i_coor = relevant_modes.ivector[ i_mode ];
	  
	//	fprintf( stdout, "   i_coor = %d \n", i_coor );

	/* set to zero */
	dummy_sum = CMPLX_ZERO;

	// 1st term

	dummy_index_tmp1.ivector[ i_mode ] -= 1;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode ] ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp1.ivector[ i_mode ] += 1;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	// 2nd term

	dummy_index_tmp2.ivector[ i_mode ] -= 1;

	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( dummy_rho_index2.ivector[ i_mode ] ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp2.ivector[ i_mode ] += 1;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// final

	dummy = CMPLX_PRODUCT( ( delta_F_matrix_reduced.vector[ i_coor ] ), CMPLX( ar_p->rvector[ i_coor ] ) );

	dummy_sum = CMPLX_PRODUCT( dummy, dummy_sum ); //WARNING: overwtiting

	//
	dummy_mode = CMPLX_SUM( dummy_mode, dummy_sum ); //WARNING: overwtiting

    } /* i_mode loop */

      // global

      dummy = CMPLX( sqrt( 0.5e0 ) );

      dummy_mode = CMPLX_PRODUCT( dummy, dummy_mode ); //WARNING: overwtiting

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [3] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

      effective_hamiltonian_p->matrix[ index ] = CMPLX_DIF( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting

      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /*----------------------
	fourth part
	----------------------*/

      //      fprintf( stdout, "  ---------- Fourth part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;

      for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

	i_coor = relevant_modes.ivector[ i_mode ];
	  
	K_matrix_index = COORDINATE_INDEX( i_coor, i_coor );

	//	fprintf( stdout, "   i_coor = %d \n", i_coor );

	//	fprintf( stdout, "   K_matrix_index = %d \n", K_matrix_index );

	/* set to zero */
	dummy_sum = CMPLX_ZERO;

	// 1st term

	dummy_index_tmp1.ivector[ i_mode ] -= 2;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] -1 ) *( dummy_rho_index1.ivector[ i_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp1.ivector[ i_mode ] += 2;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// 2nd term

	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( 2 *dummy_rho_index1.ivector[ i_mode ] +1 );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */


	// 3rd term

	dummy_index_tmp2.ivector[ i_mode ] -= 2;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] -1 ) *( dummy_rho_index2.ivector[ i_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp2.ivector[ i_mode ] += 2;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// final

	dummy = CMPLX_PRODUCT( K_matrix_reduced.matrix[ K_matrix_index ], CMPLX( ar_p->rvector[ i_coor ] *ar_p->rvector[ i_coor ] ) );

	dummy_sum = CMPLX_PRODUCT( dummy, dummy_sum ); //WARNING: overwtiting

	//
	dummy_mode = CMPLX_SUM( dummy_mode, dummy_sum ); //WARNING: overwtiting

    } /* i_mode loop */

      // global

      dummy = CMPLX( 0.25e0 );

      dummy_mode = CMPLX_PRODUCT( dummy, dummy_mode ); //WARNING: overwtiting

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [4] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

      effective_hamiltonian_p->matrix[ index ] = CMPLX_SUM( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting

      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /*----------------------
	Fifth part
	----------------------*/

      //      fprintf( stdout, "  ---------- Fifth part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;

      for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

	for( j_mode=i_mode+1; j_mode<N_coor_red; j_mode++ ){

	  i_coor = relevant_modes.ivector[ i_mode ];
	  
	  j_coor = relevant_modes.ivector[ j_mode ];

	  K_matrix_index = COORDINATE_INDEX( i_coor, j_coor );
	  
	 
	  //	    fprintf( stdout, "   i_coor = %d, j_coor = %d \n", i_coor, j_coor );
	    
	  //	    fprintf( stdout, "   K_matrix_index = %d \n", K_matrix_index );
	

	  /* set to zero */
	  dummy_sum = CMPLX_ZERO;
	  
	  // 1st term
	  
	  dummy_index_tmp1.ivector[ i_mode ] -= 1;
	  
	  dummy_index_tmp1.ivector[ j_mode ] -= 1;

	  if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] ) *( dummy_rho_index1.ivector[ j_mode ] ) ) );
	    
	    /*
	      fprintf( stdout, "dummy\n" );
	      if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	      fprintf( stdout, "\n" );
	    */

	    dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	  } /* compare if */

	  dummy_index_tmp1.ivector[ i_mode ] += 1;
	  
	  dummy_index_tmp1.ivector[ j_mode ] += 1;
	  
#ifdef __DEBUG__
	  
	  if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;
	
	  if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	  // 2nd term
	  
	  dummy_index_tmp1.ivector[ i_mode ] -= 1;
	  
	  dummy_index_tmp2.ivector[ j_mode ] -= 1;

	  if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] ) *( dummy_rho_index2.ivector[ j_mode ] ) ) );
	    
	    /*
	      fprintf( stdout, "dummy\n" );
	      if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	      fprintf( stdout, "\n" );
	    */

	    dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	  } /* compare if */

	  dummy_index_tmp1.ivector[ i_mode ] += 1;
	  
	  dummy_index_tmp2.ivector[ j_mode ] += 1;
	  
#ifdef __DEBUG__
	  
	  if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	  if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	  // 3rd term
	  
	  dummy_index_tmp1.ivector[ j_mode ] -= 1;
	  
	  dummy_index_tmp2.ivector[ i_mode ] -= 1;

	  if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ j_mode ] ) *( dummy_rho_index2.ivector[ i_mode ] ) ) );
	    
	    /*
	      fprintf( stdout, "dummy\n" );
	      if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	      fprintf( stdout, "\n" );
	    */

	    dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	  } /* compare if */

	  dummy_index_tmp1.ivector[ j_mode ] += 1;
	  
	  dummy_index_tmp2.ivector[ i_mode ] += 1;

#ifdef __DEBUG__
	  
	  if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;
	
	  if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	  // 4th term
	  
	  dummy_index_tmp2.ivector[ i_mode ] -= 1;
	  
	  dummy_index_tmp2.ivector[ j_mode ] -= 1;

	  if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] ) *( dummy_rho_index2.ivector[ j_mode ] ) ) );
	    
	    /*
	      fprintf( stdout, "dummy\n" );
	      if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	      fprintf( stdout, "\n" );
	    */

	    dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	  } /* compare if */

	  dummy_index_tmp2.ivector[ i_mode ] += 1;
	  
	  dummy_index_tmp2.ivector[ j_mode ] += 1;
	  
#ifdef __DEBUG__
	  
	  if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	  if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	  // final
	  dummy = CMPLX_PRODUCT( K_matrix_reduced.matrix[ K_matrix_index ], CMPLX( ar_p->rvector[ i_coor ] *ar_p->rvector[ j_coor ] ) );

	  dummy_sum = CMPLX_PRODUCT( dummy, dummy_sum ); //WARNING: overwtiting

	  //
	  dummy_mode = CMPLX_SUM( dummy_mode, dummy_sum ); //WARNING: overwtiting

	} /* j_mode loop */

      } /* i_mode loop */

      // global

      dummy = CMPLX( 0.5e0 );

      dummy_mode = CMPLX_PRODUCT( dummy, dummy_mode ); //WARNING: overwtiting

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [5] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

      effective_hamiltonian_p->matrix[ index ] = CMPLX_SUM( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting

      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */


#ifdef __DEBUG_PLUS__
      fprintf( stdout, "\n" );
      fprintf( stdout, "-----------\n" );
      fprintf( stdout, "effective_hamiltonian[ %d ]\n", index );
      if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
      fprintf( stdout, "\n" );
      fprintf( stdout, "-----------\n" );
      fprintf( stdout, "\n" );
#endif /* __DEBUG_PLUS__ */


      /* conditioning */

      //      if( CMPLX_NORM(  effective_hamiltonian_p->matrix[ index ] ) < EPS_LOC )  effective_hamiltonian_p->matrix[ index ] = CMPLX_ZERO;

      
    } /* index loop */

  }

  //----------


  return info;
    
}

//------------------------------------------

/* effective_hamiltonian_update_aux */

// BUGFIX: this was originally DEF
int effective_hamiltonian_update_aux( const constants constants, state_p state_p, config_p config_p, vector electronic_state, matrix_p effective_hamiltonian_p ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  int             N_coor;
  int             N_coor_red;
  int             max_rho_index;
  ivector         relevant_modes;
  rvector         masses_aux;
  /* state */
  matrix_p        mu00_p;
  matrix_p        H_matrix_p;
  matrix_array_p  delta_F_matrix_p;
  matrix_array_p  K_matrix_p;
  matrix_p        dummy_matrix1_p;
  rvector_p       ar_p;
  rvector_p       ap_p;
  /* dummies */
  int             i, j;
  int             index;
  int             i_mode;
  int             j_mode;
  int             i_coor;
  int             j_coor;
  int             K_matrix_index;
  int             info=0;
  complex         H_matrix_reduced;
  complex         dummy;
  complex         dummy_sum;
  complex         dummy_mode;
#ifdef __DEBUG_PLUS__
  double          norm;
#endif /* __DEBUG_PLUS__ */


  flag_normal_mode_expansion =  constants.flag_normal_mode_expansion;
  N_coor                     =  constants.N_coor;
  N_coor_red                 =  constants.N_coor_red;
  max_rho_index              =  constants.max_rho_index;
  relevant_modes             =  constants.relevant_modes;
  masses_aux                 =  config_p->atoms.masses_aux;

  mu00_p                     = &(config_p->electrons.mu00);
  H_matrix_p                 = &(config_p->electrons.H_matrix);
  delta_F_matrix_p           = &(config_p->electrons.delta_F_matrix);
  K_matrix_p                 = &(config_p->electrons.K_matrix);
  dummy_matrix1_p            = &(state_p->dummy_matrix1);
  ar_p                       = &(state_p->phonons.ar);
  ap_p                       = &(state_p->phonons.ap);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: effective_hamiltonian_update\n" );

#endif /* __DEBUG_PLUS__ */



  // WARNING: mu are needed by force calculation

  if( !info ){

    /* update mu */
    if( MU_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }


  if( !info && flag_normal_mode_expansion ){

    /* transform mu */
    if( TRANSFORM_MU( constants, *state_p, *config_p ) ) info=1;    //WARNING: is it really needed?

  }

  
  if( !info ){

    /* update distance */
    if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;
    
  }


  if( !info ){
    
    /* initialise Ehrenfest frame */
    if( INITIALISE_EHRENFEST_FRAME( constants, *state_p, *config_p ) ) info=1;
    
  }


  if( !info ){

    /* update hamiltonian */
    if( HAMILTONIAN_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;
    
  }


  if( !info && flag_normal_mode_expansion ){

    /* transform hamiltonian */
    if( TRANSFORM_HAMILTONIAN( constants, *state_p, *config_p, config_p->electrons.delta_F_matrix, config_p->electrons.K_matrix ) ) info=1;    //WARNING: is it really needed?

  }


  // H_matrix_reduced
  if( MATRIX_MATRIX_PRODUCT( *H_matrix_p, *mu00_p, *dummy_matrix1_p ) ) info=1;
    
  H_matrix_reduced = MATRIX_TRACE( *dummy_matrix1_p );
  
  /*
    fprintf( stdout, "H_matrix_reduced [AUX]\n" );
    if( CMPLX_PRINT_PLUS( stdout, H_matrix_reduced ) ) info=1;
    fprintf( stdout, "\n\n" );
  */

  /*
    fprintf( stdout, "mu00 [AUX]\n" );
    if( MATRIX_PRINT_PLUS( stdout, *mu00_p ) ) info=1;
    fprintf( stdout, "\n\n" );

    fprintf( stdout, "delta_F_matrix [AUX]\n" );
    if( MATRIX_ARRAY_PRINT_PLUS( stdout, *delta_F_matrix_p ) ) info=1;
    fprintf( stdout, "\n\n" );
  */


  // delta_F_matrix_reduced
  for( i=0; i<N_coor; i++ ){

    if( MATRIX_MATRIX_PRODUCT( delta_F_matrix_p->array[ i ], *mu00_p, *dummy_matrix1_p ) ) info=1;
      
    delta_F_matrix_reduced.vector[ i ] = MATRIX_TRACE( *dummy_matrix1_p );

  } /* i loop */

  /*
    fprintf( stdout, "delta_F_matrix_reduced [AUX]\n" );
    if( VECTOR_PRINT_PLUS( stdout, delta_F_matrix_reduced ) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "forces [AUX]\n" );
    if( RVECTOR_PRINT_PLUS( stdout, state_p->atoms.forces ) ) info=1;
    fprintf( stdout, "\n" );
  */


  // K_matrix_reduced
  for( i=0; i<N_coor; i++ ){

    for( j=0; j<N_coor; j++ ){

      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ COORDINATE_INDEX( i, j ) ], *mu00_p, *dummy_matrix1_p ) ) info=1;
      
      K_matrix_reduced.matrix[ COORDINATE_INDEX( i, j ) ] = MATRIX_TRACE( *dummy_matrix1_p );

    } /* j loop */

  } /* i loop */


#ifndef __NO_HESSIAN_CORRECTIONS__

  /*
    fprintf( stdout, "K_matrix_reduced [AUX, before corrections]\n" );
    if( MATRIX_PRINT_PLUS( stdout, K_matrix_reduced ) ) info=1;
    fprintf( stdout, "\n" );
  */


  if( flag_normal_mode_expansion ){

    /* transform F_matrix as well*/
    if( TRANSFORM_F_MATRIX( constants, *state_p, *config_p, config_p->electrons.F_matrix ) ) info=1;

  }


  /* correction to the bare K_matrix_reduced */
  if( HESSIAN_CORRECTIONS_AUX( constants, *state_p, *config_p, electronic_state, K_matrix_reduced ) ) info=1;


#endif /* __NO_HESSIAN_CORRECTIONS__ */

  /*
    fprintf( stdout, "K_matrix_reduced [AUX]\n" );
    if( MATRIX_PRINT_PLUS( stdout, K_matrix_reduced ) ) info=1;
    fprintf( stdout, "\n" );
  */

  //----------
  /* updating */

  if( !info ){

    /* main loop */
    for( index=0; index<max_rho_index; index++ ){

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "->index = %d\n", index );
#endif /* __DEBUG_PLUS__ */


      /* indeces initilisation */
      if( RHO_INDEX_INVERSE( index, dummy_rho_index1, dummy_rho_index2 ) ) info=1;

      /*
	fprintf( stdout, "  dummy_rho_indices\n" );
	if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
	if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
      */

      if( IVECTOR_COPY( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

      if( IVECTOR_COPY( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

      /*
	fprintf( stdout, "  dummy_rho_indices\n" );
	if( IVECTOR_PRINT_PLUS( stdout, dummy_index_tmp1 ) ) info=1;
	if( IVECTOR_PRINT_PLUS( stdout, dummy_index_tmp2 ) ) info=1;
      */

      /* set to zero */
      effective_hamiltonian_p->matrix[ index ] = CMPLX_ZERO;

      /*----------------------
	First part
	----------------------*/

      //      fprintf( stdout, "  ---------- First part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;
      
      for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

	i_coor = relevant_modes.ivector[ i_mode ];
	  
	//	fprintf( stdout, "   i_coor = %d, i_atom = %d \n", i_coor, i_atom );

	/* set to zero */
	dummy_sum = CMPLX_ZERO;

	// 1st term

	dummy_index_tmp1.ivector[ i_mode ] -= 2;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] -1 ) *( dummy_rho_index1.ivector[ i_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp1.ivector[ i_mode ] += 2;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// 2nd term

	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( 2 *dummy_rho_index1.ivector[ i_mode ] +1 );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_DIF( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */


	// 3rd term

	dummy_index_tmp2.ivector[ i_mode ] -= 2;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] -1 ) *( dummy_rho_index2.ivector[ i_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp2.ivector[ i_mode ] += 2;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// final

	dummy = CMPLX( ( ap_p->rvector[ i_coor ] ) *( ap_p->rvector[ i_coor ] ) /( masses_aux.rvector[ i_coor ] ) );

	dummy_sum = CMPLX_PRODUCT( dummy, dummy_sum ); //WARNING: overwtiting

	//
	dummy_mode = CMPLX_SUM( dummy_mode, dummy_sum ); //WARNING: overwtiting

    } /* i_mode loop */

      // global

      dummy = CMPLX( 0.25e0 );

      dummy_mode = CMPLX_PRODUCT( dummy, dummy_mode ); //WARNING: overwtiting

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [1] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

      effective_hamiltonian_p->matrix[ index ] = CMPLX_DIF( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting

      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /*----------------------
	Second part
	----------------------*/

      //      fprintf( stdout, "  ---------- Second part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;

      if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	for( i_mode=0; i_mode<N_coor_red; i_mode++ ){
	  
	  //	  fprintf( stdout, "-->i_mode = %d\n", i_mode );
	  	  
	  dummy = H_matrix_reduced;
	  
	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  dummy_mode = CMPLX_SUM( dummy_mode, dummy ); //WARNING: overwtiting

	} /* i_mode loop */

	// final

	dummy_mode = CMPLX_PRODUCT( CMPLX( 0.5e0 ), dummy_mode ); //WARNING: overwtiting

      }

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [2] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */
      
      effective_hamiltonian_p->matrix[ index ] = CMPLX_SUM( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting
	
      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */


      /*----------------------
	Third part
	----------------------*/

      //      fprintf( stdout, "  ---------- Third part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;

      for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

	i_coor = relevant_modes.ivector[ i_mode ];
	  
	//	fprintf( stdout, "   i_coor = %d \n", i_coor );

	/* set to zero */
	dummy_sum = CMPLX_ZERO;

	// 1st term

	dummy_index_tmp1.ivector[ i_mode ] -= 1;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode ] ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp1.ivector[ i_mode ] += 1;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	// 2nd term

	dummy_index_tmp2.ivector[ i_mode ] -= 1;

	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( dummy_rho_index2.ivector[ i_mode ] ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp2.ivector[ i_mode ] += 1;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// final

	dummy = CMPLX_PRODUCT( ( delta_F_matrix_reduced.vector[ i_coor ] ), CMPLX( ar_p->rvector[ i_coor ] ) );

	dummy_sum = CMPLX_PRODUCT( dummy, dummy_sum ); //WARNING: overwtiting

	//
	dummy_mode = CMPLX_SUM( dummy_mode, dummy_sum ); //WARNING: overwtiting

    } /* i_mode loop */

      // global

      dummy = CMPLX( sqrt( 0.5e0 ) );

      dummy_mode = CMPLX_PRODUCT( dummy, dummy_mode ); //WARNING: overwtiting

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [3] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

      effective_hamiltonian_p->matrix[ index ] = CMPLX_DIF( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting

      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /*----------------------
	fourth part
	----------------------*/

      //      fprintf( stdout, "  ---------- Fourth part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;

      for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

	i_coor = relevant_modes.ivector[ i_mode ];
	  
	K_matrix_index = COORDINATE_INDEX( i_coor, i_coor );

	//	fprintf( stdout, "   i_coor = %d \n", i_coor );

	//	fprintf( stdout, "   K_matrix_index = %d \n", K_matrix_index );

	/* set to zero */
	dummy_sum = CMPLX_ZERO;

	// 1st term

	dummy_index_tmp1.ivector[ i_mode ] -= 2;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] -1 ) *( dummy_rho_index1.ivector[ i_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp1.ivector[ i_mode ] += 2;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// 2nd term

	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( 2 *dummy_rho_index1.ivector[ i_mode ] +1 );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */


	// 3rd term

	dummy_index_tmp2.ivector[ i_mode ] -= 2;


	if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] -1 ) *( dummy_rho_index2.ivector[ i_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	} /* compare if */

	dummy_index_tmp2.ivector[ i_mode ] += 2;

#ifdef __DEBUG__

	if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */

	// final

	dummy = CMPLX_PRODUCT( K_matrix_reduced.matrix[ K_matrix_index ], CMPLX( ar_p->rvector[ i_coor ] *ar_p->rvector[ i_coor ] ) );

	dummy_sum = CMPLX_PRODUCT( dummy, dummy_sum ); //WARNING: overwtiting

	//
	dummy_mode = CMPLX_SUM( dummy_mode, dummy_sum ); //WARNING: overwtiting

    } /* i_mode loop */

      // global

      dummy = CMPLX( 0.25e0 );

      dummy_mode = CMPLX_PRODUCT( dummy, dummy_mode ); //WARNING: overwtiting

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [4] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

      effective_hamiltonian_p->matrix[ index ] = CMPLX_SUM( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting

      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /*----------------------
	Fifth part
	----------------------*/

      //      fprintf( stdout, "  ---------- Fifth part -----------\n" );

      /* set to zero */
      dummy_mode = CMPLX_ZERO;

      for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

	for( j_mode=i_mode+1; j_mode<N_coor_red; j_mode++ ){

	  i_coor = relevant_modes.ivector[ i_mode ];
	  
	  j_coor = relevant_modes.ivector[ j_mode ];

	  K_matrix_index = COORDINATE_INDEX( i_coor, j_coor );
	  
	 
	  //	    fprintf( stdout, "   i_coor = %d, j_coor = %d \n", i_coor, j_coor );
	    
	  //	    fprintf( stdout, "   K_matrix_index = %d \n", K_matrix_index );
	

	  /* set to zero */
	  dummy_sum = CMPLX_ZERO;
	  
	  // 1st term
	  
	  dummy_index_tmp1.ivector[ i_mode ] -= 1;
	  
	  dummy_index_tmp1.ivector[ j_mode ] -= 1;

	  if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] ) *( dummy_rho_index1.ivector[ j_mode ] ) ) );
	    
	    /*
	      fprintf( stdout, "dummy\n" );
	      if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	      fprintf( stdout, "\n" );
	    */

	    dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	  } /* compare if */

	  dummy_index_tmp1.ivector[ i_mode ] += 1;
	  
	  dummy_index_tmp1.ivector[ j_mode ] += 1;
	  
#ifdef __DEBUG__
	  
	  if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;
	
	  if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	  // 2nd term
	  
	  dummy_index_tmp1.ivector[ i_mode ] -= 1;
	  
	  dummy_index_tmp2.ivector[ j_mode ] -= 1;

	  if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] ) *( dummy_rho_index2.ivector[ j_mode ] ) ) );
	    
	    /*
	      fprintf( stdout, "dummy\n" );
	      if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	      fprintf( stdout, "\n" );
	    */

	    dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	  } /* compare if */

	  dummy_index_tmp1.ivector[ i_mode ] += 1;
	  
	  dummy_index_tmp2.ivector[ j_mode ] += 1;
	  
#ifdef __DEBUG__
	  
	  if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	  if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	  // 3rd term
	  
	  dummy_index_tmp1.ivector[ j_mode ] -= 1;
	  
	  dummy_index_tmp2.ivector[ i_mode ] -= 1;

	  if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ j_mode ] ) *( dummy_rho_index2.ivector[ i_mode ] ) ) );
	    
	    /*
	      fprintf( stdout, "dummy\n" );
	      if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	      fprintf( stdout, "\n" );
	    */

	    dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	  } /* compare if */

	  dummy_index_tmp1.ivector[ j_mode ] += 1;
	  
	  dummy_index_tmp2.ivector[ i_mode ] += 1;

#ifdef __DEBUG__
	  
	  if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;
	
	  if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	  // 4th term
	  
	  dummy_index_tmp2.ivector[ i_mode ] -= 1;
	  
	  dummy_index_tmp2.ivector[ j_mode ] -= 1;

	  if( !IVECTOR_COMPARE_MUTE( dummy_index_tmp1, dummy_index_tmp2 ) ){
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] ) *( dummy_rho_index2.ivector[ j_mode ] ) ) );
	    
	    /*
	      fprintf( stdout, "dummy\n" );
	      if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	      fprintf( stdout, "\n" );
	    */

	    dummy_sum = CMPLX_SUM( dummy_sum, dummy ); //WARNING: overwtiting
	  
	  } /* compare if */

	  dummy_index_tmp2.ivector[ i_mode ] += 1;
	  
	  dummy_index_tmp2.ivector[ j_mode ] += 1;
	  
#ifdef __DEBUG__
	  
	  if( IVECTOR_COMPARE( dummy_index_tmp1, dummy_rho_index1 ) ) info=1;

	  if( IVECTOR_COMPARE( dummy_index_tmp2, dummy_rho_index2 ) ) info=1;

#endif /* __DEBUG__ */


	  // final
	  dummy = CMPLX_PRODUCT( K_matrix_reduced.matrix[ K_matrix_index ], CMPLX( ar_p->rvector[ i_coor ] *ar_p->rvector[ j_coor ] ) );

	  dummy_sum = CMPLX_PRODUCT( dummy, dummy_sum ); //WARNING: overwtiting

	  //
	  dummy_mode = CMPLX_SUM( dummy_mode, dummy_sum ); //WARNING: overwtiting

	} /* j_mode loop */

      } /* i_mode loop */

      // global

      dummy = CMPLX( 0.5e0 );

      dummy_mode = CMPLX_PRODUCT( dummy, dummy_mode ); //WARNING: overwtiting

#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      if( CMPLX_PRINT_PLUS( stdout, dummy_mode ) ) info=1;
      fprintf( stdout, "\n" );
      norm = CMPLX_NORM( dummy_mode );
      fprintf( stdout, "   Term's norm [5] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

      effective_hamiltonian_p->matrix[ index ] = CMPLX_SUM( effective_hamiltonian_p->matrix[ index ], dummy_mode ); //WARNING: overwtiting

      /*
	if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */


#ifdef __DEBUG_PLUS__
      fprintf( stdout, "\n" );
      fprintf( stdout, "-----------\n" );
      fprintf( stdout, "effective_hamiltonian[ %d ]\n", index );
      if( CMPLX_PRINT_PLUS( stdout, effective_hamiltonian_p->matrix[ index ]  ) ) info=1;
      fprintf( stdout, "\n" );
      fprintf( stdout, "-----------\n" );
      fprintf( stdout, "\n" );
#endif /* __DEBUG_PLUS__ */


      /* conditioning */

      //      if( CMPLX_NORM(  effective_hamiltonian_p->matrix[ index ] ) < EPS_LOC )  effective_hamiltonian_p->matrix[ index ] = CMPLX_ZERO;

      
    } /* index loop */

  }

  //----------


  return info;
    
}

//------------------------------------------

/* utilities */

//------------------------------------------

int allocate_effective_hamiltonian_update_workspace( const constants constants ){

  /* constants */
  int N_coor;
  int N_coor_red;
  /* dummies */
  int info=0;


  N_coor     = constants.N_coor;
  N_coor_red = constants.N_coor_red;


  if( !info ){

    /* allocation */
    if( ALLOCATE_COMPUTE_INITIAL_RHO_ELECTRON_WORKSPACE( constants) ) info=1;

    if( VECTOR_ALLOCATE( N_coor, delta_F_matrix_reduced ) ) info=1;
    
    if( MATRIX_ALLOCATE( N_coor, N_coor, K_matrix_reduced ) ) info=1;
    
    if( IVECTOR_ALLOCATE( N_coor_red, dummy_index_tmp1 ) ) info=1;  
    
    if( IVECTOR_ALLOCATE( N_coor_red, dummy_index_tmp2 ) ) info=1;  

  }


  return info;

}

//------------------------------------------

int free_effective_hamiltonian_update_workspace( void ){

  /* dummies */
  int info=0;


  if( !info ){

    /* deallocation */
    if( IVECTOR_FREE( dummy_index_tmp2 ) ) info=1;
    
    if( IVECTOR_FREE( dummy_index_tmp1 ) ) info=1;
    
    if( MATRIX_FREE( K_matrix_reduced ) ) info=1;
    
    if( VECTOR_FREE( delta_F_matrix_reduced ) ) info=1;

    if( FREE_COMPUTE_INITIAL_RHO_ELECTRON_WORKSPACE() ) info=1;

  }


  return info;

}

//------------------------------------------
