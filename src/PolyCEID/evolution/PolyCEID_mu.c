
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
// Sia benedetto S.Patrizio!
#include "PolyCEID_mu.h"



/*********************
  FUNCTIONS & MACROS
*********************/

/* mu_update */

int mu_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  /* dummies */
  int info=0;


  flag_normal_mode_expansion =  constants.flag_normal_mode_expansion;


  if( !info ){

    if( MU00_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }


  if( !info && flag_normal_mode_expansion ){

    if( IRRELEVANT_MU_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }


  if( !info ){

    if( MU01_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }

  if( !info ){

    if( MU10_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }

  if( !info ){

    if( MU02_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }

  if( !info ){
    
    if( MU11_UPDATE( constants, *state_p, *config_p ) ) info=1;
    
  }

  if( !info ){

    if( MU20_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }


  return info;

}

//------------------------------------------
//------------------------------------------

/* mu00_update */

int mu00_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             sqrt_max_rho_index;
  /* state */
  matrix_array_p  rho_p;
  matrix_p        mu00_p;
  /* dummies */
  int             i;
  int             index;
  int             info=0;


  sqrt_max_rho_index =  constants.sqrt_max_rho_index;

  rho_p              = &(config_p->electrons.rho);
  mu00_p             = &(config_p->electrons.mu00);


  /* set mu00 to zero */
  if( MATRIX_ZERO( *mu00_p ) ) info=1;

  if( !info ){

    for( i=0; i<sqrt_max_rho_index; i++ ){

      if( info ) break;

      index = diagonal_rho_indices.ivector[ i ];

      if( MATRIX_SUM( *mu00_p, rho_p->array[ index ], *mu00_p ) ) info=1; //WARNING: overwriting!

    }

  }


  return info;

}

//------------------------------------------

/* mu01_update */

int mu01_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor_red;
  ivector         relevant_modes;
  int             max_rho_index;
  int             sqrt_max_rho_index;
  rvector_p       ap_p;
  /* state */
  matrix_array_p  rho_p;
  matrix_array_p  mu01_p;
  matrix_p   dummy_matrix1_p;
  /* dummies */
  int             i,j;
  int             i_mode;
  int             index, index_aux;
  int             info=0;
  complex         dummy;


  N_coor_red         =  constants.N_coor_red;
  relevant_modes     =  constants.relevant_modes;
  max_rho_index      =  constants.max_rho_index;
  sqrt_max_rho_index =  constants.sqrt_max_rho_index;

  rho_p              = &(config_p->electrons.rho);
  mu01_p             = &(config_p->electrons.mu01);
  dummy_matrix1_p    = &(state_p->dummy_matrix1);
  ap_p               = &(state_p->phonons.ap); 


  for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

    if( info ) break;

    i = relevant_modes.ivector[ i_mode ]; 
    
    /* set mu01[ i ] to zero */
    if( MATRIX_ZERO( mu01_p->array[ i ] ) ) info=1;
    
    for( j=0; j<sqrt_max_rho_index; j++ ){
	
      if( info ) break;
      
      /* set dummy_matrix1 to zero */
      if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;
      
      index = diagonal_rho_indices.ivector[ j ];
	
      if( RHO_INDEX_INVERSE( index, dummy_rho_index1, dummy_rho_index2 ) ) info=1;

#ifdef __DEBUG__

      if( IVECTOR_ARE_EQUAL( dummy_rho_index1, dummy_rho_index2 ) ){
	
	fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	fflush( stderr );
	
	info=1;
	
	break;
	  
      }

#endif /* __DEBUG__ */

      // 1st part;

      index_aux = RHO_INDEX_CHANGE1_MINUS1_FIRST( index, i_mode );
	
      //      fprintf( stdout, "index_aux = %d\n", index_aux );

      if( index_aux < max_rho_index ){

	dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode ] ) ); // Why not by table lookup?
	
	/*
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/
	
	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index_aux ], dummy, *dummy_matrix1_p  ) ) info=1;
	
	if( MATRIX_SUM( mu01_p->array[ i ], *dummy_matrix1_p, mu01_p->array[ i ] ) ) info=1; //WARNING: overwriting!
	
      }

    } /* end j loop */
    
    /* the observable is hermitian */
    if( MATRIX_ADJOINT( mu01_p->array[ i ], *dummy_matrix1_p  ) ) info=1;
      
    if( MATRIX_DIF( mu01_p->array[ i ], *dummy_matrix1_p, mu01_p->array[ i ] ) ) info=1; //WARNING: overwriting!

    /* global factor */
    dummy = CMPLX( OSQRT2 *( ap_p->rvector[ i ] ) );
    
    dummy = CMPLX_PRODUCT( CMPLX_I, dummy ); //WARNING: overwriting!
    
    if( MATRIX_SCALAR_PRODUCT( mu01_p->array[ i ], dummy, mu01_p->array[ i ] ) ) info=1; //WARNING: overwriting!
    
  } /* end i loop */

  
  return info;

}

//------------------------------------------

/* mu10_updated */

int mu10_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor_red;
  ivector         relevant_modes;
  int             max_rho_index;
  int             sqrt_max_rho_index;
  rvector_p       ar_p;
  /* state */
  matrix_array_p  rho_p;
  matrix_array_p  mu10_p;
  matrix_p   dummy_matrix1_p;
  /* dummies */
  int             i,j;
  int             i_mode;
  int             index, index_aux;
  int             info=0;
  complex         dummy;


  N_coor_red         =  constants.N_coor_red;
  relevant_modes     =  constants.relevant_modes;
  max_rho_index      =  constants.max_rho_index;
  sqrt_max_rho_index =  constants.sqrt_max_rho_index;

  rho_p              = &(config_p->electrons.rho);
  mu10_p             = &(config_p->electrons.mu10);
  dummy_matrix1_p    = &(state_p->dummy_matrix1);
  ar_p               = &(state_p->phonons.ar); 


  for( i_mode=0; i_mode<N_coor_red; i_mode++ ){

    if( info ) break;

    i = relevant_modes.ivector[ i_mode ];

    /* set mu10[ i ] to zero */
    if( MATRIX_ZERO( mu10_p->array[ i ] ) ) info=1;

    for( j=0; j<sqrt_max_rho_index; j++ ){
	
      if( info ) break;
      
      /* set dummy_matrix1 to zero */
      if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;
      
      index = diagonal_rho_indices.ivector[ j ];
      
      if( RHO_INDEX_INVERSE( index, dummy_rho_index1, dummy_rho_index2 ) ) info=1;
      
#ifdef __DEBUG__

      if( IVECTOR_ARE_EQUAL( dummy_rho_index1, dummy_rho_index2 ) ){

	fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	fflush( stderr );
	  
	info=1;
	  
	break;
	  
      }
      
#endif /* __DEBUG__ */
	
      // 1st part;
	
      index_aux = RHO_INDEX_CHANGE1_MINUS1_FIRST( index, i_mode );

      //      fprintf( stdout, "index_aux = %d\n", index_aux );
	
      if( index_aux < max_rho_index ){
	  
	dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode ] ) ); // Why not by table lookup?
	  
	/*
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/
	  
	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index_aux ], dummy, *dummy_matrix1_p  ) ) info=1;
	
	if( MATRIX_SUM( mu10_p->array[ i ], *dummy_matrix1_p, mu10_p->array[ i ] ) ) info=1; //WARNING: overwriting!
	  
      }
	
    } /* end j loop */
      
    /* the observable is hermitian */
    if( MATRIX_ADJOINT( mu10_p->array[ i ], *dummy_matrix1_p  ) ) info=1;

    if( MATRIX_SUM( mu10_p->array[ i ], *dummy_matrix1_p, mu10_p->array[ i ] ) ) info=1; //WARNING: overwriting!

    /* global factor */
    dummy = CMPLX( OSQRT2 *( ar_p->rvector[ i ] ) );
	
    if( MATRIX_SCALAR_PRODUCT( mu10_p->array[ i ], dummy, mu10_p->array[ i ] ) ) info=1; //WARNING: overwriting!

  } /* end i loop */


  return info;

}

//------------------------------------------

/* mu02_updated */

int mu02_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor_red;
  ivector         relevant_modes;
  int             max_rho_index;
  int             sqrt_max_rho_index;
  rvector_p       ap_p;
  /* state */
  matrix_array_p  rho_p;
  matrix_array_p  mu02_p;
  matrix_p   dummy_matrix1_p;
  /* dummies */
  int             i1,i2,j;
  int             i_mode1, i_mode2;
  int             index, index_aux;
  int             index_mu02;
  int             info=0;
  complex         dummy;


  N_coor_red         =  constants.N_coor_red;
  relevant_modes     =  constants.relevant_modes;
  max_rho_index      =  constants.max_rho_index;
  sqrt_max_rho_index =  constants.sqrt_max_rho_index;

  rho_p              = &(config_p->electrons.rho);
  mu02_p             = &(config_p->electrons.mu02);
  dummy_matrix1_p    = &(state_p->dummy_matrix1);
  ap_p               = &(state_p->phonons.ap); 


  for( i_mode1=0; i_mode1<N_coor_red; i_mode1++ ){

    if( info ) break;

    i1 = relevant_modes.ivector[ i_mode1 ];

    for( i_mode2=0; i_mode2<N_coor_red; i_mode2++ ){

      if( info ) break;

      i2 = relevant_modes.ivector[ i_mode2 ];
	  
      index_mu02 = COORDINATE_INDEX( i1, i2 );

      //      fprintf( stdout, "index_mu02 = %d\n", index_mu02 );

      /* set mu02[ i ] to zero */
      if( MATRIX_ZERO( mu02_p->array[ index_mu02 ] ) ) info=1;
      
      for( j=0; j<sqrt_max_rho_index; j++ ){

	if( info ) break;
	
	/* set dummy_matrix1 to zero */
	if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;
	
	index = diagonal_rho_indices.ivector[ j ];
	
	//	fprintf( stdout, "diagonal_index = %d\n", index );

	if( RHO_INDEX_INVERSE( index, dummy_rho_index1, dummy_rho_index2 ) ) info=1;

#ifdef __DEBUG__
	    
	if( IVECTOR_ARE_EQUAL( dummy_rho_index1, dummy_rho_index2 ) ){
	      
	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  info=1;
	  
	  break;
	  
	}
	
#endif /* __DEBUG__ */
	
	// 1st part;
	
	//	fprintf( stdout, "first part\n");
	
	if( i1 != i2 ){

	  index_aux = RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( index, i_mode2, i_mode1 ); //WARNING: note the index inversion

	  if( index_aux < max_rho_index ){
	    
	    //	fprintf( stdout, "index_aux = %d\n", index_aux );
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode1 ] ) *( dummy_rho_index1.ivector[ i_mode2 ] +1 ) ) ); // Why not by table lookup?
	    
	  }
	  else{

	    dummy=CMPLX_ZERO;

	  }

	}
	else{

	  index_aux = index;

	  //	fprintf( stdout, "index_aux = %d\n", index_aux );
	    
	  dummy = CMPLX( (double) dummy_rho_index1.ivector[ i_mode1 ] + 0.5e0 ); // Why not by table lookup?
	  
	}

	/*
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/
	
	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index_aux ], dummy, *dummy_matrix1_p  ) ) info=1;
	
	//	  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p) ) info=1;
	
	if( MATRIX_SUM( mu02_p->array[ index_mu02 ], *dummy_matrix1_p, mu02_p->array[ index_mu02 ] ) ) info=1; //WARNING: overwriting!
      

	// 2nd part;

	//	fprintf( stdout, "second part\n");

	index_aux = RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index, i_mode1, i_mode2 );
	
	//       	fprintf( stdout, "index_aux = %d\n", index_aux );
	    
	if( index_aux < max_rho_index ){
	  
	  if( i1 != i2 ){

	    dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode1 ] *dummy_rho_index1.ivector[ i_mode2 ] ) ); // Why not by table lookup?

	  }
	  else{
	    
	    dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode1 ] *( dummy_rho_index1.ivector[ i_mode1 ] -1 ) ) ); // Why not by table lookup?
		
	  }
	      
	  /*
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index_aux ], dummy, *dummy_matrix1_p  ) ) info=1;

	  //	  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p) ) info=1;
	  
	  if( MATRIX_DIF( mu02_p->array[ index_mu02 ], *dummy_matrix1_p, mu02_p->array[ index_mu02 ] ) ) info=1; //WARNING: overwriting!
	    
	}

	//	fprintf( stdout, "-------------\n" );      
	    
      } /* end j loop */
    
      /*
	fprintf( stdout, "index_mu02 = %d\n", index_mu02 );

	fprintf( stdout, "semifinal\n" );

	if( MATRIX_PRINT_PLUS( stdout,  mu02_p->array[ index_mu02 ] ) ) info=1;
	fprintf( stdout, "\n" );
      */
      
      /* the observable is hermitian */
      if( MATRIX_ADJOINT( mu02_p->array[ index_mu02 ], *dummy_matrix1_p  ) ) info=1;

      if( MATRIX_SUM( mu02_p->array[ index_mu02 ], *dummy_matrix1_p, mu02_p->array[ index_mu02 ] ) ) info=1; //WARNING: overwriting!

      /*
	if( MATRIX_PRINT( stdout,  mu02_p->array[ index_mu02 ] ) ) info=1;
	fprintf( stdout, "\n" );
      */
	  
      /* global factor */
      dummy = CMPLX( 0.5e0 *( ap_p->rvector[ i1 ] ) *( ap_p->rvector[ i2 ] ) );
	  
      //      fprintf( stdout, "ap_p[ %d ] = %le, ap_p[ %d ] = %le.\n", i1, ap_p->rvector[ i1 ], i2, ap_p->rvector[ i2 ]);
	  
      if( MATRIX_SCALAR_PRODUCT( mu02_p->array[ index_mu02 ], dummy, mu02_p->array[ index_mu02 ] ) ) info=1; //WARNING: overwriting!

      /*	  
		 fprintf( stdout, "final\n" );

		 if( MATRIX_PRINT_PLUS( stdout,  mu02_p->array[ index_mu02 ] ) ) info=1;

		 fprintf( stdout, "mu02 norm = %le\n", MATRIX_NORM( mu02_p->array[ index_mu02 ] ) );
		 fprintf( stdout, "\n" );
      
		 fprintf( stdout, "-------------\n" );
      */

    } /* end i2 loop */

  } /* end i1 loop */

  
  return info;

}

//------------------------------------------

/* mu11_updated */

int mu11_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor_red;
  ivector         relevant_modes;
  int             max_rho_index;
  int             sqrt_max_rho_index;
  rvector_p       ar_p;
  rvector_p       ap_p;
  /* state */
  matrix_array_p  rho_p;
  matrix_array_p  mu11_p;
  matrix_p   dummy_matrix1_p;
  /* dummies */
  int             i1,i2,j;
  int             i_mode1, i_mode2;
  int             index, index_aux;
  int             index_mu11;
  int             info=0;
  complex         dummy;


  N_coor_red         =  constants.N_coor_red;
  relevant_modes     =  constants.relevant_modes;
  max_rho_index      =  constants.max_rho_index;
  sqrt_max_rho_index =  constants.sqrt_max_rho_index;

  rho_p              = &(config_p->electrons.rho);
  mu11_p             = &(config_p->electrons.mu11);
  dummy_matrix1_p    = &(state_p->dummy_matrix1);
  ar_p               = &(state_p->phonons.ar); 
  ap_p               = &(state_p->phonons.ap); 


  for( i_mode1=0; i_mode1<N_coor_red; i_mode1++ ){

    if( info ) break;
    
    i1 = relevant_modes.ivector[ i_mode1 ];
    
    for( i_mode2=0; i_mode2<N_coor_red; i_mode2++ ){

      if( info ) break;

      i2 = relevant_modes.ivector[ i_mode2 ];

      index_mu11 = COORDINATE_INDEX( i1, i2 );

      /* set mu11[ i ] to zero */
      if( MATRIX_ZERO( mu11_p->array[ index_mu11 ] ) ) info=1;

      for( j=0; j<sqrt_max_rho_index; j++ ){
	
	if( info ) break;
	
	/* set dummy_matrix1 to zero */
	if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;
	
	index = diagonal_rho_indices.ivector[ j ];
	
	if( RHO_INDEX_INVERSE( index, dummy_rho_index1, dummy_rho_index2 ) ) info=1;

#ifdef __DEBUG__
	
	if( IVECTOR_ARE_EQUAL( dummy_rho_index1, dummy_rho_index2 ) ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  info=1;
	  
	  break;
	  
	}
	
#endif /* __DEBUG__ */

	// 1st part;

	index_aux = RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index, i_mode1, i_mode2 );

	//	fprintf( stdout, "index_aux = %d\n", index_aux );
	
	if( index_aux < max_rho_index ){ // WARNING: Useless

	  if( i1 != i2 ){

	    dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode1 ] *dummy_rho_index1.ivector[ i_mode2 ] ) ); // Why not by table lookup?

	  }
	  else{

	    dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode1 ] *( dummy_rho_index1.ivector[ i_mode1 ] -1 ) ) ); // Why not by table lookup?
	  
	  }

	  /*
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index_aux ], dummy, *dummy_matrix1_p  ) ) info=1;
	  
	  if( MATRIX_SUM( mu11_p->array[ index_mu11 ], *dummy_matrix1_p, mu11_p->array[ index_mu11 ] ) ) info=1; //WARNING: overwriting!
      
	}

      } /* end j loop */
    
      //      fprintf( stdout, "index_mu11 = %d\n", index_mu11 );

      /*
	if( MATRIX_PRINT( stdout,  mu11_p->array[ index_mu11 ] ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /* the observable is hermitian */
      if( MATRIX_ADJOINT( mu11_p->array[ index_mu11 ], *dummy_matrix1_p  ) ) info=1;
      
      if( MATRIX_DIF( mu11_p->array[ index_mu11 ], *dummy_matrix1_p, mu11_p->array[ index_mu11 ] ) ) info=1; //WARNING: overwriting!

      /*
	if( MATRIX_PRINT( stdout,  mu11_p->array[ index_mu11 ] ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /* global factor */
      dummy = CMPLX( 0.5e0 *( ar_p->rvector[ i1 ] ) *( ap_p->rvector[ i2 ] ) ); // WARNING: a bit weird: where has the i1<-->i2 symmetry gone?

      //      fprintf( stdout, "ar_p[ %d ] = %le, ar_p[ %d ] = %le.\n", i1, ar_p->rvector[ i1 ], i2, ar_p->rvector[ i2 ]);

      dummy = CMPLX_PRODUCT( CMPLX_I, dummy ); //WARNING: overwriting!
      
      if( MATRIX_SCALAR_PRODUCT( mu11_p->array[ index_mu11 ], dummy, mu11_p->array[ index_mu11 ] ) ) info=1; //WARNING: overwriting!

      /*
	if( MATRIX_PRINT( stdout,  mu11_p->array[ index_mu11 ] ) ) info=1;
	fprintf( stdout, "\n" );
      */

    } /* end i2 loop */

  } /* end i1 loop */

  
  return info;

}

//------------------------------------------

/* mu20_updated */

int mu20_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor_red;
  ivector         relevant_modes;
  int             max_rho_index;
  int             sqrt_max_rho_index;
  rvector_p       ar_p;
  /* state */
  matrix_array_p  rho_p;
  matrix_array_p  mu20_p;
  matrix_p   dummy_matrix1_p;
  /* dummies */
  int             i1,i2,j;
  int             i_mode1, i_mode2;
  int             index, index_aux;
  int             index_mu20;
  int             info=0;
  complex         dummy;


  N_coor_red         =  constants.N_coor_red;
  relevant_modes     =  constants.relevant_modes;
  max_rho_index      =  constants.max_rho_index;
  sqrt_max_rho_index =  constants.sqrt_max_rho_index;

  rho_p              = &(config_p->electrons.rho);
  mu20_p             = &(config_p->electrons.mu20);
  dummy_matrix1_p    = &(state_p->dummy_matrix1);
  ar_p               = &(state_p->phonons.ar); 


  for( i_mode1=0; i_mode1<N_coor_red; i_mode1++ ){

    if( info ) break;

    i1 = relevant_modes.ivector[ i_mode1 ];

    for( i_mode2=0; i_mode2<N_coor_red; i_mode2++ ){

      if( info ) break;

      i2 = relevant_modes.ivector[ i_mode2 ];

      index_mu20 = COORDINATE_INDEX( i1, i2 );

      /* set mu20[ i ] to zero */
      if( MATRIX_ZERO( mu20_p->array[ index_mu20 ] ) ) info=1;

      for( j=0; j<sqrt_max_rho_index; j++ ){
	
	if( info ) break;
	
	/* set dummy_matrix1 to zero */
	if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;
	
	index = diagonal_rho_indices.ivector[ j ];
	
	if( RHO_INDEX_INVERSE( index, dummy_rho_index1, dummy_rho_index2 ) ) info=1;

#ifdef __DEBUG__
	
	if( IVECTOR_ARE_EQUAL( dummy_rho_index1, dummy_rho_index2 ) ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  info=1;
	  
	  break;
	  
	}
	
#endif /* __DEBUG__ */

	// 1st part;

	if( i1 != i2 ){

	  index_aux = RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( index, i_mode2, i_mode1 ); //WARNING: note index inversion

	  if( index_aux < max_rho_index ){
	    
	    //	fprintf( stdout, "index_aux = %d\n", index_aux );
	    
	    dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode1 ] ) *( dummy_rho_index1.ivector[ i_mode2 ] +1 ) ) ); // Why not by table lookup?

	  }
	  else{

	    dummy = CMPLX_ZERO;

	  }

	}
	else{
	  
	  index_aux = index;

	  //	fprintf( stdout, "index_aux = %d\n", index_aux );
	    
	  dummy = CMPLX( (double) dummy_rho_index1.ivector[ i_mode1 ] + 0.5e0 ); // Why not by table lookup?
	  
	}

	/*
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/
	
	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index_aux ], dummy, *dummy_matrix1_p  ) ) info=1;
	  
	if( MATRIX_SUM( mu20_p->array[ index_mu20 ], *dummy_matrix1_p, mu20_p->array[ index_mu20 ] ) ) info=1; //WARNING: overwriting!
      

	// 2nd part;

	index_aux = RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index, i_mode1, i_mode2 );

	//	fprintf( stdout, "index_aux = %d\n", index_aux );
	
	if( index_aux < max_rho_index ){

	  if( i1 != i2 ){

	    dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode1 ] *dummy_rho_index1.ivector[ i_mode2 ] ) ); // Why not by table lookup?

	  }
	  else{

	    dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode1 ] *( dummy_rho_index1.ivector[ i_mode1 ] -1 ) ) ); // Why not by table lookup?
	  
	  }

	  /*
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */
	  
	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index_aux ], dummy, *dummy_matrix1_p  ) ) info=1;
	  
	  if( MATRIX_SUM( mu20_p->array[ index_mu20 ], *dummy_matrix1_p, mu20_p->array[ index_mu20 ] ) ) info=1; //WARNING: overwriting!
      
	}

      } /* end j loop */
    
      //      fprintf( stdout, "index_mu20 = %d\n", index_mu20 );

      /*
	if( MATRIX_PRINT( stdout,  mu20_p->array[ index_mu20 ] ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /* the observable is hermitian */
      if( MATRIX_ADJOINT( mu20_p->array[ index_mu20 ], *dummy_matrix1_p  ) ) info=1;
      
      if( MATRIX_SUM( mu20_p->array[ index_mu20 ], *dummy_matrix1_p, mu20_p->array[ index_mu20 ] ) ) info=1; //WARNING: overwriting!

      /*
	if( MATRIX_PRINT( stdout,  mu20_p->array[ index_mu20 ] ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /* global factor */
      dummy = CMPLX( 0.5e0 *( ar_p->rvector[ i1 ] ) *( ar_p->rvector[ i2 ] ) );

      //      fprintf( stdout, "ar_p[ %d ] = %le, ar_p[ %d ] = %le.\n", i1, ar_p->rvector[ i1 ], i2, ar_p->rvector[ i2 ]);
      
      if( MATRIX_SCALAR_PRODUCT( mu20_p->array[ index_mu20 ], dummy, mu20_p->array[ index_mu20 ] ) ) info=1; //WARNING: overwriting!

      /*
	if( MATRIX_PRINT( stdout,  mu20_p->array[ index_mu20 ] ) ) info=1;
	fprintf( stdout, "\n" );
      */

    } /* end i2 loop */

  } /* end i1 loop */


  return info;

}

//------------------------------------------
 
/* utilities */

//------------------------------------------

/* irrelevant_mu_update */

int irrelevant_mu_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  int             N_coor_red;
  ivector         relevant_modes;
  rvector_p       ar_p;
  rvector_p       ap_p;
  /* state */
  matrix_p        mu00_p;
  matrix_array_p  mu01_p;
  matrix_array_p  mu10_p;
  matrix_array_p  mu02_p;
  matrix_array_p  mu11_p;
  matrix_array_p  mu20_p;
  /* dummies */
  int             i;
  int             i_mode;
  int             index;
  int             info=0;
  complex         dummy;


  N_coor             =  constants.N_coor;
  N_coor_red         =  constants.N_coor_red;
  relevant_modes     =  constants.relevant_modes;


  mu00_p             = &(config_p->electrons.mu00);
  mu01_p             = &(config_p->electrons.mu01);
  mu10_p             = &(config_p->electrons.mu10);
  mu02_p             = &(config_p->electrons.mu02);
  mu11_p             = &(config_p->electrons.mu11);
  mu20_p             = &(config_p->electrons.mu20);
  ar_p               = &(state_p->phonons.ar); 
  ap_p               = &(state_p->phonons.ap); 


  /* set mu to zero */
  if( MATRIX_ARRAY_ZERO( *mu01_p ) ) info=1;
  if( MATRIX_ARRAY_ZERO( *mu10_p ) ) info=1;
  if( MATRIX_ARRAY_ZERO( *mu02_p ) ) info=1;
  if( MATRIX_ARRAY_ZERO( *mu11_p ) ) info=1;
  if( MATRIX_ARRAY_ZERO( *mu20_p ) ) info=1;


  /* set to zero */
  i_mode=0;

  for( i=0; i<N_coor; i++ ){

    if( ( i_mode < N_coor_red ) && ( i == relevant_modes.ivector[ i_mode ] ) ){

      i_mode++;

    }
    else{

      index = COORDINATE_INDEX( i, i );

      //
      dummy = CMPLX( 0.5e0 *( ap_p->rvector[ i ] ) *( ap_p->rvector[ i ] ) );

      if( MATRIX_SCALAR_PRODUCT( *mu00_p, dummy, mu02_p->array[ index ] ) ) info=1;

      //
      dummy = CMPLX( 0.5e0 *( ar_p->rvector[ i ] ) *( ar_p->rvector[ i ] ) );

      if( MATRIX_SCALAR_PRODUCT( *mu00_p, dummy, mu20_p->array[ index ] ) ) info=1;

    }

  } /* end i loop */


  return info;

} 

//------------------------------------------
