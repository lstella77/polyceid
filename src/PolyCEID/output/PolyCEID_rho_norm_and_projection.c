
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
#include "PolyCEID_rho_norm_and_projection.h"



/*********************
  FUNCTIONS & MACROS
*********************/

/* rho_norm_and_projection_update */

//BUGFIX: this was originally DEF
int rho_norm_and_projection_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             max_rho_index;
  /* state */
  double*         rho_norm_p; 
  double*         rho_dot_norm_p; 
  double*         rho_projection_p; 
  matrix_array_p  rho_dot_p;
  matrix_array_p  rho_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  int             index;
  int             index_aux;
  int             info=0;
  complex         dummy;
  complex         projection;


#ifdef __DEBUG__

  if( !info ){
    
    if( RHO_HERMITICITY_TEST( constants, *state_p, *config_p ) ) info=1;
    
    if( RHO_DOT_HERMITICITY_TEST( constants, *state_p, *config_p ) ) info=1;
    
  }

#endif /* __DEBUG__ */


  max_rho_index         =  constants.max_rho_index;

  rho_norm_p            = &(state_p->rho_norm);
  rho_dot_norm_p        = &(state_p->rho_dot_norm);
  rho_projection_p      = &(state_p->rho_projection);
  rho_p                 = &(config_p->electrons.rho);
  rho_dot_p             = &(state_p->rho_dot);
  dummy_matrix1_p       = &(state_p->dummy_matrix1);


  /* set to zero */
  *rho_norm_p           = 0.0e0;
  *rho_dot_norm_p       = 0.0e0;
  projection            = CMPLX_ZERO;
  //  *rho_projection_p = 0.0e0; 


  for( index=0; index<max_rho_index; index++ ){

    index_aux = rho_index_conjugate.ivector[ index ];
    
    
    /* rho_norm */
    if( MATRIX_MATRIX_PRODUCT( rho_p->array[ index_aux ], rho_p->array[ index ], *dummy_matrix1_p ) ) info=1;
    
    dummy = MATRIX_TRACE( *dummy_matrix1_p );

    /*
    if( fabs( dummy.z[1] ) > EPS ){

      fprintf( stderr, "ERROR: unexpected imaginary part found at time %le, %d, %le].\n", config_p->time, index, dummy.z[1] );
      fflush( stderr );
      
      info=1;

      break;

    }
    */ 
   
    *rho_norm_p += REAL( dummy );

      
    /* rho_dot_norm */
    if( MATRIX_MATRIX_PRODUCT( rho_dot_p->array[ index_aux ], rho_dot_p->array[ index ], *dummy_matrix1_p ) ) info=1;
    
    dummy = MATRIX_TRACE( *dummy_matrix1_p );

    /*
    if( fabs( dummy.z[1] ) > EPS ){

      fprintf( stderr, "ERROR: unexpected imaginary part found at time %le, %d, %le].\n", config_p->time, index, dummy.z[1] );
      fflush( stderr );
      
      info=1;

      break;

    }
    */

    *rho_dot_norm_p += REAL( dummy );


    /* rho_projection */
    if( MATRIX_MATRIX_PRODUCT( rho_dot_p->array[ index_aux ], rho_p->array[ index ], *dummy_matrix1_p ) ) info=1;

    dummy = MATRIX_TRACE( *dummy_matrix1_p );

    projection = CMPLX_SUM( projection, dummy ); //WARNING: overwriting


  } /* end index loop */


  /* sqrt rho_norm */
  if( *rho_norm_p <0.0e0 ){

    fprintf( stderr, "ERROR: unexpected sqrt of a negative value found at time %le [%le].\n", config_p->time, *rho_norm_p );
    fflush( stderr );

    info=1;

  }
  else{

    *rho_norm_p = sqrt( *rho_norm_p );

  }


  /* sqrt rho_dot_norm */
  if( *rho_dot_norm_p <0.0e0 ){

    fprintf( stderr, "ERROR: unexpected sqrt of a negative value found at time %le [%le].\n", config_p->time, *rho_dot_norm_p );
    fflush( stderr );

    info=1;

  }
  else{

    *rho_dot_norm_p = sqrt( *rho_dot_norm_p );

  }


  /* checking projection*/
  if( fabs( projection.z[1] ) > EPS *( *rho_norm_p ) *( *rho_dot_norm_p ) && fabs( projection.z[1] ) > EPS ){
    
    fprintf( stderr, "ERROR: unexpected imaginary part found at time %le, %le].\n", config_p->time, projection.z[1] );
    fflush( stderr );
    
    info=1;
    
  }
  else{

    *rho_projection_p = REAL( projection );

  }


  if( fabs( *rho_projection_p ) > EPS *( *rho_norm_p ) *( *rho_dot_norm_p ) && fabs( *rho_projection_p ) > EPS ){

    fprintf( stderr, "ERROR: the projection of rho_dot on rho is not zero at time %le, %le]\n", config_p->time, *rho_projection_p );
    fflush( stderr ); 

    info=1;

  }


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

#ifdef __DEBUG__

/* rho_hermiticity_test */

//BUGFIX: this was originally DEF
int rho_hermiticity_test( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             max_rho_index;
  /* state */
  matrix_array_p  rho_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  int             index;
  int             index_aux;
  int             info=0;
  double          dummy;
  

  max_rho_index   =  constants.max_rho_index;

  rho_p           = &(config_p->electrons.rho);
  dummy_matrix1_p = &(state_p->dummy_matrix1);


  for( index=0; index<max_rho_index; index++ ){

    if( !info ){

      if( MATRIX_ADJOINT( rho_p->array[ index ], *dummy_matrix1_p ) ) info=1;

      index_aux = rho_index_conjugate.ivector[ index ];

      if( MATRIX_DIF( *dummy_matrix1_p, rho_p->array[ index_aux ], *dummy_matrix1_p ) ) info=1; // WARNING: overwriting

      dummy = MATRIX_NORM( *dummy_matrix1_p );

      if( dummy > EPS ){

	fprintf( stderr, "ERROR: Hermiticity violation in rho at time %le.\n", config_p->time );
	fflush( stderr );

	info=1;

	break;

      }

    }

  } /* end index loop */


  return info;

}

//------------------------------------------

/* rho_dot_hermiticity_test */

int rho_dot_hermiticity_test( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             max_rho_index;
  /* state */
  matrix_array_p  rho_dot_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  int             index;
  int             index_aux;
  int             info=0;
  double          dummy;

  
  max_rho_index   =  constants.max_rho_index;

  rho_dot_p       = &(state_p->rho_dot);
  dummy_matrix1_p = &(state_p->dummy_matrix1);


  for( index=0; index<max_rho_index; index++ ){

    if( !info ){

      if( MATRIX_ADJOINT( rho_dot_p->array[ index ], *dummy_matrix1_p ) ) info=1;

      index_aux = rho_index_conjugate.ivector[ index ];

      if( MATRIX_DIF( *dummy_matrix1_p, rho_dot_p->array[ index_aux ], *dummy_matrix1_p ) ) info=1; // WARNING: overwriting

      dummy = MATRIX_NORM( *dummy_matrix1_p );

      if( dummy > EPS ){

	fprintf( stderr, "ERROR: Hermiticity violation in rho_dot at time %le.\n", config_p->time );
	fflush( stderr );

	info=1;

	break;

      }

    }

  } /* end index loop */


  return info;

}

//------------------------------------------

#endif /* __DEBUG__*/
