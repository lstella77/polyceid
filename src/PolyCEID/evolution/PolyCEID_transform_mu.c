
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
#include "PolyCEID_transform_mu.h"



/*********************
  FUNCTIONS & MACROS
*********************/

/* transform_mu */

int transform_mu( constants constants, state_p state_p, config_p config_p){

  /* dummies */
  int  info=0;


  if( !info ){

    if( TRANSFORM_MU01( constants, *state_p, *config_p ) ) info=1;

  }

  if( !info ){

    if( TRANSFORM_MU10( constants, *state_p, *config_p ) ) info=1;

  }

  if( !info ){

    if( TRANSFORM_MU02( constants, *state_p, *config_p ) ) info=1;

  }

  if( !info ){
    
    if( TRANSFORM_MU11( constants, *state_p, *config_p ) ) info=1;
    
  }

  if( !info ){

    if( TRANSFORM_MU20( constants, *state_p, *config_p ) ) info=1;

  }


  return info;

}

//------------------------------------------
//------------------------------------------

/* transform_mu01f */

int transform_mu01( constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  matrix_array_p  mu01_p;
  vector_p        dummy_vector_coor1_p;
  vector_p        dummy_vector_coor2_p;
  matrix_p        symplectic_transform_P_p;
  /* dummies */
  int             i, r;
  int             dim;
  int             info=0;


  N_coor                    =  constants.N_coor;
  N_levels_many             =  constants.N_levels_many;

  mu01_p                    = &(config_p->electrons.mu01);
  dummy_vector_coor1_p      = &(state_p->dummy_vector_coor1);
  dummy_vector_coor2_p      = &(state_p->dummy_vector_coor2);
  symplectic_transform_P_p  = &(state_p->phonons.symplectic_transform_P);


  /* matrix_dim*/
  dim = N_levels_many *N_levels_many;


  for( i=0; i<dim; i++ ){

    /* create the vector */
    for( r=0; r<N_coor; r++ ){

      dummy_vector_coor1_p->vector[ r ] = mu01_p->array[ r ].matrix[i]; // it takes the entriy i of each matrix

    }


    if( MATRIX_VECTOR_PRODUCT( *symplectic_transform_P_p, *dummy_vector_coor1_p, *dummy_vector_coor2_p ) ) info=1;


    /* copy the vector */
    for( r=0; r<N_coor; r++ ){

      mu01_p->array[ r ].matrix[i] = dummy_vector_coor2_p->vector[ r ]; // it copies the entriy i of each matrix

    }

  } /* end i loop */


  return info;

}

//------------------------------------------

/* transform_mu10 */

int transform_mu10( constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  matrix_array_p  mu10_p;
  vector_p        dummy_vector_coor1_p;
  vector_p        dummy_vector_coor2_p;
  matrix_p        symplectic_transform_R_p;
  /* dummies */
  int             i, r;
  int             dim;
  int             info=0;


  N_coor                    =  constants.N_coor;
  N_levels_many             =  constants.N_levels_many;

  mu10_p                    = &(config_p->electrons.mu10);
  dummy_vector_coor1_p      = &(state_p->dummy_vector_coor1);
  dummy_vector_coor2_p      = &(state_p->dummy_vector_coor2);
  symplectic_transform_R_p  = &(state_p->phonons.symplectic_transform_R);


  /* matrix_dim*/
  dim = N_levels_many *N_levels_many;


  for( i=0; i<dim; i++ ){

    /* create the vector */
    for( r=0; r<N_coor; r++ ){

      dummy_vector_coor1_p->vector[ r ] = mu10_p->array[ r ].matrix[i]; // it takes the entriy i of each matrix

    }


    if( MATRIX_VECTOR_PRODUCT( *symplectic_transform_R_p, *dummy_vector_coor1_p, *dummy_vector_coor2_p ) ) info=1;


    /* copy the vector */
    for( r=0; r<N_coor; r++ ){

      mu10_p->array[ r ].matrix[i] = dummy_vector_coor2_p->vector[ r ]; // it copies the entriy i of each matrix

    }

  } /* end i loop */


  return info;

}

//------------------------------------------

/* transform_mu02 */

int transform_mu02( constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  matrix_array_p  mu02_p;
  matrix_p        dummy_matrix_coor1_p;
  matrix_p        dummy_matrix_coor2_p;
  matrix_p        symplectic_transform_P_p;
  matrix_p        conjugate_symplectic_transform_P_p;
  /* dummies */
  int             i, r, s;
  int             index;
  int             dim;
  int             info=0;


  N_coor                              =  constants.N_coor;
  N_levels_many                       =  constants.N_levels_many;

  mu02_p                              = &(config_p->electrons.mu02);
  dummy_matrix_coor1_p                = &(state_p->dummy_matrix_coor1);
  dummy_matrix_coor2_p                = &(state_p->dummy_matrix_coor2);
  symplectic_transform_P_p            = &(state_p->phonons.symplectic_transform_P);
  conjugate_symplectic_transform_P_p  = &(state_p->phonons.conjugate_symplectic_transform_P);

  /* matrix_dim*/
  dim = N_levels_many *N_levels_many;


  for( i=0; i<dim; i++ ){

    /* create the vector */
    for( r=0; r<N_coor; r++ ){

      for( s=0; s<N_coor; s++ ){

	index = COORDINATE_INDEX( r,s );

	dummy_matrix_coor1_p->matrix[ index ] = mu02_p->array[ index ].matrix[ i ]; // it takes the entriy i of each matrix

      }

    }


    if( MATRIX_MATRIX_PRODUCT( *symplectic_transform_P_p, *dummy_matrix_coor1_p, *dummy_matrix_coor2_p ) ) info=1;

    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_coor2_p, *conjugate_symplectic_transform_P_p, *dummy_matrix_coor1_p ) ) info=1;


    /* copy the vector */
    for( r=0; r<N_coor; r++ ){

      for( s=0; s<N_coor; s++ ){

	index = COORDINATE_INDEX( r,s );

	mu02_p->array[ index ].matrix[ i ] = dummy_matrix_coor1_p->matrix[ index ]; // it copies the entriy i of each matrix

      }

    }

  } /* end i loop */


  return info;

}

//------------------------------------------

/* transform_mu11 */

int transform_mu11( constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  matrix_p        mu00_p;
  matrix_array_p  mu11_p;
  matrix_p        dummy_matrix_coor1_p;
  matrix_p        dummy_matrix_coor2_p;
  matrix_p        symplectic_transform_R_p;
  matrix_p        conjugate_symplectic_transform_P_p;
  /* dummies */
  int             i, r, s;
  int             index;
  int             dim;
  int             info=0;
  complex         dummy;


  N_coor                              =  constants.N_coor;
  N_levels_many                       =  constants.N_levels_many;

  mu00_p                              = &(config_p->electrons.mu00);
  mu11_p                              = &(config_p->electrons.mu11);
  dummy_matrix_coor1_p                = &(state_p->dummy_matrix_coor1);
  dummy_matrix_coor2_p                = &(state_p->dummy_matrix_coor2);
  symplectic_transform_R_p            = &(state_p->phonons.symplectic_transform_R);
  conjugate_symplectic_transform_P_p  = &(state_p->phonons.conjugate_symplectic_transform_P);


  /* matrix_dim*/
  dim = N_levels_many *N_levels_many;


  for( i=0; i<dim; i++ ){

    /* create the vector */
    for( r=0; r<N_coor; r++ ){

      for( s=0; s<N_coor; s++ ){

	index = COORDINATE_INDEX( r, s );

	if( r != s ){

	  dummy_matrix_coor1_p->matrix[ index ] = mu11_p->array[ index ].matrix[ i ]; // it takes the entriy i of each matrix

	}
	else{

	  dummy_matrix_coor1_p->matrix[ index ] = mu11_p->array[ index ].matrix[ i ]; // it takes the entriy i of each matrix

	  dummy = CMPLX_PRODUCT( IHBAR, mu00_p->matrix[ i ] );

	  dummy = CMPLX_PRODUCT( CMPLX( 0.5e0 ), dummy ); //WARNING: overwriting

	  dummy_matrix_coor1_p->matrix[ index ] = CMPLX_DIF( dummy_matrix_coor1_p->matrix[ index ], dummy );  //WARNING: overwriting

	}

      }

    }


    if( MATRIX_MATRIX_PRODUCT( *symplectic_transform_R_p, *dummy_matrix_coor1_p, *dummy_matrix_coor2_p ) ) info=1;

    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_coor2_p, *conjugate_symplectic_transform_P_p, *dummy_matrix_coor1_p ) ) info=1;


    /* copy the vector */
    for( r=0; r<N_coor; r++ ){

      for( s=0; s<N_coor; s++ ){

	index = COORDINATE_INDEX( r, s );

	if( r != s ){

	  mu11_p->array[ index ].matrix[ i ] = dummy_matrix_coor1_p->matrix[ index ]; // it copies the entriy i of each matrix

	}
	else{

	  mu11_p->array[ index ].matrix[ i ] = dummy_matrix_coor1_p->matrix[ index ]; // it copies the entriy i of each matrix

	  dummy = CMPLX_PRODUCT( IHBAR, mu00_p->matrix[ i ] );

	  dummy = CMPLX_PRODUCT( CMPLX( 0.5e0 ), dummy ); //WARNING: overwriting

	  mu11_p->array[ index ].matrix[ i ] = CMPLX_SUM( mu11_p->array[ index ].matrix[ i ], dummy ); //WARNING: overwriting

	}

      }

    }

  } /* end i loop */


  return info;

}

//------------------------------------------

/* transform_mu20 */

int transform_mu20( constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  matrix_array_p  mu20_p;
  matrix_p        dummy_matrix_coor1_p;
  matrix_p        dummy_matrix_coor2_p;
  matrix_p        symplectic_transform_R_p;
  matrix_p        conjugate_symplectic_transform_R_p;
  /* dummies */
  int             i, r, s;
  int             index;
  int             dim;
  int             info=0;


  N_coor                              =  constants.N_coor;
  N_levels_many                       =  constants.N_levels_many;

  mu20_p                              = &(config_p->electrons.mu20);
  dummy_matrix_coor1_p                = &(state_p->dummy_matrix_coor1);
  dummy_matrix_coor2_p                = &(state_p->dummy_matrix_coor2);
  symplectic_transform_R_p            = &(state_p->phonons.symplectic_transform_R);
  conjugate_symplectic_transform_R_p  = &(state_p->phonons.conjugate_symplectic_transform_R);

  /* matrix_dim*/
  dim = N_levels_many *N_levels_many;


  for( i=0; i<dim; i++ ){

    /* create the vector */
    for( r=0; r<N_coor; r++ ){

      for( s=0; s<N_coor; s++ ){

	index = COORDINATE_INDEX( r,s );

	dummy_matrix_coor1_p->matrix[ index ] = mu20_p->array[ index ].matrix[ i ]; // it takes the entriy i of each matrix

      }

    }


    if( MATRIX_MATRIX_PRODUCT( *symplectic_transform_R_p, *dummy_matrix_coor1_p, *dummy_matrix_coor2_p ) ) info=1;

    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_coor2_p, *conjugate_symplectic_transform_R_p, *dummy_matrix_coor1_p ) ) info=1;


    /* copy the vector */
    for( r=0; r<N_coor; r++ ){

      for( s=0; s<N_coor; s++ ){

	index = COORDINATE_INDEX( r,s );

	mu20_p->array[ index ].matrix[ i ] = dummy_matrix_coor1_p->matrix[ index ]; // it copies the entriy i of each matrix

      }

    }

  } /* end i loop */


  return info;

}

//------------------------------------------
