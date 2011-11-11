
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
#include "PolyCEID_transform_hamiltonian.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* transform_hamiltonian */

int transform_hamiltonian( constants constants, state_p state_p, config_p config_p, matrix_array_p delta_F_matrix_new_p, matrix_array_p K_matrix_new_p ){

  /* dummies */
  int  info=0;


  if( !info ){

    if( TRANSFORM_DELTA_F_MATRIX( constants, *state_p, *config_p, *delta_F_matrix_new_p ) ) info=1;

  }

  if( !info ){

    if( TRANSFORM_K_MATRIX( constants, *state_p, *config_p, *K_matrix_new_p ) ) info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "delta_F_matrix [transformed ]\n" );
  if( MATRIX_ARRAY_PRINT_PLUS( stdout, *delta_F_matrix_new_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "K_matrix [transformed ]\n" );
  if( MATRIX_ARRAY_PRINT_PLUS( stdout, *K_matrix_new_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------
//------------------------------------------

/* transform_delta_F_matrix */

int transform_delta_F_matrix( constants constants, state_p state_p, config_p config_p, matrix_array_p delta_F_matrix_new_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  matrix_array_p  delta_F_matrix_p;
  vector_p        dummy_vector_coor1_p;
  vector_p        dummy_vector_coor2_p;
  matrix_p        conjugate_symplectic_transform_R_p;
  /* dummies */
  int             i, r;
  int             dim;
  int             info=0;

  
  N_coor                              =  constants.N_coor;
  N_levels_many                       =  constants.N_levels_many;

  delta_F_matrix_p                    = &(config_p->electrons.delta_F_matrix);
  dummy_vector_coor1_p                = &(state_p->dummy_vector_coor1);
  dummy_vector_coor2_p                = &(state_p->dummy_vector_coor2);
  conjugate_symplectic_transform_R_p  = &(state_p->phonons.conjugate_symplectic_transform_R);


  /* matrix_dim*/
  dim = N_levels_many *N_levels_many;


  for( i=0; i<dim; i++ ){

    /* create the vector */
    for( r=0; r<N_coor; r++ ){

      dummy_vector_coor1_p->vector[ r ] = delta_F_matrix_p->array[ r ].matrix[i]; // it takes the entriy i of each matrix

    }


    if( MATRIX_VECTOR_PRODUCT( *conjugate_symplectic_transform_R_p, *dummy_vector_coor1_p, *dummy_vector_coor2_p ) ) info=1;


    /* copy the vector */
    for( r=0; r<N_coor; r++ ){

      delta_F_matrix_new_p->array[ r ].matrix[i] = dummy_vector_coor2_p->vector[ r ]; // it copies the entriy i of each matrix

    }

  } /* end i loop */


  return info;

}

//------------------------------------------

/* transform_K_matrix */

int transform_K_matrix( constants constants, state_p state_p, config_p config_p, matrix_array_p K_matrix_new_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  matrix_array_p  K_matrix_p;
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

  K_matrix_p                          = &(config_p->electrons.K_matrix);
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

	dummy_matrix_coor1_p->matrix[ index ] = K_matrix_p->array[ index ].matrix[ i ]; // it takes the entriy i of each matrix

      }

    }


    if( MATRIX_MATRIX_PRODUCT( *conjugate_symplectic_transform_R_p, *dummy_matrix_coor1_p, *dummy_matrix_coor2_p ) ) info=1;

    if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_coor2_p, *symplectic_transform_R_p, *dummy_matrix_coor1_p ) ) info=1;


    /* copy the vector */
    for( r=0; r<N_coor; r++ ){

      for( s=0; s<N_coor; s++ ){

	index = COORDINATE_INDEX( r,s );

	K_matrix_new_p->array[ index ].matrix[ i ] = dummy_matrix_coor1_p->matrix[ index ]; // it copies the entriy i of each matrix

      }

    }

  } /* end i loop */

  
  return info;

}

//------------------------------------------
