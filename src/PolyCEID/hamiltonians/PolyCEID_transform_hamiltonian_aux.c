
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
#include "PolyCEID_transform_hamiltonian_aux.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* transform_hamiltonian_aux */

int transform_hamiltonian_aux( constants constants, state_p state_p, config_p config_p, matrix_array_p F_matrix_new_p, rvector_p momenta_new_p ){

  /* dummies */
  int  info=0;


  if( !info ){

    if( TRANSFORM_F_MATRIX( constants, *state_p, *config_p, *F_matrix_new_p ) ) info=1;

  }

  if( !info ){

    if( TRANSFORM_MOMENTA( constants, *state_p, *config_p, *momenta_new_p ) ) info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "F_matrix [transformed ]\n" );
  if( MATRIX_ARRAY_PRINT_PLUS( stdout, *F_matrix_new_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "momenta [transformed ]\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *momenta_new_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------
//------------------------------------------

/* transform_F_matrix */

int transform_F_matrix( constants constants, state_p state_p, config_p config_p, matrix_array_p F_matrix_new_p ){

  /* constants */
  int             N_coor;
  int             N_levels_many;
  /* state */
  matrix_array_p  F_matrix_p;
  vector_p        dummy_vector_coor1_p;
  vector_p        dummy_vector_coor2_p;
  matrix_p        conjugate_symplectic_transform_R_p;
  /* dummies */
  int             i, r;
  int             dim;
  int             info=0;

  
  N_coor                              =  constants.N_coor;
  N_levels_many                       =  constants.N_levels_many;

  F_matrix_p                          = &(config_p->electrons.F_matrix);
  dummy_vector_coor1_p                = &(state_p->dummy_vector_coor1);
  dummy_vector_coor2_p                = &(state_p->dummy_vector_coor2);
  conjugate_symplectic_transform_R_p  = &(state_p->phonons.conjugate_symplectic_transform_R);


  /* matrix_dim*/
  dim = N_levels_many *N_levels_many;


  for( i=0; i<dim; i++ ){

    /* create the vector */
    for( r=0; r<N_coor; r++ ){

      dummy_vector_coor1_p->vector[ r ] = F_matrix_p->array[ r ].matrix[i]; // it takes the entriy i of each matrix

    }


    if( MATRIX_VECTOR_PRODUCT( *conjugate_symplectic_transform_R_p, *dummy_vector_coor1_p, *dummy_vector_coor2_p ) ) info=1;


    /* copy the vector */
    for( r=0; r<N_coor; r++ ){

      F_matrix_new_p->array[ r ].matrix[i] = dummy_vector_coor2_p->vector[ r ]; // it copies the entriy i of each matrix

    }

  } /* end i loop */


  return info;

}

//------------------------------------------

/* transform_momenta */

int transform_momenta( constants constants, state_p state_p, config_p config_p, rvector_p momenta_new_p ){

  /* constants */
  int             N_coor;
  /* state */
  rvector_p       momenta_p;
  vector_p        dummy_vector_coor1_p;
  vector_p        dummy_vector_coor2_p;
  matrix_p        conjugate_symplectic_transform_R_p;
  /* dummies */
  int             i;
  int             info=0;


  N_coor                              =  constants.N_coor;

  momenta_p                           = &(config_p->atoms.momenta);
  dummy_vector_coor1_p                = &(state_p->dummy_vector_coor1);
  dummy_vector_coor2_p                = &(state_p->dummy_vector_coor2);
  conjugate_symplectic_transform_R_p  = &(state_p->phonons.conjugate_symplectic_transform_R);


  /* initial copy */
  for( i=0; i<N_coor; i++){

    dummy_vector_coor1_p->vector[ i ] = CMPLX( momenta_p->rvector[ i ] );

  }


  if( MATRIX_VECTOR_PRODUCT( *conjugate_symplectic_transform_R_p, *dummy_vector_coor1_p, *dummy_vector_coor2_p ) ) info=1;


  /* final copy */
  for( i=0; i<N_coor; i++){

    momenta_new_p->rvector[ i ] = REAL( dummy_vector_coor2_p->vector[ i ] );

  }


  return info;

}

//------------------------------------------

/* transform_forces */

int transform_forces( constants constants, state_p state_p, config_p config_p, rvector_p forces_new_p ){

  /* constants */
  int             N_coor;
  /* state */
  rvector_p       forces_p;
  vector_p        dummy_vector_coor1_p;
  vector_p        dummy_vector_coor2_p;
  matrix_p        conjugate_symplectic_transform_R_p;
  /* dummies */
  int             i;
  int             info=0;


  N_coor                              =  constants.N_coor;

  forces_p                            = &(config_p->atoms.forces);
  dummy_vector_coor1_p                = &(state_p->dummy_vector_coor1);
  dummy_vector_coor2_p                = &(state_p->dummy_vector_coor2);
  conjugate_symplectic_transform_R_p  = &(state_p->phonons.conjugate_symplectic_transform_R);


  /* initial copy */
  for( i=0; i<N_coor; i++){

    dummy_vector_coor1_p->vector[ i ] = CMPLX( forces_p->rvector[ i ] );

  }


  if( MATRIX_VECTOR_PRODUCT( *conjugate_symplectic_transform_R_p, *dummy_vector_coor1_p, *dummy_vector_coor2_p ) ) info=1;


  /* final copy */
  for( i=0; i<N_coor; i++){

    forces_new_p->rvector[ i ] = REAL( dummy_vector_coor2_p->vector[ i ] );

  }


  return info;

}

//------------------------------------------

