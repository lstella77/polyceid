
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
#include "PolyCEID_one_body_electronic_density_matrix.h"


/*********************
  FUNCTIONS & MACROS
*********************/

//BUGFIX: this was originally DEF
int compute_one_body_electronic_density_matrix( const constants constants, state_p state_p, config_p config_p ){

  /* constanst */

  /* state */
  matrix_p   one_body_electronic_density_matrix_p;
  matrix_p   one_body_electronic_density_matrix_Ehrenfest_p;
  matrix_p   Ehrenfest_frame_p;
  matrix_p   dummy_matrix_single1_p;
  matrix_p   dummy_matrix_single2_p;
  /* dummies */
  int        info=0;


  one_body_electronic_density_matrix_p           = &(state_p->one_body_electronic_density_matrix);
  one_body_electronic_density_matrix_Ehrenfest_p = &(state_p->one_body_electronic_density_matrix_Ehrenfest);
  Ehrenfest_frame_p                              = &(config_p->electrons.Ehrenfest_frame);
  dummy_matrix_single1_p                         = &(state_p->dummy_matrix_single1);
  dummy_matrix_single2_p                         = &(state_p->dummy_matrix_single2);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_one_body_electronic_density_matrix\n" );

#endif /* __DEBUG_PLUS__ */


  if( COMPUTE_ONE_BODY_ELECTRONIC_DENSITY_MATRIX_EHRENFEST( constants, *state_p, *config_p ) ) info=1;

  /*
    fprintf( stdout, "one_body_electronic_density_matrix_Ehrenfest\n" );
    if( MATRIX_PRINT_PLUS( stdout, *one_body_electronic_density_matrix_Ehrenfest_p ) ) info=1;
    fprintf( stdout, "\n" );
  */

  if( MATRIX_ADJOINT( *Ehrenfest_frame_p, *dummy_matrix_single1_p ) ) info=1;

  /*
    fprintf( stdout, "Ehrenfest_frame\n" );
    if( MATRIX_PRINT_PLUS( stdout, *Ehrenfest_frame_p ) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "dummy_matrix_single1\n" );
    if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
    fprintf( stdout, "\n" );
  */

  if( MATRIX_MATRIX_PRODUCT( *one_body_electronic_density_matrix_Ehrenfest_p, *dummy_matrix_single1_p, *dummy_matrix_single2_p ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( *Ehrenfest_frame_p, *dummy_matrix_single2_p, *one_body_electronic_density_matrix_p ) ) info=1;

  /*
    fprintf( stdout, "one_body_electronic_density_matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, *one_body_electronic_density_matrix_p ) ) info=1;
    fprintf( stdout, "\n" );
  */


  return info;

}


//------------------------------------------

/* utilities */

//------------------------------------------

//BUGFIX: this was originally DEF
int compute_one_body_electronic_density_matrix_Ehrenfest( const constants constants, state_p state_p, config_p config_p ){

  /* constanst */
  int        N_levels_single;
  int        N_levels_many;
  imatrix    single_level_occupation;
  imatrix    many_body_hybridisation_table;
  vector     symmetry_coefficient;
  /* state */
  matrix_p   mu00_p;
  matrix_p   one_body_electronic_density_matrix_Ehrenfest_p;
  /* dummies */
  complex    dummy;
  complex    symmetry_coefficient_loc;
  int        level1, level2;
  int        counter;
  int        i, j, k;
  int        index;
  int        info=0;


  N_levels_single                                 =  constants.N_levels_single;
  N_levels_many                                   =  constants.N_levels_many;
  single_level_occupation                         =  constants.single_level_occupation;
  many_body_hybridisation_table                   =  constants.many_body_hybridisation_table;
  symmetry_coefficient                            =  constants.symmetry_coefficient;

  one_body_electronic_density_matrix_Ehrenfest_p  = &(state_p->one_body_electronic_density_matrix_Ehrenfest);
  mu00_p                                          = &(config_p->electrons.mu00);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_one_body_electronic_density_matrix_Ehrenfest\n" );

#endif /* __DEBUG_PLUS__ */


  /* set matrix to zero */
  if( MATRIX_ZERO( *one_body_electronic_density_matrix_Ehrenfest_p ) ) info=1; // Important


  /* set counter to zero */
  counter=0;

  for( i=0; i<N_levels_many; i++ ){


    // diagonal terms
    for( k=0; k<N_levels_single*SPIN_DEG; k++ ){

      //      fprintf( stdout, "single_level %d, single_level_occupation %d\n", k, single_level_occupation.imatrix[i][ k ] );

      if( single_level_occupation.imatrix[i][ k ] ){

	level1 = k /SPIN_DEG;

	//	fprintf( stdout, "OK, level1 = %d\n", level1 );

	index = ELECTRON_SINGLE_INDEX( level1, level1 );

	/* writing */
	one_body_electronic_density_matrix_Ehrenfest_p->matrix[ index ] = CMPLX_SUM( one_body_electronic_density_matrix_Ehrenfest_p->matrix[ index ], mu00_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ] );

      } /* end conditional */

    } /* end k loop */


    // off-diagonal terms
    for( j=(i+1); j<N_levels_many; j++ ){


      if( many_body_hybridisation_table.imatrix[counter][0] != i ||
	  many_body_hybridisation_table.imatrix[counter][1] != j ){

	fprintf( stderr, "ERROR: inconsistency found in many_body_hybridisation_table ordering [%d,%d]\n", i, j );
	fflush( stderr );

	info=1;

      }
      else{

	level1 = many_body_hybridisation_table.imatrix[counter][2];
	level2 = many_body_hybridisation_table.imatrix[counter][4];
	symmetry_coefficient_loc = symmetry_coefficient.vector[counter];

	counter++;


	if( level1 != N_levels_single && level2 != N_levels_single ){

	  /* writing */
	  index = ELECTRON_SINGLE_INDEX( level2, level1 ); //WARNING: reverse ordering!

	  dummy = CMPLX_PRODUCT( symmetry_coefficient_loc, mu00_p->matrix[ ELECTRON_MANY_INDEX( j, i ) ] );

	  one_body_electronic_density_matrix_Ehrenfest_p->matrix[ index ] = CMPLX_SUM( one_body_electronic_density_matrix_Ehrenfest_p->matrix[ index ], dummy );

	  //
	  index = ELECTRON_SINGLE_INDEX( level1, level2 ); //WARNING: reverse ordering!

	  dummy = CMPLX_PRODUCT( CONJ( symmetry_coefficient_loc ), mu00_p->matrix[ ELECTRON_MANY_INDEX( i, j ) ] ); //WARNING: CONJ here is redondant

	  one_body_electronic_density_matrix_Ehrenfest_p->matrix[ index ] = CMPLX_SUM( one_body_electronic_density_matrix_Ehrenfest_p->matrix[ index ], dummy );

	  
	}

      } /* end conditional */

    } /* end j loop */

  } /* end i loop */


  return info;

}

//------------------------------------------
