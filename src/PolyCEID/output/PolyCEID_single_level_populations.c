
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
#include "PolyCEID_single_level_populations.h"


/*********************
  FUNCTIONS & MACROS
*********************/

//BUGFIX; this was originally DEF
int compute_single_level_populations( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        N_levels_single;
  /* state */
  matrix     one_body_electronic_density_matrix;
  matrix     one_body_electronic_density_matrix_Ehrenfest;
  rvector_p  adiabatic_PES_single_p;
  matrix_p   adiabatic_states_single_p;
  rvector_p  single_level_populations_adiabatic_p;
  rvector_p  single_level_populations_Ehrenfest_p;
  matrix_p   dummy_matrix_single1_p;
  matrix_p   dummy_matrix_single2_p;
  /* dummies */
  int     i;
  int     index;
  int     info=0;


  N_levels_single                               =   constants.N_levels_single;

  one_body_electronic_density_matrix            =   state_p->one_body_electronic_density_matrix;
  one_body_electronic_density_matrix_Ehrenfest  =   state_p->one_body_electronic_density_matrix_Ehrenfest;
  adiabatic_PES_single_p                        = &(state_p->adiabatic_PES_single);
  adiabatic_states_single_p                     = &(state_p->adiabatic_states_single);
  single_level_populations_adiabatic_p          = &(state_p->single_level_populations_adiabatic);
  single_level_populations_Ehrenfest_p          = &(state_p->single_level_populations_Ehrenfest);

  dummy_matrix_single1_p                        = &(state_p->dummy_matrix_single1);
  dummy_matrix_single2_p                        = &(state_p->dummy_matrix_single2);


 if( H_MATRIX_SINGLE_UPDATE( constants, *state_p, *config_p, *dummy_matrix_single1_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "H_matrix_single\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix_single1_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    if( DIAGONALISATION( *dummy_matrix_single1_p, *adiabatic_states_single_p, *adiabatic_PES_single_p ) ) info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "adiabatic_states_single\n" );
  if( MATRIX_PRINT_PLUS( stdout, *adiabatic_states_single_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "adiabatic_PES_single\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *adiabatic_PES_single_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  // single_level_populations_Ehrenfest
  for( i=0; i<N_levels_single; i++ ){

    index = ELECTRON_SINGLE_INDEX( i, i );

    single_level_populations_Ehrenfest_p->rvector[ i ] = REAL( one_body_electronic_density_matrix_Ehrenfest.matrix[ index ] );

  } /* end i_single_E loop */


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "single_level_populations_Ehrenfest\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *single_level_populations_Ehrenfest_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  // chenge of frame
  if( MATRIX_ADJOINT( *adiabatic_states_single_p, *dummy_matrix_single1_p ) ) info=1;


  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single1_p, one_body_electronic_density_matrix, *dummy_matrix_single2_p ) ) info=1;

  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix_single2_p, *adiabatic_states_single_p, *dummy_matrix_single1_p ) )            info=1;


  // single_level_populations_adiabatic
  for( i=0; i<N_levels_single; i++ ){

    index = ELECTRON_SINGLE_INDEX( i, i );

    single_level_populations_adiabatic_p->rvector[ i ] = REAL( dummy_matrix_single1_p->matrix[ index ] );

  } /* end i_single_A loop */


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "single_level_populations_adiabatic\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *single_level_populations_Ehrenfest_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}
