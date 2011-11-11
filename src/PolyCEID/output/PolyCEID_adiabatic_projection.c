
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
#include "PolyCEID_adiabatic_projection.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int adiabatic_projection_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        N_levels_many;
  /* state */
  matrix_p   H_matrix_p;
  rvector_p  adiabatic_populations_p;
  matrix_p   adiabatic_projection_p;
  matrix_p   adiabatic_density_matrix_p;
  matrix_p   adiabatic_states_many_p;
  rvector_p  adiabatic_PES_many_p;
  rvector_p  electronic_density_eigenvalues_p;
  matrix_p   electronic_density_eigenvectors_p;
  matrix_p   mu00_p;
  matrix_p   dummy_matrix1_p;
  matrix_p   dummy_matrix2_p;
  /* dummies */
  int        i;
  int        info=0;


  N_levels_many                      =  constants.N_levels_many;
                            
  H_matrix_p                         = &(config_p->electrons.H_matrix);
  adiabatic_projection_p             = &(state_p->adiabatic_projection);
  adiabatic_populations_p            = &(state_p->adiabatic_populations);
  adiabatic_density_matrix_p         = &(state_p->adiabatic_density_matrix);
  mu00_p                             = &(config_p->electrons.mu00);
  electronic_density_eigenvalues_p   = &(state_p->electronic_density_eigenvalues);
  electronic_density_eigenvectors_p  = &(state_p->electronic_density_eigenvectors);
  adiabatic_states_many_p            = &(state_p->adiabatic_states_many);
  adiabatic_PES_many_p               = &(state_p->adiabatic_PES_many);
  dummy_matrix1_p                    = &(state_p->dummy_matrix1);
  dummy_matrix2_p                    = &(state_p->dummy_matrix2);


  // adiabatic frame
  if( DIAGONALISATION( *H_matrix_p, *adiabatic_states_many_p, *adiabatic_PES_many_p ) ) info=1;

  /* 
    fprintf( stdout, "adiabatic_PES_many\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *adiabatic_PES_many_p) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "adiabatic_states_many\n" );
    if( MATRIX_PRINT_PLUS( stdout, *adiabatic_states_many_p) ) info=1;
    fprintf( stdout, "\n" );
  */


  /* mu00 diagonalisation */
  if( DIAGONALISATION( *mu00_p, *electronic_density_eigenvectors_p, *electronic_density_eigenvalues_p ) )  info=1;

  /* 
    fprintf( stdout, "electronic_density_eigenvalues\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *electronic_density_eigenvalues_p) ) info=1;
    fprintf( stdout, "\n" );

    fprintf( stdout, "electronic_density_eigenvectors\n" );
    if( MATRIX_PRINT_PLUS( stdout, *electronic_density_eigenvectors_p) ) info=1;
    fprintf( stdout, "\n" );
  */


  /* unitary transform */
  if( MATRIX_ADJOINT( *adiabatic_states_many_p, *dummy_matrix1_p ) )                                            info=1;

  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, *electronic_density_eigenvectors_p, *adiabatic_projection_p ) )  info=1;

  /*
    fprintf( stdout, "electronic_density_matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, *electronic_density_eigenvectors_p) ) info=1;
    fprintf( stdout, "\n" );
  */


  /* adiabatic_density_matrix */
  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, *mu00_p, *dummy_matrix2_p ) )  info=1;

  if( MATRIX_MATRIX_PRODUCT( *dummy_matrix2_p, *adiabatic_states_many_p, *adiabatic_density_matrix_p ) )  info=1;

  /*
    fprintf( stdout, "adiabatic_density_matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, *adiabatic_density_matrix_p) ) info=1;
    fprintf( stdout, "\n" );
  */

  // compute adiabatic populations
  for( i=0; i<N_levels_many; i++ ){

    adiabatic_populations_p->rvector[ i ] = REAL( adiabatic_density_matrix_p->matrix[ ELECTRON_MANY_INDEX(i, i) ] );

  }


  return info;

}

//------------------------------------------
