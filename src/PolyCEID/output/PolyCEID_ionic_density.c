
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
#include "PolyCEID_ionic_density.h"


/*********************
  FUNCTIONS & MACROS
*********************/


//BUGFIX: this was originally DEF
int ionic_density_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int              sqrt_max_rho_index;
  /* state */
  matrix_array_p   rho_p;
  matrix_p         ionic_density_p;
  rvector_p        ionic_density_eigenvalues_p;
  matrix_p         ionic_density_eigenvectors_p;
  /* dummies */
  int              i, j;
  int              index;
  int              info=0;


  sqrt_max_rho_index           =  constants.sqrt_max_rho_index;

  rho_p                        = &(config_p->electrons.rho);
  ionic_density_p              = &(state_p->ionic_density);
  ionic_density_eigenvalues_p  = &(state_p->ionic_density_eigenvalues);
  ionic_density_eigenvectors_p = &(state_p->ionic_density_eigenvectors);


  /* initialise ionic density */
  for( i=0; i<sqrt_max_rho_index; i++ ){
  
    for( j=0; j<sqrt_max_rho_index; j++ ){

      index=RHO_INDEX_AUX( i, j );
    
      ionic_density_p->matrix[ index ] = MATRIX_TRACE( rho_p->array[ index ] );

    } /* j loop */

  } /* i loop */
 


  /* ionic_density diagonalisation */
  if( DIAGONALISATION( *ionic_density_p, *ionic_density_eigenvectors_p, *ionic_density_eigenvalues_p ) )  info=1;


  return info;

}
