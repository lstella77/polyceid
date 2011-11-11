
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
#include "PolyCEID_hamiltonian_many_aux.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* hamiltonian_many_aux update */

int hamiltonian_many_aux_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  /* states */
  /* dummies */
  int                 comp;
  int                 info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: hamiltonian_many_aux_update\n" );

#endif /* __DEBUG_PLUS__*/


#if defined( __DEBUG__ ) && defined( __HAMILTONIAN_CONSISTENCY__ )

  /* check consistency */
  if( HAMILTONIAN_SINGLE_CHECK( constants, *state_p, *config_p ) ) info=1;

#endif /* defined( __DEBUG__ ) && defined( __HAMILTONIAN_CONSISTENCY__ ) */


  /* distances update */
  // WARNING: and it is needed by the next
  if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;


  /* centre of mass update */
  // WARNING: and it is needed by the next
  if( COMPUTE_CENTRE_OF_MASS( constants, config_p->atoms ) ) info=1;


  /* H_matrix_many update */
  if( H_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;
      
    
  /* F_matrix_many update */
  if( F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p  ) ) info=1;
  

#if defined( __DEBUG__ ) && defined( __HAMILTONIAN_CONSISTENCY__ )
    
  /* check consistency */
  if( HAMILTONIAN_MANY_CHECK( constants, *state_p, *config_p ) ) info=1;

#endif /* defined( __DEBUG__ ) && defined( __HAMILTONIAN_CONSISTENCY__ ) */


  // WARNING: here there should be the loop over cartesian components...
  comp = 0;

  /* K_matrix_many update */
  if( CLASSICAL_DIPOLE_MANY_UPDATE( constants, *state_p, *config_p, comp ) ) info=1;


  /* forces update */
  // WARNING: it must be here because it depends on the previous matrices
  // and it is needed by the next
  if( FORCES_UPDATE( constants, *state_p, *config_p ) ) info=1;


  /* delta_F_matrix_many update */
  if( DELTA_F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;
  

  return info;

}

//------------------------------------------
