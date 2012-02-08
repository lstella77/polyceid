
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
#include "PolyCEID_distances.h"


/*********************
  FUNCTIONS & MACROS
*********************/


/* distances update */
/* WARNING: BLAS should be used instead */

int distances_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        N_atoms;
  int        sdim;
  /* state */
  rvector_p  positions_p;
  rvector_p  distances_p;
  /* dummies */
  int       index;
  int       i, j, comp;
  int       info=0;
  double    dummy;


  N_atoms         =  constants.N_atoms;
  sdim            =  constants.spacial_dimension; 
  positions_p     = &config_p->atoms.positions;
  distances_p     = &config_p->atoms.distances;


  if( !info ){

    for( i=0; i<N_atoms; i++ ){

      for( j=0; j<N_atoms; j++ ){

	//index = MATRIX_ARG( i, j, N_atoms, N_atoms );
	index = ATOM_INDEX( i, j );

	if( j != i ){
	  
	  distances_p->rvector[ index ] = 0.0e0;

	  for( comp=0; comp<sdim; comp++ ){

	    dummy = positions_p->rvector[ comp +i *sdim ] -positions_p->rvector[ comp +j *sdim ];

	    if( constants.flag_periodic_boundary_condition ){
                
              dummy = PBC_DIFF( dummy, constants.cell_dim.rvector[ comp ] );

	    }

	    dummy = dummy *dummy;

	    distances_p->rvector[ index ] += dummy;

	  }

	  distances_p->rvector[ index ] = sqrt( distances_p->rvector[ index ] ); // WARNING: is it possible to avoid it?

	  if( distances_p->rvector[ index ] < EPS ){

	    fprintf( stderr, "ERROR: atoms %d and %d are superimposed\n", i, j );
	    fprintf( stderr, "positions: %f  %f\n", positions_p->rvector[ i ],  positions_p->rvector[ j ] );
	    fprintf( stderr, "distance: %f\n", distances_p->rvector[ index ] );
	    fflush( stderr );

	    info=1;

	  }

	}
	else{

	  distances_p->rvector[ index ] = 0.0e0;

	}

	if( info ) break;

      }

      if( info ) break;

    }

  }


  return info;

}

//------------------------------------------
