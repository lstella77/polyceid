
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
#include "../utils/PolyCEID_sorting.h"


#define SWAP( a, b, c ) (c)=(a); (a)=(b); (b)=(c) //WARNING: the last variable is a dummy one


/* __my_qsort */

int __my_rvector_qsort( const rvector rvector, ivector_p ivector_p ){

  /* dummies */
  int i;
  int dim;
  int info=0;


  dim = rvector.rvector_dim;


  if( dim != ivector_p->ivector_dim ){

    fprintf( stderr, "ERROR: value and index vectors have different dimensions.\n");

    info=1;

  }


  if( !info ){

    /* initialisation index vector */
    for( i=0; i<dim; i++ ){ 

      ivector_p->ivector[ i ] = i;

    }


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "# original vector\n");
      
      for( i=0; i<dim; i++ ){

	fprintf( stdout, "# %d %le\n", ivector_p->ivector[ i ], rvector.rvector[ ivector_p->ivector[ i ] ] );

	
      }

    fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */


    if( RVECTOR_QSORT_AUX( rvector, *ivector_p, 0, dim-1 ) ) info=1;

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "# sorted vector\n");

  for( i=0; i<dim; i++ ){

    fprintf( stdout, "# %d %le\n", ivector_p->ivector[ i ], rvector.rvector[ ivector_p->ivector[ i ] ] );


  }

  fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int __my_rvector_qsort_aux( const rvector rvector, ivector_p ivector_p, int ini, int fin ){

  /* dummies */
  int i, last;
  int dummy;
  int info=0;


  if( rvector.rvector_dim != ivector_p->ivector_dim ){

    fprintf( stderr, "ERROR: value and index vectors have different dimensions.\n");

    info=1;

  }


  if( !info && ini < fin ){


    /* move partitioning elevent to the first position */
    SWAP( ivector_p->ivector[ ini ], ivector_p->ivector[ (ini+fin)/2 ], dummy );

    last = ini;

    for( i=ini+1; i<=fin; i++ ){ //WARNING: it skips the first element because it's the partitioning one

      if( rvector.rvector[ ivector_p->ivector[ i ] ] < rvector.rvector[ ivector_p->ivector[ ini ] ] ){

	last++;

	SWAP( ivector_p->ivector[ last ], ivector_p->ivector[ i ], dummy );

      } 

    }

    /* restore partioning element in the right position */
    SWAP( ivector_p->ivector[ ini ], ivector_p->ivector[ last ], dummy );

    /* recurrency */
    if( RVECTOR_QSORT_AUX( rvector, *ivector_p, ini, last-1 ) ) info=1;

    if( RVECTOR_QSORT_AUX( rvector, *ivector_p, last+1, fin ) ) info=1;

  }


  return info;

}

//------------------------------------------


