
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
#include "PolyCEID_combinatorics.h"

#define EPS_LOC 1.0e-8


unsigned long int  __my_factorial( const int n ){

  int    i; 
  double factorial=1.0e0;


  if( n<0 ){

    fprintf( stderr, "invalid condition in factorial\n" );
    fflush( stderr );

    factorial=0.0e0;

  }


  for( i=0; i<(n-1); i++ ){

    factorial *= (double)(n-i);

  }

  return (unsigned long int) factorial;

}

unsigned long int __my_double_factorial( const int n ){

  int    i; 
  double double_factorial=1.0e0;


  if( n<0 ){
    
    fprintf( stderr, "invalid condition in double_factorial\n" );
    fflush( stderr );
    
    double_factorial=0.0e0;

  }


  for( i=0; i<(n-1); i +=2 ){ 

    double_factorial *= (double)(n-i);

  }


  return (unsigned long int) double_factorial;

}

unsigned long int __my_choose( const int n, const int k ){

  double dummy=0.0e0;
  int    k_opt;
  int    static info=0;


  if( info || n < 0 || k < 0 || k > n ){

    fprintf( stderr, "ERROR: invalid conditions in choose\n");
    fflush( stderr );
 
    dummy = 0.0e0;

    info =1;

  }



  if( !info ){


    if( 0==k || 0==n || k==n ){

      dummy = 1.0e0;

    }
    else{

      if( k>n/2 ){

	k_opt = n-k; /* WARNING: faster */

      }
      else{

	k_opt = k;

      }


      dummy = (double)( CHOOSE( n-1, k_opt-1 ) + CHOOSE( n-1, k_opt ) ); /* WARNING: recursion */

    }


  } /* end info condistional */

  return (unsigned long int) dummy+EPS_LOC;

}
