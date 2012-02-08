
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
#include "PolyCEID_utilities.h"
#include <stdio.h>

/*********************
  FUNCTIONS & MACROS
*********************/

int pbc_int( int site, int periodicity ){

  int site_pbc;


  site_pbc = site %periodicity;

  if( site_pbc < 0)
    site_pbc = pbc_int( site + periodicity, periodicity );


  return site_pbc;

}

double pbc_diff( double diff, double periodicity ){

  double diff_pbc;

  if( diff < -0.5 *periodicity ){

    diff_pbc = pbc_diff( diff + periodicity, periodicity );

  }
  else if( diff >= 0.5 *periodicity ){
  
    diff_pbc = pbc_diff( diff - periodicity, periodicity );

  }        
  else{

    diff_pbc = diff;

  }        


  return diff_pbc;

}

