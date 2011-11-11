
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

#ifndef __PolyCEID_SORTING__
#define __PolyCEID_SORTING__


#include "../algebra/PolyCEID_vector.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* __my_qsort */

int __my_rvector_qsort( const rvector, ivector_p );

#define RVECTOR_QSORT( rvec, ivec ) FUNCTION_CHECK( __my_rvector_qsort( (rvec), &(ivec) ), __my_rvector_qsort ) 

//------------------------------------------

/* utilities */

//------------------------------------------

int __my_rvector_qsort_aux( const rvector, ivector_p, int, int );

#define RVECTOR_QSORT_AUX( rvec, ivec, ini, fin ) FUNCTION_CHECK( __my_rvector_qsort_aux( (rvec), &(ivec), (ini), (fin) ), __my_rvector_qsort_aux ) 

//------------------------------------------


#endif /* __PolyCEID_SORTING__ */
