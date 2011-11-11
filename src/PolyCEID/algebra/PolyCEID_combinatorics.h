
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

#ifndef __PolyCEID_COMBINATORICS__
#define __PolyCEID_COMBINATORICS__

#include <stdio.h>

/*
  Notice that outputs are doubles 
  because of internal precision
  requirements.
*/


unsigned long int __my_factorial( const int n );

#define FACTORIAL( n ) __my_factorial( (n) )

unsigned long int  __my_double_factorial( const int n );

#define DOUBLE_FACTORIAL( n ) __my_double_factorial( (n) )

unsigned long int __my_choose( const int n, const int k );

#define CHOOSE( n, k )  __my_choose( (n), (k) )


#endif /* __PolyCEID_COMBINATORICS__ */
