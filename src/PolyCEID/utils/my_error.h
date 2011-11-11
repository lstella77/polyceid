
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

#ifndef __MY_ERROR__
#define __MY_ERROR__

#include <stdlib.h>
#include <stdio.h>


int   __my_file_check( void*, char*, char*, int );

#ifdef __DEBUG__

#define FILE_CHECK( p, string )    __my_file_check( (void*)(p), #string, __FILE__, __LINE__  )

#else /* no debug */ 

#define FILE_CHECK( p, string )    0 //((void)0)

#endif


int   __my_pointer_check( void*, char*,  char*, int );

#ifdef __DEBUG__

#define POINTER_CHECK( p, string )  __my_pointer_check( (void*)(p), #string, __FILE__, __LINE__ )

#else /* no debug */ 

#define POINTER_CHECK( p, string )  0 //((void)0)

#endif


int   __my_function_check( int, char*, char*, int );

#ifdef __DEBUG__

#define FUNCTION_CHECK( p, string )  __my_function_check( (p), #string, __FILE__, __LINE__ )

#else /* no debug */ 

#define FUNCTION_CHECK( p, string )  (p)

#endif

#endif /* __MY_ERROR__ */
