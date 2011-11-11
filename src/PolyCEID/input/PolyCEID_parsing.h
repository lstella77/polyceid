
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

#ifndef __PolyCEID_PARSING__
#define __PolyCEID_PARSING__


#include "../utils/my_error.h"
#include <ctype.h>
#include <string.h>


int __my_parsing_field_search( FILE* );

#define FIELD_SEARCH( fp )  FUNCTION_CHECK(  __my_parsing_field_search( (fp) ),  __my_parsing_field_search )

int __my_parsing_field_read( FILE*, char* );

#define FIELD_READ( fp, s )  FUNCTION_CHECK(  __my_parsing_field_read( (fp), (s) ),  __my_parsing_field_read )

int __my_parsing_multifield_read( FILE*, int, char** );

#define MULTIFIELD_READ( fp, n, sp )  FUNCTION_CHECK(  __my_parsing_multifield_read( (fp), (n), (sp) ),  __my_parsing_multifield_read )

int __my_parsing_field_print( FILE*, char* );

#define FIELD_PRINT( fp, s )  FUNCTION_CHECK(  __my_parsing_field_print( (fp), (s) ),  __my_parsing_field_print )

int __my_parsing_multifield_print( FILE*, int, char** );

#define MULTIFIELD_PRINT( fp, n, sp )  FUNCTION_CHECK(  __my_parsing_multifield_print( (fp), (n), (sp) ),  __my_parsing_multifield_print )

int __my_parsing_field_check( FILE*, char*, char* );

#define FIELD_CHECK( fp, s, b ) FUNCTION_CHECK( __my_parsing_field_check( (fp), (s), (b) ), __my_parsing_field_check )

int __my_parsing_multifield_check( FILE*, char*, int, char** );

#define MULTIFIELD_CHECK( fp, s, n, bp ) FUNCTION_CHECK( __my_parsing_multifield_check( (fp), (s), (n), (bp) ), __my_parsing_multifield_check )

/* utilities */

int __my_skip_blanks( FILE* );

#define SKIP_BLANKS( fp ) FUNCTION_CHECK( __my_skip_blanks( (fp) ), __my_skip_blanks )

int __my_skip_line( FILE* );

#define SKIP_LINE( fp ) FUNCTION_CHECK( __my_skip_line( (fp) ), __my_skip_line )

int __my_get_string( FILE*, char* );

#define GET_STRING( fp, s ) FUNCTION_CHECK( __my_get_string( (fp), (s) ), __my_get_string )

#endif /* __PolyCEID_PARSING__ */
