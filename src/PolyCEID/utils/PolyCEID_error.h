
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

/* FILE: PolyCEID_error.c
 *
 * contains the error handling
 *
 */

#ifndef _PolyCEID_ERROR_
#define _PolyCEID_ERROR_

// included headers
#include <stdlib.h>
#include <stdio.h>

/*****************************
 ************MACROS***********
 *****************************/

#define   PolyCEID_STDERR    stderr

#define   PolyCEID_STDOUT    stdout

#define   PASSING_HERE       fprintf( PolyCEID_STDOUT, "# [TRACING] passing line %d of file %s\n", __LINE__, __FILE__ )

#ifdef _DEBUG_

#define   S_ECHO( msg, x )   fprintf( PolyCEID_STDOUT, "# [S_ECHO] " msg "\n" ); \
                             fprintf( PolyCEID_STDOUT, "# [S_ECHO] " #x  " = %s\n", x ); \
                             fprintf( PolyCEID_STDOUT, "# [S_ECHO] at line %d of file %s\n", __LINE__, __FILE__  ); \
			     fprintf( PolyCEID_STDOUT, "\n" );

#define   D_ECHO( msg, x )   fprintf( PolyCEID_STDOUT, "# [D_ECHO] " msg "\n" ); \
                             fprintf( PolyCEID_STDOUT, "# [D_ECHO] " #x  " = %d\n", x ); \
                             fprintf( PolyCEID_STDOUT, "# [D_ECHO] at line %d of file %s\n", __LINE__, __FILE__  ); \
			     fprintf( PolyCEID_STDOUT, "\n" );

#define   F_ECHO( msg, x )   fprintf( PolyCEID_STDOUT, "# [F_ECHO] " msg "\n" ); \
                             fprintf( PolyCEID_STDOUT, "# [F_ECHO] " #x  " = %f\n", x ); \
                             fprintf( PolyCEID_STDOUT, "# [F_ECHO] at line %d of file %s\n", __LINE__, __FILE__  ); \
			     fprintf( PolyCEID_STDOUT, "\n" );

#else /* _DEBUG_ */

#define   S_ECHO( msg, x )

#define   D_ECHO( msg, x )

#define   F_ECHO( msg, x )

#endif /* _DEBUG */

/*****************************
 **********FUNCTIONS**********
 *****************************/

/* FUNCTION: _pointer_check
 *
 * ARGUMENTS: void* p        pointer
 *
 *            char* file     file name
 *
 *            int   line     file line
 *
 * check a pointer
 *
 */
int _pointer_check( void* p, char* file, int line );

#ifndef _NO_CHECK_

#define POINTER_CHECK_NEW( p )    _pointer_check( (p), (__FILE__), (__LINE__) )

#else /* _NO_CHECK_ */

#define POINTER_CHECK_NEW( p )   (0)

#endif /* _NO_CHECK_ */

/* FUNCTION: _function_check
 *
 * ARGUMENTS: int   p        function returned value
 *
 *            char* string   pointer name
 *
 *            char* file     file name
 *
 *            int   line     file line
 *
 * check a pointer
 *
 */
int _function_check( int p, char* string, char* file, int line );

#ifndef _NO_CHECK_

#define FUNCTION_CHECK_NEW( p )    ( !CHECK_GLOBAL_INFO ? _function_check( (p), (#p), (__FILE__), (__LINE__) ) : 1 )

#else /* _NO_CHECK_ */

#define FUNCTION_CHECK_NEW( p )   (p)

#endif /* _NO_CHECK_ */

/* FUNCTION: _force_global_info
 *
 * force global_info to 1
 *
 */
void _force_global_info( void );

#define FORCE_GLOBAL_INFO    _force_global_info()

/* FUNCTION: _reset_global_info
 *
 * reset global_info to 0
 *
 */
void _reset_global_info( void );

#define RESET_GLOBAL_INFO    _reset_global_info()

/* FUNCTION: _check_global_info
 *
 * return global_info
 *
 */
int _check_global_info( void );

#define CHECK_GLOBAL_INFO    _check_global_info()

#endif /* _PolyCEID_ERROR_ */
