
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
/* FILE: PolyCEID_error.c
 *
 * contains the error handling
 *
 */



// header
#include "PolyCEID_error.h"



// global variable
int global_info=0;   //default value, all is fine



/*****************************
 **********FUNCTIONS**********
 *****************************/



/* FUNCTION: _pointer_check
 *
 * ARGUMENTS: void* p        pointer
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
int _pointer_check( void* p, char* file, int line ){


  /* dummies */
  int info=0;


  if( !(p) ) {

    fprintf( PolyCEID_STDERR, "ERROR: pointer %p pointing to NULL in file '%s', line %d\n", p, file, line );

    fflush( PolyCEID_STDERR );

    global_info=1;

    info=1;

  }

  
  return info;

}




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
int _function_check( int p, char* string, char* file, int line ){


  /* dummies */
  int info=0; // WARNING: default value, no errors


  info=p;

  if( info ) {

    fprintf( PolyCEID_STDERR, "ERROR: function '%s' returns %d in file '%s', line %d\n", string, info, file, line );

    fflush( PolyCEID_STDERR );

    global_info=1;;  

  }


  return info;

}



/* FUNCTION: _force_global_info
 *
 * force global_info to 1
 *
 */
void _force_global_info( void ){

  global_info=1;

}



/* FUNCTION: _reset_global_info
 *
 * reset global_info to 0
 *
 */
void _reset_global_info( void ){

  global_info=0;

}	



/* FUNCTION: _check_global_info
 *
 * return global_info
 *
 */
int _check_global_info( void ){

  return global_info;

}	
