
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
#include "PolyCEID_start_file.h"


FILE* start_file_p;  /* private variable */


//------------------------------------------

/* OPEN */

int open_start_file( const constants constants ){

  /* constants */
  char*     output_label;
  /* dummies */
  char      buffer[MAX_STRING_LENGTH];
  int       info=0;


  output_label = (char*)constants.output_label;

  /* file name */
  if( !strcpy(  buffer, "start_file_" ) ) info=1;
  if( !strcat(  buffer, output_label )  ) info=1;
  if( !strcat(  buffer, ".dat" )        ) info=1;


  start_file_p = fopen( buffer, "wt");

  if( FILE_CHECK( start_file_p, output_start_file ) ){
    fprintf( stderr, "ERROR occured when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }


  return info;

}

//------------------------------------------

/* CLOSE */

int close_start_file( const constants constants ){

  /* constants */
  char*     output_label;
  /* dummies */
  char      buffer[MAX_STRING_LENGTH];
  int       info=0;


  output_label = (char*)constants.output_label;

  /* file name */
  if( !strcpy(  buffer, "start_file_" ) ) info=1;
  if( !strcat(  buffer, output_label )  ) info=1;
  if( !strcat(  buffer, ".dat" )        ) info=1;


  if( FILE_CHECK( start_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occured when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  fclose( start_file_p );


  return info;

}

//------------------------------------------

/* READ */

int read_start_file( constants constants, state state, config config ){

  int info=0;


  if( CONSTANTS_READ( start_file_p, constants ) )  info=1;
  if( STATE_READ( start_file_p, state ) )          info=1;
  if( CONFIG_READ( start_file_p, config ) )        info=1;


  return info;

}


//------------------------------------------

/* PRINT */

int print_start_file( const constants constants, const state state, const config config ){

  int info=0;


  if( CONSTANTS_PRINT( start_file_p, constants ) )  info=1;
  if( STATE_PRINT( start_file_p, state ) )          info=1;
  if( CONFIG_PRINT( start_file_p, config ) )        info=1;


  return info;

}


//------------------------------------------
