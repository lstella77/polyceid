
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
#include "PolyCEID_start_file.h"


//------------------------------------------

/* READ */

int read_start_file( constants_p constants_p, state_p state_p, config_p config_p ){

  /* dummies */
  char  buffer[MAX_STRING_LENGTH];
  FILE* start_file_p;
  int   info=0;


  /* file name */
  if( !strcpy(  buffer, "restart.dat" ) ) info=1;

  /* open */
  start_file_p = fopen( buffer, "rb");

  if( FILE_CHECK( start_file_p, output_start_file ) ){
    fprintf( stderr, "ERROR occured when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_READ( start_file_p, *constants_p ) )  info=1;
  // if( STATE_READ( start_file_p, *state_p ) )          info=1;
  if( CONFIG_READ( start_file_p, *config_p ) )        info=1;

  /* close */
  if( FILE_CHECK( start_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occured when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  fclose( start_file_p );


  return info;

}


//------------------------------------------

/* PRINT */

int print_start_file( const constants constants, const state state, const config config ){

  /* dummies */
  char  buffer[MAX_STRING_LENGTH];
  FILE* start_file_p;
  int   info=0;


  /* file name */
  if( !strcpy(  buffer, "restart.dat" ) ) info=1;

  /* open */
  start_file_p = fopen( buffer, "wb");

  if( FILE_CHECK( start_file_p, output_start_file ) ){
    fprintf( stderr, "ERROR occured when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_PRINT( start_file_p, constants ) )  info=1;
  // if( STATE_PRINT( start_file_p, state ) )          info=1;
  if( CONFIG_PRINT( start_file_p, config ) )        info=1;

  /* close */
  if( FILE_CHECK( start_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occured when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  fclose( start_file_p );


  return info;

}


//------------------------------------------
