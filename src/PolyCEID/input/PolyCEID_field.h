
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

/* FILE: PolyCEID_field.h
 *
 * It contains function for parsing
 * The input file is made by string fields
 * followed by variables
 *
 */

#ifndef _PolyCEID_FIELD_NEW_
#define _PolyCEID_FIELD_NEW_

// included headers
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "../utils/PolyCEID_error.h"

/*****************************
 ************MACROS***********
 *****************************/

// Max length string (fortran like)
#define MAX_NAME_LENGTH   80

// convert a field to a string
#define FIELD_NEW_TO_STRING( field )   (field).name  

// convert a field to a int
#define FIELD_NEW_TO_INT( field )      atoi( (field).name )  

// convert a field to a double
#define FIELD_NEW_TO_DOUBLE( field )   atof( (field).name )  

/*****************************
 **********STRUCTURES*********
 *****************************/

/* TYPE: field
 *
 */
struct _PolyCEID_field
{

  char name[MAX_NAME_LENGTH];

};

typedef  struct  _PolyCEID_field     field;

typedef  struct  _PolyCEID_field*    field_p;

typedef  struct  _PolyCEID_field**   field_pp;

/*****************************
 **********FUNCTIONS**********
 *****************************/

/* FUNCTION: _field_read
 *
 * ARGUMENTS: FILE* fp
 *            field_p field_p
 *
 * read the field name after the field
 * has been found
 *
 * NOTE: it must not be preprocessed!
 *
 */
int _field_read( FILE* fp, field_p field_p );

#define FIELD_NEW_READ( fp, field_p )   _field_read( (fp), (field_p) )

/* FUNCTION: _field_print
 *
 * ARGUMENTS: FILE* fp
 *
 *            field_p field_p
 *
 *
 */
int _field_print( FILE* fp, field_p field_p );

#define FIELD_NEW_PRINT( fp, field_p )   FUNCTION_CHECK_NEW( _field_print( (fp), (field_p) ) )

/* FUNCTION: _is_check
 *
 * ARGUMENTS: char*   name    
 *
 *            field_p field_p
 *
 * check whether the field argument is "name"
 *
 * Note: this function must not be preprocessed!
 *
 */
int _is_field( char* name, field_p field_p );

#define IS_FIELD_NEW( name, field_p )   _is_field( (name), (field_p) )

/* FUNCTION: _field_assign
 *
 * ARGUMENTS:  char*    name
 *
 *             field_p  field_p
 *
 * name --> field_p->name
 * It checks the field length, and truncates
 * if necessary.
 *
 */
int _field_assign( char* name, field_p field_p );
	
#define FIELD_NEW_ASSIGN( name, field_p )   FUNCTION_CHECK_NEW(  _field_assign( (name), (field_p) ) )	

/* FUNCTION: _field_copy
 *
 * ARGUMENTS:  field_p   field_in_p
 *
 *             field_p   field_out_p
 *
 * field_in --> field_out
 * It checks the field length, and truncates
 * if necessary.
 *
 */
int _field_copy( field_p field_in_p, field_p field_out_p );
	
#define FIELD_NEW_COPY( field_in_p, field_out_p )   FUNCTION_CHECK_NEW(  _field_copy( (field_in_p), (field_out_p) ) )	
	
/*****************************
 *********UTLITIES************
 *****************************/

/* FUNCTION: _field_search
 *
 * ARGUMENTS: FILE* fp
 *
 * search for a field in the file fp
 *
 * NOTE: it must not be preprocessed!
 *
 */
int _field_search( FILE* fp );

#define FIELD_NEW_SEARCH( fp )   _field_search( (fp) )

/* FUNCTION _get_field
 *
 * ARGUMENTS FILE* fp
 *           field_p field_p
 *
 * Put string of the file fp
 * into variable string
 */
int _get_field( FILE* fp, field_p field_p );

#define GET_FIELD_NEW( fp, field_p )   FUNCTION_CHECK_NEW( _get_field( (fp), (field_p) ) )

#endif /* _PolyCEID_FIELD_NEW_ */
