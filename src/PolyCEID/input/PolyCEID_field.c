
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
/* FILE: PolyCEID_field.c
 *
 * It contains function for parsing
 * The input file is made by string fields
 * followed by variables
 *
 */



// header
#include "PolyCEID_field.h"



/*****************************
 **********FUNCTIONS**********
 *****************************/



/* FUNCTION: _field_read
 *
 * ARGUMENTS: FILE*    fp
 *            field_p  field_p
 *
 * read the field name after the field
 * has been found
 *
 */
int _field_read( FILE* fp, field_p field_p ){

  /* dummies */
  int info=0;


  // if a field is found, get it 
  if( FIELD_NEW_SEARCH( fp ) ){

    GET_FIELD_NEW( fp, field_p ); 
    
    info=1;

  }
  

  return info;

}



/* FUNCTION: _field_print
 *
 * ARGUMENTS: FILE* fp
 *
 *            field_p field_p
 *
 *
 */
int _field_print( FILE* fp, field_p field_p ){

  /* dummies */
  int info=0;


  // check for errors
  if( POINTER_CHECK_NEW( fp ) || POINTER_CHECK_NEW( field_p->name ) ){

    fprintf( PolyCEID_STDERR, "ERROR: not able to print the field  [%s] on file [%p]!\n", field_p->name, fp );

    info=1;

  }
  else{

    fprintf( fp, "%s\n", field_p->name );

  }


  return info;

}



/* FUNCTION: _is_field
 *
 * ARGUMENTS: char*   name
 *
 *            field_p field_p
 *
 * check whether the field argument is "name"
 *
 */
int _is_field( char* name, field_p field_p ){

  /* dummies */
  int    info=1; // WARNING: default value is 1 here


  // check field name
  if( strcmp( name, field_p->name ) != 0 ){

    info=0;

  }

  return info;

}



/* FUNCTION: _field_assign
 *
 * ARGUMENTS:  char*     name
 *
 *             field_p   field_p
 *
 * name --> field_p
 * It checks the field length, and truncates
 * if necessary.
 *
 */
int _field_assign( char* name, field_p field_p ){


  /* dummies */
  int i;
  int info=0;

  // if( !strncpy( field_p->name, name, MAX_NAME_LENGTH ) ) info=1;

  for( i=0; i<MAX_NAME_LENGTH; i++ ){

    field_p->name[ i ] = name[ i ];

  }


  return info;

}



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
int _field_copy( field_p field_in_p, field_p field_out_p ){

  /* dummies */
  int info=0;


  if( !strncpy( field_out_p->name, field_in_p->name, MAX_NAME_LENGTH ) ) info=1;


  return info;

}



/*****************************
 *********UTLITIES************
 *****************************/



/* FUNCTION: _field_search
 *
 * ARGUMENTS: FILE* fp
 *
 * search for a field in the file fp
 *
 */
int _field_search( FILE* fp ){

  /* dummies */
  int  c;
  int  info=0;

  
  // pointer check
  if( POINTER_CHECK_NEW( fp ) ){ 
  
    fprintf( PolyCEID_STDERR, "ERROR: file pointer [%p] not allocated yet!\n", fp );

    fflush( PolyCEID_STDERR );

  }
  else{

    // till the end of the file is reached
    while( EOF != ( c = getc( fp ) ) ){

      // skip blanks
      while( isspace( c ) ){

        if( EOF == ( c = getc( fp ) ) ) break; 
	      
      }

      // if EOF exit
      if( EOF == c ){

	break;

      }
      else if( '#' == c || '!' == c ){ 
	
	// skip comments
        while( EOF != ( c = getc( fp ) ) ){

          if( '\n' == c || '\0' == c ) break;

        }

      }
      else{

        // The field has been found
        // 1 char rewind WARNING: important!
        fseek( fp, -1, SEEK_CUR );

        info=1;

        break;
    
      }

    } /* end while*/

  } /* end pointer check */


  return info;  

}



/* FUNCTION _get_field
 *
 * ARGUMENTS FILE* fp
 *           char* string
 *
 * Put string of the file fp
 * into variable string
 */
int _get_field( FILE* fp, field_p field_p ){

  /* dummies */
  char  c;
  int   i;
  int   length=0; 
  int   info=0;


  // check if the file exist
  if( POINTER_CHECK_NEW( fp ) ){
	  
   info=1;

  }
  else{
    
    // till the end of the file is reached
    while( EOF != ( c = getc( fp ) ) ){

      // field's name ends when a blanck char is found
      if( isspace( c ) ){
	      
	break;

      }  
      else{

        if( length < MAX_NAME_LENGTH ){

          // char copy
          field_p->name[ length ] = (char)c;

          // string length is increased
          length++;    

        }
        else if( MAX_NAME_LENGTH == length ){

          // ignore the rest of the field
          while( EOF != ( c = getc( fp ) ) ){

            if( isspace( c ) ) break;          

	  }        

          // trim if it is needed
          for( i=length; i<MAX_NAME_LENGTH; i++ ){

	    field_p->name[ i ] = '\0';

	  }

          fprintf( PolyCEID_STDERR, "WARNING [GET]: field name is longer than %d chars.\n", MAX_NAME_LENGTH );
          fprintf( PolyCEID_STDERR, "WARNING [GET]: field name truncated to '%s'.\n", field_p->name );

          fflush( PolyCEID_STDERR );

          break;

        } /* end length check */   

      } /* end of isspace check */

    } /* end while */

  }

  if( length < MAX_NAME_LENGTH ){

    // trim if it is needed
    for( i=length; i<MAX_NAME_LENGTH; i++ ){

      field_p->name[ i ] = '\0';

    }

  }  


  return info;

}
