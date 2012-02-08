
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
/* FILE: PolyCEID_field.h
 *
 * It contains function for parsing
 * The input file is made by string fields
 * followed by variables
 *
 */                

#include "../input/PolyCEID_assign_variable.h"

/*****************************
 **********FUNCTIONS**********
 *****************************/

/* FUNCTION: _assign_int_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 *            int     default value
 *
 * RETURN:    int
 *
 * search for the int variable in the list and assign 
 * its value
 *
 */
int _assign_int_variable( list_p list_p, char* name, int n, int* value_p, int* default_value_p ){

  /* dummies */
  entry_p  entry_found_p=NULL;  
  entry_p  dummy_entry_p=NULL;  
  int      i;
  int      info=0;


  // search for the variable
  entry_found_p = LIST_SEARCH( list_p, name );

  if( !entry_found_p ){

    if( !default_value_p ){

      fprintf( PolyCEID_STDERR, "ERROR: variable %s must be explicitly decalred\n", name );

      fflush( PolyCEID_STDERR );

      info=1;

    }        
    else{

      S_ECHO( "Using default value for variable:", name );

      for( i=0; i<n; i++ ){

        value_p[ i ] = default_value_p[ i ];      

      }

    }

  }
  else{

    S_ECHO( "Extracted value for variable:", name );
    
    if( !entry_found_p->sublist_p ){

      fprintf( PolyCEID_STDERR, "ERROR: no value to extract from variable %s\n", entry_found_p->field.name );

      fflush( PolyCEID_STDERR );

      info=1;

    }        
    else{

      if( entry_found_p->sublist_p->N_entries < n ){

        fprintf( PolyCEID_STDERR, "ERROR: too few values to extract from variable %s\n", entry_found_p->field.name );

        fflush( PolyCEID_STDERR );

        info=1;

      }
      else{

        dummy_entry_p = entry_found_p->sublist_p->first_entry_p;

        for( i=0; i<n; i++ ){

          value_p[ i ] = atoi( dummy_entry_p->field.name );

          dummy_entry_p = dummy_entry_p->next_entry_p;

        }        

        // cut the entry
        ENTRY_CUT( entry_found_p );

      }        

    }

  }

  
  return info;

}

/* FUNCTION: _assign_double_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 *            double     default value
 *
 * RETURN:    int
 *
 * search for the double variable in the list and assign 
 * its value
 *
 */
int _assign_double_variable( list_p list_p, char* name, int n, double* value_p, double* default_value_p ){

  /* dummies */
  entry_p  entry_found_p=NULL;  
  entry_p  dummy_entry_p=NULL;  
  int      i;
  int      info=0;


  // search for the variable
  entry_found_p = LIST_SEARCH( list_p, name );

  if( !entry_found_p ){

    if( !default_value_p ){

      fprintf( PolyCEID_STDERR, "ERROR: variable %s must be explicitly decalred\n", name );

      fflush( PolyCEID_STDERR );

      info=1;

    }        
    else{

      S_ECHO( "Using default value for variable:", name );

      for( i=0; i<n; i++ ){

        value_p[ i ] = default_value_p[ i ];      

      }

    }

  }
  else{

    S_ECHO( "Extracted value for variable:", name );
    
    if( !entry_found_p->sublist_p ){

      fprintf( PolyCEID_STDERR, "ERROR: no value to extract from variable %s\n", entry_found_p->field.name );

      fflush( PolyCEID_STDERR );

      info=1;

    }        
    else{

      if( entry_found_p->sublist_p->N_entries < n ){

        fprintf( PolyCEID_STDERR, "ERROR: too few values to extract from variable %s\n", entry_found_p->field.name );

        fflush( PolyCEID_STDERR );

        info=1;

      }
      else{

        dummy_entry_p = entry_found_p->sublist_p->first_entry_p;

        for( i=0; i<n; i++ ){

          value_p[ i ] = atof( dummy_entry_p->field.name );

          // fprintf( stdout, "%g\n", value_p[ i ] );

          dummy_entry_p = dummy_entry_p->next_entry_p;

        }        

        // cut the entry
        ENTRY_CUT( entry_found_p );

      }        

    }

  }

  
  return info;

}


/* FUNCTION: _assign_string_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 *            char*  default value
 *
 * RETURN:    char*
 *
 * search for the char* variable in the list and assign 
 * its value
 *
 */
int _assign_string_variable( list_p list_p, char* name, char* value_p, char* default_value_p ){

  /* dummies */
  entry_p  entry_found_p=NULL;  
  entry_p  dummy_entry_p=NULL;
  int      info=0;


  // search for the variable
  entry_found_p = LIST_SEARCH( list_p, name );

  if( !entry_found_p ){

    S_ECHO( "Using default value for variable:", name );

    value_p[ 0 ] = '\0';
    strncat( value_p, dummy_entry_p->field.name, MAX_NAME_LENGTH );

  }
  else{

    S_ECHO( "Extracted value for variable:", name );

    if( !entry_found_p->sublist_p ){

      fprintf( PolyCEID_STDERR, "ERROR: no value to extract from variable %s\n", entry_found_p->field.name );

      fflush( PolyCEID_STDERR );

      info=1;

    }        
    else{      

      dummy_entry_p = entry_found_p->sublist_p->first_entry_p; 

      value_p[ 0 ] = '\0';
      strncat( value_p, dummy_entry_p->field.name, MAX_NAME_LENGTH );

      // cut the entry
      ENTRY_CUT( entry_found_p );

    }  

  }


  return info;

}

/* FUNCTION: _assign_flag_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 * RETURN:    unsigned short int
 *
 * search for the flag variable in the list and assign 
 * its value
 *
 */
unsigned short int _assign_flag_variable( list_p list_p, char* name ){

  /* dummies */
  entry_p             entry_found_p;
  unsigned short int  flag=0;


  // search for the variable
  entry_found_p = LIST_SEARCH( list_p, name );

  if( !entry_found_p ){

    flag=0; // default

  }
  else{

    S_ECHO( "Set flag", name );

    flag=1;

    // cut the entry
    ENTRY_CUT( entry_found_p );

  }


  return flag;

}
