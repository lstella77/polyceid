
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

#include "config.h"
#include "PolyCEID_parsing.h"


#define MAX_STRING_LENGTH   128


int __my_parsing_field_search( FILE* fp ){

  /* dummies */
  int  c;
  int  info=0;


  while( !info ){

    if( SKIP_BLANKS( fp ) || feof( fp ) ){
      info=1;
      break;
    }

    c = fgetc( fp );
    /* fp is advanced */
    
    if( ( c == '#' ) || ( c == '!' ) ){

      if( SKIP_LINE( fp ) ){
	info=1;
	break;
      }
      
    }
    else{
      
      fseek( fp, -1, SEEK_CUR ); /* 1 char rewind */
      break;
    
    }

  }


  return info;  

}

int __my_parsing_field_read( FILE* fp, char* string ){

  /* dummies */
  int info=0;


  if( FIELD_SEARCH( fp ) || POINTER_CHECK( string, __my_parsing_field_read ) ){
    fprintf( stderr, "No field found.\n" );
    fflush( stderr );
    info=1;
  }

  if( !info ){

    if( GET_STRING( fp, string ) ) info=1;

  }

  return info;

}

int __my_parsing_multifield_read( FILE* fp, int n, char** sp ){

  /* dummies */
  char* string;
  int   i;
  int   info=0;


  for( i=0; i<n; i++){

    if( FIELD_SEARCH( fp ) ){
      if( i==0) 
	fprintf( stderr, "1st field not found.\n" );
      else if( i==1 )
	fprintf( stderr, "2nd field not found.\n" );
      else if( i==2 )
	fprintf( stderr, "3rd field not found.\n" );
      else
	fprintf( stderr, "%ith field not found.\n", i );
      fflush( stderr );
      info=1;
    }

    //      fprintf( stderr, "I pass here (1)\n" );

    string = *sp++;

    if( POINTER_CHECK( string, __my_parsing_multifield_read ) ){
      info=1;
      break;
    }
    else{

      /* This check must be here */
      if( feof( fp ) ){
	fprintf( stderr, "unexpected EOF reached!\n" );
	info=1;
	break;
      }
      
      //	fprintf( stderr, "I pass here (3)\n" );
      
      if( GET_STRING( fp, string ) ) info=1;
      
    }

  }

  return info;

}

int __my_parsing_field_print( FILE* fp, char* string ){

  /* dummies */
  int info=0;


  if( FILE_CHECK( fp, __my_parsing_field_write ) || 
      POINTER_CHECK( string, __my_parsing_field_write ) ){
    info=1;
  }

  if( !info ){

    fprintf( fp, "%s\n", string );

  }


  return 0;

}

int __my_parsing_multifield_print( FILE* fp, int n, char** sp ){

  /* dummies */
  char* string;
  int   i;
  int   info=0;


  if( FILE_CHECK( fp, __my_parsing_multifield_print )  ){
    info=1;
  }

  if( !info ){

    for( i=0; i<n; i++){

      string = *sp++;

      if( POINTER_CHECK( string, __my_parsing_multifield_print ) ){
	info=1;
	break;
      }
      else{

	fprintf( fp, "%s ", string );

      }

    }

    fprintf( fp, "\n" );

  }


  return info;

}

int __my_parsing_field_check( FILE* file_p, char* string, char* buffer ){

  /* dummies */
  char field[MAX_STRING_LENGTH];
  int  info=0;


  if( FIELD_READ( file_p, field ) ) info=1;

  if( !info ){

    if( strcmp( field, string ) != 0 ){

      fprintf(stderr, "MATCHING ERROR: field %s expected, field %s found instead.\n", string, field );
      fflush( stderr );

      info=1;

    }
    else{

      if( FIELD_READ( file_p, buffer ) ) info=1;

    }

  }

  return info;

}

int __my_parsing_multifield_check( FILE* file_p, char* string, int n, char** buffer ){

  /* dummies */
  char field[MAX_STRING_LENGTH];
  int  info=0;


  if( FIELD_READ( file_p, field ) ) info=1;

  if( !info ){

    if( strcmp( field, string ) != 0 ){

      fprintf(stderr, "MATCHING ERROR: field %s expected, field %s found instead.\n", string, field );
      fflush( stderr );

      info=1;

    }
    else{

      if( MULTIFIELD_READ( file_p, n, buffer ) ) info=1;

    }

  }

  return info;

}

/* utilities */

int __my_skip_blanks( FILE* fp ){

  /* dummies */
  int  c;
  int  info=0;


  if( FILE_CHECK( fp, __my_skip_blanks ) ){
    fprintf( stderr, "No file found!\n");
    fflush( stderr );
    info=1;
  }

  while( !feof( fp ) && !info ){

    if( ferror( fp ) ){
      fprintf( stderr, "File error found!\n");
      fflush( stderr );
      info=1;
      break;
    }
    
    c = fgetc( fp );
    /* fp is advanced */
    
    if( !isspace( c ) ){
      fseek( fp, -1, SEEK_CUR ); /* 1 char rewind */
      break;
    }

  }


  return info;

}

int __my_skip_line( FILE* fp ){

  /* dummies */
  int  c;
  int  info=0;


  if( FILE_CHECK( fp, __my_skip_line ) ){ 
    fprintf( stderr, "No file found!\n");
    fflush( stderr );
    info=1;
  }

  while( !feof(fp) && !info ){

    if( ferror( fp ) ){
      fprintf( stderr, "File error found!\n");
      fflush( stderr );
      info=1;
      break;
    }
    
    c = fgetc( fp );
    /* fp is advanced */
    
    if( ( c == '\n' ) || ( c == '\0' ) ) break;
    
  }


  return info;

}
int __my_get_string( FILE* fp, char* string ){

  /* dummies */
  int c;
  int length=0; 
  int info=0;


  if( FILE_CHECK( fp, __my_skip_line ) || SKIP_BLANKS( fp ) ) info=1;

  c = fgetc( fp );
  /* fp is advanced */

  while( !isspace( c ) && !feof(fp) && !info ){

    if( ferror( fp ) ){
      info=1;
      break;
    }
    
    /* copy */
    string[length] = (char)c;
    /*increase length */
    length++;

    c = fgetc( fp );
    /* fp is advanced */     


  }

  /* finish the string */
  string[length] = '\0';


  return info;

}
