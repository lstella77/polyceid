
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
#include "PolyCEID_vector.h"


/* complex vectors */

int __my_vector_allocate( const int dim, vector_p vec_p ){


  vec_p->vector_dim  = dim;
  vec_p->vector      = (complex*) calloc( (size_t) dim, sizeof(complex) );

  if( (VECTOR_CHECK( *vec_p )) ) return 1;

  return 0;

}

int       __my_vector_free( vector_p vec_p){

  if( (VECTOR_CHECK( *vec_p )) ) return 1;

  vec_p->vector_dim = 0;

  free( vec_p->vector );

  vec_p->vector = NULL;

  return 0;

}

int       __my_vector_copy( vector_p vec1_p, const vector vec2  ){

  /* dummies */
  int dim, i;


  if( VECTOR_CONSINSTENCY( *vec1_p, vec2 ) ) return 1; 

  dim = vec2.vector_dim;

  for( i=0; i<dim; i++ ){
    //CMPLX_COPY( vec1_p->vector[i], vec2.vector[i] ); 
    vec1_p->vector[i] = vec2.vector[i];
  }


  return 0;

}

int       __my_vector_compare( const vector vec1, const vector vec2  ){

  /* dummies */
  int dim, i;
  int info=0;


  if( VECTOR_CONSINSTENCY( vec1, vec2 ) ) info=1; 

  dim = vec2.vector_dim;

  for( i=0; i<dim; i++ ){

    if( CMPLX_COMPARE( vec1.vector[ i ], vec2.vector[ i ] ) ){

      fprintf( stderr, "ERROR: the %d components are not equal\n", i );
      fflush( stderr );

      info=1;

      break;

    }

  }


  return info;

}

int __my_vector_read( FILE* file_p, vector_p vec_p ){

  /* dummies */
  int dim, i;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_vector_print ) ) || (VECTOR_CHECK( *vec_p ) ) ) info=1;

  if( fread( ( void* ) &dim, sizeof( dim ), 1, file_p ) < 1 ) info=1;

  if( dim != vec_p->vector_dim ) info=1;

  for( i=0; i<dim; i++){

    CMPLX_READ( file_p, vec_p->vector[ i ] );

  }


  return info;

}

int __my_vector_print( FILE* file_p, const vector vec ){

  /* dummies */
  int i, dim;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_vector_print ) ) || ( VECTOR_CHECK( vec ) ) ) info=1;

  dim = vec.vector_dim;

  if( fwrite( ( const void* ) &dim, sizeof( dim ), 1, file_p ) < 1 ) info=1;

  for( i=0; i<dim; i++ ){

    CMPLX_PRINT( file_p, vec.vector[ i ] );

  }


  return info;

}

int       __my_vector_print_plus( FILE* file_p, const vector vec ){

  int i, dim;


  if( (FILE_CHECK( file_p, __my_vector_print )) || (VECTOR_CHECK( vec )) ) return 1;

  dim = vec.vector_dim;

  fprintf( file_p, "# vector dim = %d\n", dim );

  for( i=0; i<dim; i++){

    CMPLX_PRINT_PLUS( file_p, vec.vector[i] );
    fprintf( file_p, "\n");

  }
  fflush( file_p );

  return 0;

}

int       __my_vector_check( const vector vec ){

  if( !(vec.vector_dim) || !(vec.vector) ){
    fprintf( stderr, "VECTOR CHECK ERROR.\n" );
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }

}

int       __my_vector_consinstency( const vector vec1, const vector vec2 ){

  if( (VECTOR_CHECK( vec1 )) || (VECTOR_CHECK( vec2 )) || (vec1.vector_dim != vec2.vector_dim )){
    fprintf( stderr, "VECTOR CONSISTENCY ERROR.\n" );
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }

}

int    __my_vector_set_to_zero( vector_p vec_p ){

  int i, dim;


  if( (VECTOR_CHECK( *vec_p )) ) return 1;

  dim = vec_p->vector_dim;

  for( i=0; i<dim; i++ ){
    //CMPLX_COPY( vec_p->vector[i], CMPLX_ZERO );
    //vec_p->vector[i] = CMPLX_ZERO; /* WARNING! That could be dangerous! */
    CMPLX_ZERO_VAR( vec_p->vector[i] );
  }

  return 0;

}

int    __my_vector_set_to_unit( vector_p vec_p ){

  int i, dim;


  if( (VECTOR_CHECK( *vec_p )) ) return 1;

  dim = vec_p->vector_dim;

  for( i=0; i<dim; i++ ){
    //CMPLX_COPY( vec_p->vector[i], CMPLX_UNIT );
    vec_p->vector[i] = CMPLX_UNIT; /* WARNING! That could be dangerous! */
  }

  return 0;

}

/* ----------------------------------------------- */

int     __my_vector_array_allocate( int dim, int dim_vec, vector_array_p ap ){

  /* dummies */
  int i;
  int info=0;


  ap->dim = dim;

  ap->array = (vector_p) calloc( (size_t)(dim), sizeof(vector) );

  if( POINTER_CHECK( ap->array, __my_vector_array_allocate ) ) info=1;

  if( !info ){
    
    for( i=0; i<dim; i++ ){

      if( VECTOR_ALLOCATE( dim_vec, ap->array[i] ) ) info=1;

    }

  }


  return info;

}

int     __my_vector_array_free( vector_array_p ap ){

  /* dummies */
  int i, dim;
  int info=0;


  if( POINTER_CHECK( ap->array, __my_vector_array_free ) ) info=1;

  if( !info ){

    dim = ap->dim;

    for( i=0; i<dim; i++ ){

      if( VECTOR_FREE( ap->array[i] ) ) info=1;

    }
  
    ap->dim = 0;

    free( ap->array );

    ap->array = NULL;

  }


  return info;

}

int     __my_vector_array_copy( vector_array_p ap, const vector_array az ){

  /* dummies */
  int dim, i;
  int info=0;

  if( POINTER_CHECK( ap->array, __my_vector_array_copy ) ||
      ( ap->dim != az.dim )
      ) info=1;

  dim = ap->dim;

  for( i=0; i<dim; i++){
    if( VECTOR_CONSINSTENCY( ap->array[i], az.array[i] ) ) info=1;
    if( VECTOR_COPY( ap->array[i], az.array[i] ) ) info=1;
  }
  

  return info;

}

int     __my_vector_array_compare( const vector_array az1, const vector_array az2 ){

  /* dummies */
  int dim, i;
  int info=0;


  if( POINTER_CHECK( az1.array, __my_vector_array_copy ) ||
      ( az1.dim != az2.dim )
      ) info=1;

  dim = az1.dim;

  for( i=0; i<dim; i++ ){

    if( VECTOR_CONSINSTENCY( az1.array[i], az2.array[i] ) ) info=1;

    if( VECTOR_COMPARE( az1.array[i], az2.array[i] ) ){

      fprintf( stderr, "ERROR: the %d components are not equal\n", i );
      fflush( stderr );

      info=1;

      break;

    }

  }
  

  return info;

}

int __my_vector_array_read( FILE* fp, vector_array_p ap ){

  /* dummies */
  int dim, i;
  int info=0;


  if( FILE_CHECK( fp, __my_vector_array_copy ) ) info=1;

  if( fread( ( void* ) &dim, sizeof( dim ), 1, fp ) < 1 ) info=1;

  ap->dim = dim;

  for( i=0; i<dim; i++ ){

    if( VECTOR_READ( fp, ap->array[ i ] ) ) info=1;

  }


  return info;

}

int __my_vector_array_print( FILE* fp, vector_array az ){

  /* dummies */
  int dim, i;
  int info=0;


  if( FILE_CHECK( fp, __my_vector_array_copy ) || POINTER_CHECK( az.array, __my_vector_array_copy ) ) info=1;

  dim = az.dim;

  if( fwrite( ( const void* ) &dim, sizeof( dim ), 1, fp ) < 1 ) info=1;

  for( i=0; i<dim; i++){

    if( VECTOR_PRINT( fp, az.array[ i ] ) ) info=1;

  }


  return info;

}
int     __my_vector_array_print_plus( FILE* fp, vector_array az ){

 /* dummies */
  int dim, i;
  int info=0;

  if( FILE_CHECK( fp, __my_vector_array_copy ) ||
      POINTER_CHECK( az.array, __my_vector_array_copy )
      ) info=1;

  dim = az.dim;

  fprintf( fp, "# vector_array dim  = %d\n", dim );

  for( i=0; i<dim; i++){

  fprintf( fp, "# array index       = %d\n", i );

    if( VECTOR_PRINT_PLUS( fp, az.array[i] ) ) info=1;


  }


  return info;

}

int     __my_vector_array_zero( vector_array_p vec_array_p){

  /* dummies */
  int i, dim;
  int info=0;


  dim = vec_array_p->dim;

  for( i=0; i<dim; i++ ){

    if( !info ){

      if( VECTOR_ZERO( vec_array_p->array[ i ] ) ) info=1;

    }
    else{

      break;

    }

  }

  return info;

}

/* ----------------------------------------------- */

/* int vectors */

int __my_ivector_allocate( const int dim, ivector_p vec_p ){


  vec_p->ivector_dim  = dim;
  vec_p->ivector      = (int*) calloc( (size_t) dim, sizeof(int) );

  if( (IVECTOR_CHECK( *vec_p )) ) return 1;

  return 0;

}

int       __my_ivector_free( ivector_p vec_p){

  if( (IVECTOR_CHECK( *vec_p )) ) return 1;

  vec_p->ivector_dim = 0;

  free( vec_p->ivector );

  vec_p->ivector = NULL;

  return 0;

}

int       __my_ivector_copy( ivector_p vec1_p, const ivector vec2  ){

  int dim, i;


  if( IVECTOR_CONSINSTENCY( *vec1_p, vec2 ) ) return 1; 

  dim = vec2.ivector_dim;

  for( i=0; i<dim; i++ ){
    vec1_p->ivector[i] = vec2.ivector[i];
  }


  return 0;

}

int       __my_ivector_compare( const ivector vec1, const ivector vec2  ){

  /* dummies */
  int dim, i;
  int info=0;


  if( IVECTOR_CONSINSTENCY( vec1, vec2 ) ) info=1; 

  dim = vec2.ivector_dim;

  for( i=0; i<dim; i++ ){

    if( vec1.ivector[ i ] != vec2.ivector[ i ] ){

      fprintf( stderr, "ERROR: the %d components are not equal\n", i );
      fflush( stderr );

      info=1;

      break;

    }

  }


  return info;

}

int       __my_ivector_compare_mute( const ivector vec1, const ivector vec2  ){

  /* dummies */
  int dim, i;
  int info=0;


#ifdef __DEBUG__

  if( IVECTOR_CONSINSTENCY( vec1, vec2 ) ) info=1; 

#endif /* __DEBUG__ */


  dim = vec2.ivector_dim;

  for( i=0; i<dim; i++ ){

    if( vec1.ivector[ i ] != vec2.ivector[ i ] ){

      info=1;

      break;

    }

  }


  return info;

}

int __my_ivector_read( FILE* file_p, ivector_p vec_p ){

  /* dummies */
  int dim;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_ivector_print ) ) || ( IVECTOR_CHECK( *vec_p ) ) ) info=1;

  if( fread( ( void* ) &dim, sizeof( dim ), 1, file_p ) < 1 ) info=1;

  if( dim != vec_p->ivector_dim ) info=1;

  if( fread( ( void* ) &vec_p->ivector, sizeof( vec_p->ivector ), 1, file_p ) < 1 ) info=1;


  return info;

}

int  __my_ivector_print( FILE* file_p, const ivector vec ){

  /* dummies */
  int dim;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_ivector_print ) ) || ( IVECTOR_CHECK( vec ) ) ) info=1;

  dim = vec.ivector_dim;

  if( fwrite( ( const void* ) &dim, sizeof( dim ), 1, file_p ) < 1 ) info=1;

  if( fwrite( ( const void* ) &vec.ivector, sizeof( vec.ivector ), dim, file_p ) < 1 ) info=1;


  return info;

}

int       __my_ivector_print_plus( FILE* file_p, const ivector vec ){

  int i, dim;


  if( (FILE_CHECK( file_p, __my_ivector_print )) || (IVECTOR_CHECK( vec )) ) return 1;

  dim = vec.ivector_dim;

  fprintf( file_p, "# ivector dim = %d\n", dim );

  for( i=0; i<dim; i++){

    fprintf( file_p, "# %d\n", vec.ivector[i] );

  }
  fflush( file_p );

  return 0;

}

int       __my_ivector_check( const ivector vec ){

  if( !(vec.ivector_dim) || !(vec.ivector) ){
    fprintf( stderr, "IVECTOR CHECK ERROR.\n" );
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }

}

int       __my_ivector_consinstency( const ivector vec1, const ivector vec2 ){

  if( (IVECTOR_CHECK( vec1 )) || (IVECTOR_CHECK( vec2 )) || (vec1.ivector_dim != vec2.ivector_dim )){
    fprintf( stderr, "IVECTOR CONSISTENCY ERROR.\n" );
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }

}

int    __my_ivector_set_to_zero( ivector_p vec_p ){

  int i, dim;


  if( (IVECTOR_CHECK( *vec_p )) ) return 1;

  dim = vec_p->ivector_dim;

  for( i=0; i<dim; i++ ){
    vec_p->ivector[i] = 0;
  }

  return 0;

}

int    __my_ivector_set_to_unit( ivector_p vec_p ){

  int i, dim;


  if( (IVECTOR_CHECK( *vec_p )) ) return 1;

  dim = vec_p->ivector_dim;

  for( i=0; i<dim; i++ ){
    vec_p->ivector[i] = 1;
  }

  return 0;

}

int __my_ivector_are_equal( const ivector vec1, const ivector vec2 ){

  /* dummies */
  int i, dim;
  int info=0;

#ifdef __DEBUG__

  if( IVECTOR_CONSINSTENCY( vec1, vec2 ) ) info=1;

#endif /* __DEBUG__ */

  dim = vec1.ivector_dim;

  for( i=0; i<dim; i++ ){

    if( vec2.ivector[ i ] != vec1.ivector[ i ]){

      fprintf( stderr, "ERROR: entry %d is different\n", i ); 
      fflush( stderr );
      info =1;

      break;

    }

  }


  return info;

}



/* ----------------------------------------------- */

/* real vectors */

int __my_rvector_allocate( const int dim, rvector_p vec_p ){


  vec_p->rvector_dim  = dim;
  vec_p->rvector      = (double*) calloc( (size_t) dim, sizeof(double) );

  if( (RVECTOR_CHECK( *vec_p )) ) return 1;

  return 0;

}

int       __my_rvector_free( rvector_p vec_p){

  if( (RVECTOR_CHECK( *vec_p )) ) return 1;

  vec_p->rvector_dim = 0;

  free( vec_p->rvector );

  vec_p->rvector = NULL;

  return 0;

}

int       __my_rvector_copy( rvector_p vec1_p, const rvector vec2  ){

  int dim, i;


  if( RVECTOR_CONSINSTENCY( *vec1_p, vec2 ) ) return 1; 

  dim = vec2.rvector_dim;

  for( i=0; i<dim; i++ ){
    vec1_p->rvector[i] = vec2.rvector[i];
  }

  return 0;

}

int       __my_rvector_compare( const rvector vec1, const rvector vec2  ){

  /* dummies */
  int dim, i;
  int info=0;


  if( RVECTOR_CONSINSTENCY( vec1, vec2 ) ) info=1; 

  dim = vec2.rvector_dim;

  for( i=0; i<dim; i++ ){

    if( fabs( vec1.rvector[ i ] - vec2.rvector[ i ] ) > DBL_EPSILON ){

      fprintf( stderr, "ERROR: the %d components are not equal\n", i );
      fflush( stderr );

      info=1;

      break;

    }

  }


  return info;

}

int __my_rvector_read( FILE* file_p, rvector_p vec_p ){

  /* dummies */
  int dim;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_rvector_print ) ) || ( RVECTOR_CHECK( *vec_p ) ) ) info=1;

  if( fread( ( void* ) &dim, sizeof( dim ), 1, file_p ) < 1 ) info=1;

  if( dim != vec_p->rvector_dim ) info=1;
  
  if( fread( ( void* ) &vec_p->rvector, sizeof( vec_p->rvector ), 1, file_p ) < 1 ) info=1;


  return info;

}

int __my_rvector_print( FILE* file_p, const rvector vec ){

  /* dummies */
  int dim;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_rvector_print ) ) || ( RVECTOR_CHECK( vec ) ) ) info=1;

  dim = vec.rvector_dim;

  if( fwrite( ( const void* ) &dim, sizeof( dim ), 1, file_p ) < 1 ) info=1;

  if( fwrite( ( const void* ) &vec.rvector, sizeof( vec.rvector ), 1, file_p ) < 1 ) info=1;


  return info;

}

int __my_rvector_print_plus( FILE* file_p, const rvector vec ){

  int i, dim;


  if( (FILE_CHECK( file_p, __my_rvector_print )) || (RVECTOR_CHECK( vec )) ) return 1;

  dim = vec.rvector_dim;

  fprintf( file_p, "# rvector dim = %d\n", dim );

  for( i=0; i<dim; i++){

    fprintf( file_p, "# "DOUBLE_FORMAT"\n", vec.rvector[i] );

  }
  fflush( file_p );

  return 0;

}

int       __my_rvector_check( const rvector vec ){

  if( !(vec.rvector_dim) || !(vec.rvector) ){
    fprintf( stderr, "RVECTOR CHECK ERROR.\n" );
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }

}

int       __my_rvector_consinstency( const rvector vec1, const rvector vec2 ){

  if( (RVECTOR_CHECK( vec1 )) || (RVECTOR_CHECK( vec2 )) || (vec1.rvector_dim != vec2.rvector_dim )){
    fprintf( stderr, "RVECTOR CONSISTENCY ERROR.\n" );
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }

}

int    __my_rvector_set_to_zero( rvector_p vec_p ){

  int i, dim;


  if( (RVECTOR_CHECK( *vec_p )) ) return 1;

  dim = vec_p->rvector_dim;

  for( i=0; i<dim; i++ ){
    vec_p->rvector[i] = 0.0e0;
  }

  return 0;

}

int    __my_rvector_set_to_unit( rvector_p vec_p ){

  int i, dim;


  if( (RVECTOR_CHECK( *vec_p )) ) return 1;

  dim = vec_p->rvector_dim;

  for( i=0; i<dim; i++ ){
    vec_p->rvector[i] = 1.0e0;
  }

  return 0;

}
