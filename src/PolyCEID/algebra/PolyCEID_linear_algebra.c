
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
#include "PolyCEID_linear_algebra.h"



/*----------------------------------------- */

/* private working space utilities */ 

/*----------------------------------------- */

#define WORKING_DIM 1000000 /* static working space dimension, sizeof(double). PRIVATE, it must be here */

int     int_space[ WORKING_DIM ];
double  working_space1[ WORKING_DIM ];
double  working_space2[ WORKING_DIM ];
double  working_space3[ WORKING_DIM ];

/*----------------------------------------- */

int       __my_working_space_check( const int dim ){

  if( dim > WORKING_DIM ){
    fprintf( stderr, "STATIC WORKING SPACE EXCEEDED.\n");
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }
  
} 

/*----------------------------------------- */

int       __my_vector_transcription( const vector vec, double* array_p ){

  int     real_dim, i;
  double* real_vector_p;


#ifdef __DEBUG__
  if( (VECTOR_CHECK( vec )) ) return 1;
#endif /* __DEBUG__*/

  real_dim = 2* vec.vector_dim; /* complex = 2*double */

  real_vector_p = (double*)vec.vector;
  for( i=0; (i<real_dim  && i<WORKING_DIM); i++ ){
    array_p[i] = real_vector_p[i]; /* I bet it can be coded better */
  }

  return 0;

}

/*----------------------------------------- */

int       __my_matrix_transcription( const matrix mat, double* array_p ){

  int     real_dim, i;
  double* real_matrix_p;


#ifdef __DEBUG__
  if( MATRIX_CHECK( mat ) ) return 1;
#endif /* __DEBUG__*/

  real_dim = 2* mat.matrix_dim;  /* complex = 2*double */

#ifdef __DEBUG__
  if( WORKING_SPACE_CHECK( real_dim ) ) return 1;
#endif /* __DEBUG__*/

  real_matrix_p = (double*)mat.matrix;
  for( i=0; i<real_dim; i++ ){
    array_p[i] = real_matrix_p[i]; /* I bet it can be coded better */
  }

  /*
    Possibly, a more efficienty way is by the auxiliari Lapack subroutine ZLACPY
  */

  return 0;

}

/*----------------------------------------- */

int       __my_vector_pointer_check( const vector v, const vector w ){

  if( v.vector == w.vector ){
    fprintf( stderr, "VECTOR POINTERS SUPERPOSITION ERROR.\n" );
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }

}

int       __my_matrix_pointer_check( const matrix a, const matrix b ){

  if( a.matrix == b.matrix ){
    fprintf( stderr, "MATRIX POINTERS SUPERPOSITION ERROR.\n" );
    fflush( stderr );
    return 1;
  } else{
    return 0;
  }

}

/*----------------------------------------- */

/* utilities */

int       __my_matrix_vector_product_consinstency( const matrix mat, const vector vec1, const vector vec2 ){

  if( (MATRIX_CHECK( mat )) || (VECTOR_CHECK( vec1 )) || (VECTOR_CHECK( vec2 )) || (mat.matrix_dim_column != vec1.vector_dim) 
      || (vec1.vector_dim != vec2.vector_dim) ){

    fprintf( stderr, "MATRIX_VECTOR_PRODUCT_CONSISTENCY ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;

}

/*----------------------------------------- */

int       __my_matrix_matrix_product_consinstency( const matrix mat1, const matrix mat2, const matrix mat3 ){

  if( (MATRIX_CHECK( mat1 )) || (MATRIX_CHECK( mat2 )) || (MATRIX_CHECK( mat3 )) || (mat1.matrix_dim_column != mat2.matrix_dim_row) 
      || (mat1.matrix_dim_row != mat3.matrix_dim_row) || (mat2.matrix_dim_column != mat3.matrix_dim_column) ){

    fprintf( stderr, "MATRIX_MATRIX_PRODUCT_CONSISTENCY ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;

}

/*----------------------------------------- */

/* blas & lapack */

/*----------------------------------------- */

double __my_complex_norm( const complex z ){

  int    n=1;
  int    incx=1;
  double norm;


  norm = F77_FUNC( dznrm2, DZNRM2 )( &n, (double*)(z.z), &incx );

  return norm;

}


complex __my_scalar_product( const vector a, const vector b ){

  complex z;
  int     dim;
  int     stride_a=1;
  int     stride_b=1;


#ifdef __DEBUG__
  if( VECTOR_CONSINSTENCY( a, b ) ) return CMPLX_HUGE_VAL;
#endif /* __DEBUG__ */


  dim      = a.vector_dim;

  //  z = F77_FUNC(zdotc,ZDOTC)(&dim, (double*)&a.vector[0], &stride_a, (double*)&b.vector[0], &stride_b);
  z = F77_FUNC(zdotc,ZDOTC)( &dim, (double*)a.vector, &stride_a, (double*)b.vector, &stride_b );


  return z;

}

/*----------------------------------------- */

int     __my_matrix_vector_product( const matrix a, const vector v, vector_p w_p){


  char    trans ='N';
  int     m;
  int     n;
  complex alpha;
  int     lda;
  int     incx = 1;
  complex beta;
  int     incy = 1;

  
#ifdef __DEBUG__
  if( (MATRIX_VECTOR_PRODUCT_CONSINSTENCY( a, v, *w_p )) || (VECTOR_POINTER_CHECK( *w_p, v)) ) return 1;
#endif /* __DEBUG__ */

  m     = a.matrix_dim_row;
  n     = a.matrix_dim_column;
  alpha = CMPLX_UNIT;
  beta  = CMPLX_ZERO;
  lda   = MAX( 1, m );

  F77_FUNC(zgemv,ZGEMV)( &trans, &m, &n, alpha.z, (double*)a.matrix, &lda, (double*)v.vector, &incx, beta.z, (double*)(w_p->vector), &incy );
  
  return 0;
  
}

/*----------------------------------------- */

int     __my_matrix_matrix_product( const matrix a, const matrix b, matrix_p c_p ){

  char    transa ='N';
  char    transb ='N';
  int     m;
  int     n;
  int     k;
  complex alpha;
  int     lda;
  int     ldb;
  complex beta;
  int     ldc;


#ifdef __DEBUG__
  if( (MATRIX_MATRIX_PRODUCT_CONSINSTENCY( a, b, *w_p )) || (MATRIX_POINTER_CHECK( *c_p, a )) 
      || (MATRIX_POINTER_CHECK( *c_p, b )) ) return 1;
#endif /* __DEBUG__ */

  m     = a.matrix_dim_row; 
  n     = b.matrix_dim_column; 
  k     = a.matrix_dim_column; 
  alpha = CMPLX_UNIT;
  beta  = CMPLX_ZERO;
  lda   = MAX( 1, m ); /* if transa = 'N' */
  ldb   = MAX( 1, k ); /* if transk = 'N' */
  ldc   = MAX( 1, m );

  F77_FUNC(zgemm,ZGEMM)( &transa, &transb, &m, &n, &k, alpha.z, (double*)a.matrix, &lda, (double*)b.matrix, &ldb, beta.z, (double*)(c_p->matrix), &ldc );

  return 0;

}

/*----------------------------------------- */

double  __my_vector_norm( const vector v ){

  double norm;
  int    dim;
  int    incx=1;


#ifdef __DEBUG__
  if( VECTOR_CHECK( v ) ) return HUGE_VAL;
#endif /* _DEBUG__ */

  dim  = v.vector_dim;

  norm =  F77_FUNC( dznrm2, DZNRM2 )( &dim, (double*)v.vector, &incx );

  return norm;

}

/*----------------------------------------- */

double  __my_rvector_norm( const rvector v ){

  double norm;
  int    dim;
  int    incx=1;


#ifdef __DEBUG__
  if( RVECTOR_CHECK( v ) ) return HUGE_VAL;
#endif /* _DEBUG__ */

  dim  = v.rvector_dim;

  norm =  F77_FUNC( dnrm2, DNRM2 )( &dim, v.rvector, &incx );

  return norm;

}

/*----------------------------------------- */

double  __my_matrix_norm( const matrix a ){

  double  fr_norm ;
  char    norm ='F';
  int     m, n, lda;
  //double* work_p;   /* working space, dimension ldwork >= m*/


  m      = a.matrix_dim_row;
#ifdef __DEBUG__
  if( ( MATRIX_CHECK( a ) ) || ( WORKING_SPACE_CHECK( m ) )  ) return HUGE_VAL;
#endif /* __DEBUG__ */

  n      = a.matrix_dim_column; 
  lda    = MAX( m, 1 );
  //work_p =  working_space1; /* using static workin space */

  fr_norm =  F77_FUNC( zlange, ZLANGE )( &norm, &m, &n, (double*)a.matrix, &lda, working_space1 );

  return fr_norm;

}

/*----------------------------------------- */

int     __my_diagonalisation( const matrix a, matrix_p vr_p, rvector_p w_p ){

  char jobz ='V';
  char uplo ='U';
  int  n, lda, lwork;
  int  info =0;


#ifdef __DEBUG__
  if( (MATRIX_SQUARE_CHECK( a )) ||  (MATRIX_SQUARE_CHECK( *vr_p )) ||  
      (RVECTOR_CHECK( *w_p )) || ( a.matrix_dim_column != w_p->rvector_dim ) || 
      ( a.matrix_dim_column != vr_p->matrix_dim_column ) ) return 1;
#endif /* __DEBUG__ */

  n     = a.matrix_dim_row;
  lwork = MAX( 1, 2*n-1 );

#ifdef __DEBUG__
  if( (WORKING_SPACE_CHECK( 2*lwork )) || (WORKING_SPACE_CHECK( 2*n )) ) return 1;
#endif /* __DEBUG__ */

  if( MATRIX_COPY( *vr_p, a ) ) return 1; /* This is compulsary, since matrix a will be overwritten otherwise */

  lda   = MAX( 1, n );
  lwork = WORKING_DIM /2; /* That's for better efficiency. notice that the previous definition allows for check. */

  F77_FUNC( zheev, ZHEEV )( &jobz, &uplo, &n, (double*)(vr_p->matrix), &lda, (double*)w_p->rvector, 
			    working_space1, &lwork, working_space2, &info);

  return 0;

}

/*----------------------------------------- */

int     __my_non_symmetric_diagonalisation( const matrix a, matrix_p vl_p, matrix_p vr_p, vector_p w_p ){

  char jobvl ='N';
  char jobvr ='N';
  int  n, lda, ldvl, ldvr, lwork;
  int  info =0;

#ifdef __DEBUG__
  if( (MATRIX_SQUARE_CHECK( a )) ||  (MATRIX_SQUARE_CHECK( *vl_p )) ||  (MATRIX_SQUARE_CHECK( *vr_p )) ||  
      (VECTOR_CHECK( *w_p )) || ( a.matrix_dim_column != w_p->vector_dim ) || 
      ( a.matrix_dim_column != vl_p->matrix_dim_column ) || ( a.matrix_dim_column != vr_p->matrix_dim_column ) ) return 1;
#endif /* __DEBUG__ */

  n     = a.matrix_dim_row;
  lwork = MAX( 1, 2*n ); /* here only for checking purposes */
 
#ifdef __DEBUG__
  if( (WORKING_SPACE_CHECK( 2*lwork )) || (WORKING_SPACE_CHECK( 2*n )) ) return 1;
#endif /* __DEBUG__ */

  if( (MATRIX_TRANSCRIPTION( a, working_space1 )) ) return 1;
  
  lda   = MAX( 1, n );
  lwork = WORKING_DIM /2; /* That's for better efficiency. */
  ldvl  = n;
  ldvr  = n;

  F77_FUNC( zgeev, ZGEEV )( &jobvl, &jobvr, &n, working_space1, &lda, (double*)(w_p->vector), (double*)(vl_p->matrix),  
			    &ldvl, (double*)(vr_p->matrix), &ldvr, working_space2, &lwork, working_space3, &info );

  if( info != 0 ) return info;

  return 0;

}

/*----------------------------------------- */

int     __my_matrix_inverse( const matrix mat1 , matrix_p mat_p ){

  int      i, dim;
  int      n, lda, lwork;
  int      info=0;
  //  complex* dummy_pointer;


#ifdef __DEBUG__
  if( (MATRIX_SQUARE_CHECK( mat1 )) ||  (MATRIX_SQUARE_CHECK( *mat_p )) ||   
      ( mat1.matrix_dim_row != mat_p->matrix_dim_row ) ) return 1; 
#endif /* __DEBUG__ */

  n     = mat1.matrix_dim_row;

#ifdef __DEBUG__
  if( WORKING_SPACE_CHECK( n ) ) return 1;
#endif /* __DEBUG__ */

  if( (MATRIX_TRANSCRIPTION( mat1, working_space1 )) ) return 1;

  lda   = MAX( 1, n );

  F77_FUNC( zgetrf, ZGETRF )( &n, &n, working_space1, &lda, int_space, &info );

  if( info != 0 ){
    fprintf( stderr, "ERROR: ZGETRF returned non-zero value.\n");
    fflush( stderr );
    return info;
  }

  lwork = MAX( 1, n ); /* here only for checking purposes */

#ifdef __DEBUG__
  if( WORKING_SPACE_CHECK( 2*lwork ) ) return 1;
#endif /* __DEBUG__ */

  /* lda already defined */
  lwork = WORKING_DIM /2; /* That's for better efficiency. */

  F77_FUNC( zgetri, ZGETRI )( &n, working_space1, &lda, int_space, working_space2, &lwork, &info);
 
  if( info != 0 ){
    fprintf( stderr, "ERROR: ZGETRI returned non-zero value.\n");
    fflush( stderr );
    return info;
  }


  /* final copy */

  dim = mat_p->matrix_dim;

  /*
  dummy_pointer = (complex*)working_space1;

  the line above produces the following gcc warning:
  warning: dereferencing type-punned pointer will break strict-aliasing rules

  for( i=0; i<dim; i++ ){
    mat_p->matrix[ i ] = dummy_pointer[ i ];
  }
  */
  
  for( i=0; i<dim; i++ ){

    mat_p->matrix[ i ].z[ 0 ] = working_space1[ 2 *i    ];
    mat_p->matrix[ i ].z[ 1 ] = working_space1[ 2 *i +1 ];
    
  }


  return 0;

}

/*----------------------------------------- */

int     __my_rvector_increase( rvector_p vec1_p, const rvector vec2, double alpha ){

  /* dummies */
  int  dim;
  int  incx=1;
  int  incy=1;
  int  info=0;


#ifdef __DEBUG__

  if( RVECTOR_CONSINSTENCY( *vec1_p, vec2 ) ){

    info=1;

  }

#endif /* __DEBUG__ */

  dim = vec1_p->rvector_dim;

  F77_FUNC( daxpy, DAXPY )( &dim, &alpha, vec2.rvector, &incx, vec1_p->rvector, &incy );


  return info;

}

/*----------------------------------------- */

/* non_blas & non-lapack functions */

/*----------------------------------------- */

int     __my_matrix_adjoint( const matrix mat1, matrix_p mat_p ){

  int     i, j, dim;


#ifdef __DEBUG__
  if( (MATRIX_SQUARE_CHECK( mat1 )) ||  (MATRIX_SQUARE_CHECK( *mat_p )) ||   
      ( mat1.matrix_dim_row != mat_p->matrix_dim_row ) ) return 1; 
#endif /* __DEBUG__ */

  dim = mat1.matrix_dim_row;

  for( i=0; i<dim; i++ ){
    for( j=0; j<dim; j++ ){

      mat_p->matrix[ MATRIX_ARG( i, j, dim, dim ) ] = CONJ( mat1.matrix[ MATRIX_ARG( j, i, dim, dim ) ] );

    }
  }

  return 0;

}

/*----------------------------------------- */

int     __my_matrix_conjugate( const matrix mat1, matrix_p mat_p ){

  int     i, j, dim;


#ifdef __DEBUG__
  if( (MATRIX_SQUARE_CHECK( mat1 )) ||  (MATRIX_SQUARE_CHECK( *mat_p )) ||   
      ( mat1.matrix_dim_row != mat_p->matrix_dim_row ) ) return 1; 
#endif /* __DEBUG__ */

  dim = mat1.matrix_dim_row;

  for( i=0; i<dim; i++ ){
    for( j=0; j<dim; j++ ){

      mat_p->matrix[ MATRIX_ARG( i, j, dim, dim ) ] = CONJ( mat1.matrix[ MATRIX_ARG( i, j, dim, dim ) ] );

    }
  }

  return 0;

}

/*----------------------------------------- */

int     __my_matrix_transpose( const matrix mat1, matrix_p mat_p ){

  int     i, j, dim;


#ifdef __DEBUG__
  if( (MATRIX_SQUARE_CHECK( mat1 )) ||  (MATRIX_SQUARE_CHECK( *mat_p )) ||   
      ( mat1.matrix_dim_row != mat_p->matrix_dim_row ) ) return 1; 
#endif /* __DEBUG__ */

  dim = mat1.matrix_dim_row;

  for( i=0; i<dim; i++ ){
    for( j=0; j<dim; j++ ){

      mat_p->matrix[ MATRIX_ARG( i, j, dim, dim ) ] = mat1.matrix[ MATRIX_ARG( j, i, dim, dim ) ];

    }
  }

  return 0;

}

/*----------------------------------------- */

complex __my_matrix_trace( const matrix mat ){

  int     i, dim;
  complex trace;


#ifdef __DEBUG__
  if( MATRIX_SQUARE_CHECK( mat ) ){
    fprintf( stderr, "ERROR in  __my_matrix_trace: this is not a square matrix!\n" );
    fflush( stderr );
    return CMPLX_HUGE_VAL;
  }
#endif /* _DEBUG__ */

  dim = mat.matrix_dim_row;

  trace = CMPLX_ZERO; /* set the trace to zero */
  for( i=0; i<dim; i++ ){

    trace.z[0] += mat.matrix[ MATRIX_ARG( i, i, dim, dim) ].z[0];
    trace.z[1] += mat.matrix[ MATRIX_ARG( i, i, dim, dim) ].z[1];

  }

  return trace;

}

/*----------------------------------------- */

int     __my_matrix_scalar_product( const matrix mat1, const complex scalar, matrix_p mat_p ){

  /* dummies */
  int i, dim;
  int info=0;

#ifdef __DEBUG__
  if( MATRIX_CONSINSTENCY( mat1, *mat_p ) ) info=1;
#endif /* __DEBUG__ */

  dim = mat_p->matrix_dim;

  for( i=0; i<dim; i++ ){

    mat_p->matrix[ i ] = CMPLX_PRODUCT( mat1.matrix[ i ], scalar );

  }

  return info;

}

/*----------------------------------------- */

int     __my_matrix_scalar_division( const matrix mat1, const complex scalar, matrix_p mat_p ){

  /* dummies */
  int i, dim;
  int info=0;


#ifdef __DEBUG__
  if( MATRIX_CONSINSTENCY( mat1, *mat_p ) ) info=1;
#endif /* __DEBUG__ */

  dim = mat_p->matrix_dim;


  for( i=0; i<dim; i++ ){

    mat_p->matrix[ i ] = CMPLX_DIVISION( mat1.matrix[ i ], scalar );

  }

  return info;

}

/*----------------------------------------- */

int     __my_matrix_addition( const matrix mat1, const matrix mat2, matrix_p mat_p ){

  int i, dim;


#ifdef __DEBUG__
  if( MATRIX_CONSINSTENCY( mat1, *mat_p ) ||  MATRIX_CONSINSTENCY( mat2, *mat_p ) ) return 1;
#endif /* __DEBUG__ */

  dim = mat_p->matrix_dim;

  for( i=0; i<dim; i++ ){
    mat_p->matrix[ i ] = CMPLX_SUM( mat1.matrix[ i ], mat2.matrix[ i ] );
  }

  return 0;

}

/*----------------------------------------- */

int     __my_matrix_subtraction( const matrix mat1, const matrix mat2, matrix_p mat_p ){

  int i, dim;


#ifdef __DEBUG__
  if( MATRIX_CONSINSTENCY( mat1, *mat_p ) ||  MATRIX_CONSINSTENCY( mat2, *mat_p ) ) return 1;
#endif /* __DEBUG__ */

  dim = mat_p->matrix_dim;

  for( i=0; i<dim; i++ ){
    mat_p->matrix[ i ] = CMPLX_DIF( mat1.matrix[ i ], mat2.matrix[ i ] );
  }

  return 0;

}

/*----------------------------------------- */

int __my_matrix_commutator( const matrix a, const matrix b, matrix_p c_p ){

  /*
    WARNING: not optimized!
  */

  //  static int counter=1;
  int     i, j, k, dim;
  int     index0, index1, index2;
  complex dummy1, dummy2, dummy3;


#ifdef __DEBUG__
  if( (MATRIX_MATRIX_PRODUCT_CONSINSTENCY( a, b, *c_p )) || (a.matrix_dim_column != c_p->matrix_dim_row) || 
      (b.matrix_dim_row != c_p->matrix_dim_column) ) return 1;
#endif /* __DEBUG__ */

  dim = c_p->matrix_dim_row;

  /*
  fprintf( stdout, "# --------------\n");
  fprintf( stdout, "# Commutator of:\n");

  MATRIX_PRINT( stdout, a );

  fprintf( stdout, "# and \n");

  MATRIX_PRINT( stdout, b );

  fprintf( stdout, "# gives  \n");
  */

  for( i=0; i<dim; i++){
    for( j=0; j<dim; j++){

      c_p->matrix[ MATRIX_ARG( i, j, dim, dim) ] = CMPLX_ZERO;
      
      index0 = MATRIX_ARG( i, j, dim, dim );

      for( k=0; k<dim; k++){
	
	index1 = MATRIX_ARG( i, k, dim, dim );

	index2 = MATRIX_ARG( k, j, dim, dim );

	dummy1 = CMPLX_PRODUCT( a.matrix[ index1 ], b.matrix[ index2 ] );
	dummy2 = CMPLX_PRODUCT( b.matrix[ index1 ], a.matrix[ index2 ] );

	dummy3 = CMPLX_DIF( dummy1, dummy2 );

	c_p->matrix[ index0 ].z[0] += dummy3.z[0];
	c_p->matrix[ index0 ].z[1] += dummy3.z[1];
	
      }
    }
  }

  /*
  MATRIX_PRINT( stdout, *c_p );
  fprintf( stdout, "# --------------\n");
  */

  /*
    fprintf( stdout, "# __my_matrix_commutator called %d times \n", counter );
    counter++;
  */


  return 0;

}

/*----------------------------------------- */

int __my_matrix_commutator_alt( const matrix a, const matrix b, matrix_p c_p ){

  /* dummies */
  char    transa ='N';
  char    transb ='N';
  int     m;
  int     n;
  int     k;
  complex alpha;
  int     lda;
  int     ldb;
  complex beta;
  int     ldc;
  int     i;
  int     info=0;


#ifdef __DEBUG__

  if( (MATRIX_MATRIX_PRODUCT_CONSINSTENCY( a, b, *c_p )) || (a.matrix_dim_column != c_p->matrix_dim_row) || 
      (b.matrix_dim_row != c_p->matrix_dim_column) ){

    info=1;

  }

#endif /* __DEBUG__ */


  m     = a.matrix_dim_row; 
  n     = b.matrix_dim_column; 
  k     = a.matrix_dim_column; 
  alpha = CMPLX_UNIT;
  beta  = CMPLX_ZERO;
  lda   = MAX( 1, m ); /* if transa = 'N' */
  ldb   = MAX( 1, k ); /* if transk = 'N' */
  ldc   = MAX( 1, m );


  /* a times b */
  F77_FUNC(zgemm,ZGEMM)( &transa, &transb, &m, &n, &k, alpha.z, (double*)a.matrix, &lda, (double*)b.matrix, &ldb, beta.z, working_space1, &ldc );

  /* b times a */
  F77_FUNC(zgemm,ZGEMM)( &transa, &transb, &m, &n, &k, alpha.z, (double*)b.matrix, &lda, (double*)a.matrix, &ldb, beta.z, working_space2, &ldc );


  /* final difference */
  for( i=0; i<( k *k ); i++ ){

    c_p->matrix[ i ].z[ 0 ] = working_space1[ 2 *i    ] -working_space2[ 2 *i    ];

    c_p->matrix[ i ].z[ 1 ] = working_space1[ 2 *i +1 ] -working_space2[ 2 *i +1 ];

  }


  return info;

}

/*----------------------------------------- */

int __my_matrix_anticommutator( const matrix a, const matrix b, matrix_p c_p ){

  /*
    WARNING: not optimized!
  */

  int     i, j, k, dim;
  int     index0, index1, index2;
  complex dummy1, dummy2, dummy3;


#ifdef __DEBUG__
  if( (MATRIX_MATRIX_PRODUCT_CONSINSTENCY( a, b, *c_p )) || (a.matrix_dim_column != c_p->matrix_dim_row) || 
      (b.matrix_dim_row != c_p->matrix_dim_column) ) return 1;
#endif /* __DEBUG__ */

  dim = c_p->matrix_dim_row;

  /*
  fprintf( stdout, "# --------------\n");
  fprintf( stdout, "# Anticommutator of:\n");

  MATRIX_PRINT( stdout, a );

  fprintf( stdout, "# and \n");

  MATRIX_PRINT( stdout, b );

  fprintf( stdout, "# gives  \n");
  */

  for( i=0; i<dim; i++){
    for( j=0; j<dim; j++){

      c_p->matrix[ MATRIX_ARG( i, j, dim, dim) ] = CMPLX_ZERO;
      
      index0 = MATRIX_ARG( i, j, dim, dim );

      for( k=0; k<dim; k++){
	
	index1 = MATRIX_ARG( i, k, dim, dim );

	index2 = MATRIX_ARG( k, j, dim, dim );

	dummy1 = CMPLX_PRODUCT( a.matrix[ index1 ], b.matrix[ index2 ] );
	dummy2 = CMPLX_PRODUCT( b.matrix[ index1 ], a.matrix[ index2 ] );

	dummy3 = CMPLX_SUM( dummy1, dummy2 );

	c_p->matrix[ index0 ].z[0] += dummy3.z[0];
	c_p->matrix[ index0 ].z[1] += dummy3.z[1];
	
      }
    }
  }

  /*
  MATRIX_PRINT( stdout, *c_p );
  fprintf( stdout, "# --------------\n");
  */

  return 0;

}

/*----------------------------------------- */

int __my_matrix_anticommutator_alt( const matrix a, const matrix b, matrix_p c_p ){

  /* dummies */
  char    transa ='N';
  char    transb ='N';
  int     m;
  int     n;
  int     k;
  complex alpha;
  int     lda;
  int     ldb;
  complex beta;
  int     ldc;
  int     i;
  int     info=0;


#ifdef __DEBUG__

  if( (MATRIX_MATRIX_PRODUCT_CONSINSTENCY( a, b, *c_p )) || (a.matrix_dim_column != c_p->matrix_dim_row) || 
      (b.matrix_dim_row != c_p->matrix_dim_column) ){

    info=1;

  }

#endif /* __DEBUG__ */


  m     = a.matrix_dim_row; 
  n     = b.matrix_dim_column; 
  k     = a.matrix_dim_column; 
  alpha = CMPLX_UNIT;
  beta  = CMPLX_ZERO;
  lda   = MAX( 1, m ); /* if transa = 'N' */
  ldb   = MAX( 1, k ); /* if transk = 'N' */
  ldc   = MAX( 1, m );


  /* a times b */
  F77_FUNC(zgemm,ZGEMM)( &transa, &transb, &m, &n, &k, alpha.z, (double*)a.matrix, &lda, (double*)b.matrix, &ldb, beta.z, working_space1, &ldc );

  /* b times a */
  F77_FUNC(zgemm,ZGEMM)( &transa, &transb, &m, &n, &k, alpha.z, (double*)b.matrix, &lda, (double*)a.matrix, &ldb, beta.z, working_space2, &ldc );


  /* final difference */
  for( i=0; i<( k *k ); i++ ){

    c_p->matrix[ i ].z[ 0 ] = working_space1[ 2 *i    ] +working_space2[ 2 *i    ];

    c_p->matrix[ i ].z[ 1 ] = working_space1[ 2 *i +1 ] +working_space2[ 2 *i +1 ];

  }


  return info;

}

/*----------------------------------------- */

int __my_matrix_commutator_hermitian( const matrix a, const matrix b, matrix_p c_p ){

  /* dummies */
  int     i, j, dim;
  int     index1, index2;
  int     info=0;


#ifdef __DEBUG__

  if( (MATRIX_MATRIX_PRODUCT_CONSINSTENCY( a, b, *c_p )) || (a.matrix_dim_column != c_p->matrix_dim_row) || 
      (b.matrix_dim_row != c_p->matrix_dim_column) ) info=1;

#endif /* __DEBUG__ */

  dim = c_p->matrix_dim_row;

  /*
  fprintf( stdout, "# --------------\n");
  fprintf( stdout, "# Commutator of:\n");

  MATRIX_PRINT_PLUS( stdout, a );

  fprintf( stdout, "# and \n");

  MATRIX_PRINT_PLUS( stdout, b );

  fprintf( stdout, "# gives  \n");
  */

  if( MATRIX_MATRIX_PRODUCT( a, b, *c_p ) ) info=1;

  for( i=0; i<dim; i++){

    for( j=i+1; j<dim; j++){

      index1=MATRIX_ARG( i, j, dim, dim );

      index2=MATRIX_ARG( j, i, dim, dim );


      c_p->matrix[ index1 ].z[0] -= c_p->matrix[ index2 ].z[0];
      
      c_p->matrix[ index1 ].z[1] += c_p->matrix[ index2 ].z[1];


      c_p->matrix[ index2 ].z[0]  = -c_p->matrix[ index1 ].z[0];

      c_p->matrix[ index2 ].z[1]  =  c_p->matrix[ index1 ].z[1];
	
    }

    index1=MATRIX_ARG( i, i, dim, dim );

    c_p->matrix[ index1 ].z[0]  = 0.0e0;

    c_p->matrix[ index1 ].z[1] *= 2.0e0;

  }

  /*
  MATRIX_PRINT( stdout, *c_p );
  fprintf( stdout, "# --------------\n");
  */

  return info;

}

/*----------------------------------------- */

int __my_matrix_anticommutator_hermitian( const matrix a, const matrix b, matrix_p c_p ){

  /* dummies */
  int     i, j, dim;
  int     index1, index2;
  int     info=0;


#ifdef __DEBUG__

  if( (MATRIX_MATRIX_PRODUCT_CONSINSTENCY( a, b, *c_p )) || (a.matrix_dim_column != c_p->matrix_dim_row) || 
      (b.matrix_dim_row != c_p->matrix_dim_column) ) info=1;

#endif /* __DEBUG__ */

  dim = c_p->matrix_dim_row;

  /*
  fprintf( stdout, "# --------------\n");
  fprintf( stdout, "# Anticommutator of:\n");

  MATRIX_PRINT_PLUS( stdout, a );

  fprintf( stdout, "# and \n");

  MATRIX_PRINT_PLUS( stdout, b );

  fprintf( stdout, "# gives  \n");
  */

  if( MATRIX_MATRIX_PRODUCT( a, b, *c_p )) info=1;

  for( i=0; i<dim; i++){

    for( j=i+1; j<dim; j++){

      index1=MATRIX_ARG( i, j, dim, dim );

      index2=MATRIX_ARG( j, i, dim, dim );


      c_p->matrix[ index1 ].z[0] += c_p->matrix[ index2 ].z[0];
      
      c_p->matrix[ index1 ].z[1] -= c_p->matrix[ index2 ].z[1];


      c_p->matrix[ index2 ].z[0]  =  c_p->matrix[ index1 ].z[0];

      c_p->matrix[ index2 ].z[1]  = -c_p->matrix[ index1 ].z[1];
	
    }

    index1=MATRIX_ARG( i, i, dim, dim );

    c_p->matrix[ index1 ].z[0] *= 2.0e0;

    c_p->matrix[ index1 ].z[1]  = 0.0e0;

  }

  /*
  MATRIX_PRINT( stdout, *c_p );
  fprintf( stdout, "# --------------\n");
  */

  return info;

}

/*----------------------------------------- */
