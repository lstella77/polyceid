
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
#include "PolyCEID_matrix.h"


int __my_matrix_allocate( const int dim_row, const int dim_column, matrix_p mat_p ){


  int         dim, i_column;
  complex_pp  column_p;


  mat_p->matrix_dim_row     = dim_row;
  mat_p->matrix_dim_column  = dim_column;

  dim                = dim_row *dim_column;
  mat_p->matrix_dim  = dim;
  mat_p->matrix      = (complex*) calloc( (size_t) dim, sizeof(complex) );

  /* allocate a vector for every column */
  mat_p->column_p = (complex**)  calloc( (size_t) dim_column, sizeof(complex*) );

  if( (POINTER_CHECK( mat_p->column_p, __my_matrix_allocate  )) ) return 1;

  column_p =  mat_p->column_p;
  for( i_column=0; i_column<dim_column; i_column++, column_p++ ){

    *column_p = &mat_p->matrix[ i_column *dim_column ]; /* column beginning */
    if( (POINTER_CHECK( *column_p, __my_matrix_allocate  )) ) return 1;

  }

  if( (MATRIX_CHECK( *mat_p )) ) return 1;

  return 0;

}

int       __my_matrix_free( matrix_p mat_p){

  complex_pp  column_p;

  if( (MATRIX_CHECK( *mat_p )) ) return 1;

  mat_p->matrix_dim_row    = 0;
  mat_p->matrix_dim_column = 0;
  mat_p->matrix_dim        = 0;

  free( mat_p->matrix );

  mat_p->matrix = NULL;

  /* deallocate columns */
  column_p   = mat_p->column_p;

  if( (POINTER_CHECK( column_p, __my_matrix_free  )) ) return 1;
  free( column_p );

  return 0;

}

int       __my_matrix_copy( matrix_p mat1_p, const matrix mat2 ){

  int i, dim;

  if( MATRIX_CONSINSTENCY( *mat1_p, mat2 ) ) return 1;

  dim = mat2.matrix_dim;

  for( i=0; i<dim; i++){
    //CMPLX_COPY( mat1_p->matrix[i], mat2.matrix[i] );
    mat1_p->matrix[i] = mat2.matrix[i];
  }

  return 0; 

}

int       __my_matrix_compare( const matrix mat1, const matrix mat2 ){

  /* dummies */
  int i, dim;
  int info=0;


  if( MATRIX_CONSINSTENCY( mat1, mat2 ) ) info=1;

  dim = mat2.matrix_dim;

  for( i=0; i<dim; i++){

    if( CMPLX_COMPARE( mat1.matrix[ i ], mat2.matrix[ i ] ) ){

      fprintf( stderr, "ERROR: the %d components are not equal\n", i );
      fflush( stderr );

      info=1;

      break;

    } 

  }


  return info; 

}


int __my_matrix_read( FILE* file_p, matrix_p mat_p ){

  /* dummies */
  int i;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_matrix_read ) ) || ( MATRIX_CHECK( *mat_p ) ) ) info=1;

  if( fread( ( void* ) &mat_p->matrix_dim_row, sizeof( int ), 1, file_p ) < 1 ) info=1;

  if( fread( ( void* ) &mat_p->matrix_dim_column, sizeof( int ), 1, file_p ) < 1 ) info=1;

  for( i=0; i<mat_p->matrix_dim; i++ ){

    CMPLX_READ( file_p, mat_p->matrix[ i ] );

  }

  
  return info;

}

int __my_matrix_print( FILE* file_p, const matrix mat ){

  /* dummies */
  int i;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_matrix_print ) ) || ( MATRIX_CHECK( mat ) ) ) info=1;

  if( fwrite( ( const void* ) &mat.matrix_dim_row, sizeof( int ), 1, file_p ) < 1 ) info=1;

  if( fwrite( ( const void* ) &mat.matrix_dim_column, sizeof( int ), 1, file_p ) < 1 ) info=1;

  for( i=0; i<mat.matrix_dim; i++ ){

    CMPLX_PRINT( file_p, mat.matrix[ i ] );
  
  }

  
  return info;

}

int __my_matrix_print_plus( FILE* file_p, const matrix mat){

  int dim_row, dim_column, i_row, i_column;


  if( (FILE_CHECK( file_p, __my_matrix_print  )) || (MATRIX_CHECK( mat )) ) return 1;

  dim_row    = mat.matrix_dim_row;
  dim_column = mat.matrix_dim_column;

  fprintf( file_p, "# matrix dim_row    = %d\n", dim_row    );
  fprintf( file_p, "# matrix dim_column = %d\n", dim_column );

  for( i_row=0; i_row<dim_row; i_row++){

    fprintf( file_p, "# " );

    for( i_column=0; i_column<dim_column; i_column++){
      CMPLX_PRINT_PLUS_VAR( file_p, mat.matrix[ MATRIX_ARG( i_row, i_column, dim_row, dim_column ) ] );
      fprintf( file_p, " ");
    }
    fprintf( file_p, "\n");
  }
  fflush( file_p );

  return 0;

}

/* array */

int     __my_matrix_array_allocate( int dim, int dim_row, int dim_column, matrix_array_p ap ){

  /* dummies */
  int i;
  int info=0;


  ap->dim = dim;

  ap->array = (matrix_p) calloc( (size_t)(dim), sizeof(matrix) );

  if( POINTER_CHECK( ap->array, __my_matrix_array_allocate ) ) info=1;

  if( !info ){
    
    for( i=0; i<dim; i++ ){

      if( MATRIX_ALLOCATE( dim_row, dim_column, ap->array[i] )) info=1;

    }

  }


  return info;

}

int     __my_matrix_array_free( matrix_array_p ap ){

  /* dummies */
  int i, dim;
  int info=0;


  if( POINTER_CHECK( ap->array, __my_matrix_array_free ) ) info=1;

  if( !info ){

    dim = ap->dim;

    for( i=0; i<dim; i++ ){

      if( MATRIX_FREE( ap->array[i] ) ) info=1;

    }
  
    ap->dim = 0;

    free( ap->array );

    ap->array = NULL;

  }


  return info;

}

int     __my_matrix_array_copy( matrix_array_p ap, const matrix_array az ){

  /* dummies */
  int dim, i;
  int info=0;

  if( POINTER_CHECK( ap->array, __my_matrix_array_copy ) ||
      ( ap->dim != az.dim )
      ) info=1;

  dim = ap->dim;

  for( i=0; i<dim; i++){
    if( MATRIX_CONSINSTENCY( ap->array[i], az.array[i] ) ) info=1;
    if( MATRIX_COPY( ap->array[i], az.array[i] ) ) info=1;
  }
  

  return info;

}

int     __my_matrix_array_compare( const matrix_array az1, const matrix_array az2 ){

  /* dummies */
  int dim, i;
  int info=0;


  if( POINTER_CHECK( az1.array, __my_matrix_array_copy ) ||
      ( az1.dim != az2.dim )
      ) info=1;

  dim = az1.dim;

  for( i=0; i<dim; i++ ){

    if( MATRIX_CONSINSTENCY( az1.array[i], az2.array[i] ) ) info=1;

    if( MATRIX_COMPARE( az1.array[i], az2.array[i] ) ){

      fprintf( stderr, "ERROR: the %d components are not equal\n", i );
      fflush( stderr );

      info=1;

      break;

    }

  }
  

  return info;

}

int  __my_matrix_array_read( FILE* fp, matrix_array_p ap ){

  /* dummies */
  int i;
  int info=0;


  if( FILE_CHECK( fp, __my_matrix_array_copy ) ) info=1;

  if( fread( ( void* ) &ap->dim, sizeof( int ), 1, fp ) < 1 ) info=1;

  for( i=0; i<ap->dim; i++ ){

    if( MATRIX_READ( fp, ap->array[ i ] ) ) info=1;

  }


  return info;

}

int __my_matrix_array_print( FILE* fp, matrix_array az ){

  /* dummies */
  int i;
  int info=0;


  if( FILE_CHECK( fp, __my_matrix_array_copy ) || POINTER_CHECK( az.array, __my_matrix_array_copy ) ) info=1;

  if( fwrite( ( const void* ) &az.dim, sizeof( int ), 1, fp ) < 1 ) info=1;

  for( i=0; i<az.dim; i++ ){

    if( MATRIX_PRINT( fp, az.array[ i ] ) ) info=1;

  }


  return info;

}

int __my_matrix_array_print_plus( FILE* fp, matrix_array az ){

 /* dummies */
  int dim, i;
  int info=0;

  if( FILE_CHECK( fp, __my_matrix_array_copy ) ||
      POINTER_CHECK( az.array, __my_matrix_array_copy )
      ) info=1;

  dim = az.dim;

  fprintf( fp, "# matrix_array dim  = %d\n", dim );

  for( i=0; i<dim; i++){

  fprintf( fp, "# array index       = %d\n", i );

    if( MATRIX_PRINT_PLUS( fp, az.array[i] ) ) info=1;


  }


  return info;

}


/* utilities */

int __my_matrix_arg( const int i_row, const int i_column, const int dim_row, const int dim_column ){

  /* dummies */
  int index =-1; //WARNING: default values means an error occurred
  //  int info=0;


  if( i_row<0 || i_row>dim_row-1 || i_column<0 || i_column>dim_column-1 ){ //WARNING: no further checks are required

    fprintf( stderr, "ERROR: illegal matrix indeces.\n");
    fprintf( stderr, "i_row = %d, dim_row = %d, i_column = %d, d_column = %d \n", i_row, dim_row, i_column, dim_column );
    fflush( stderr );

#ifdef __DEBUG__ 

    fprintf( stdout, "ERROR: illegal matrix indeces.\n");
    fprintf( stdout, "i_row = %d, dim_row = %d, i_column = %d, d_column = %d \n", i_row, dim_row, i_column, dim_column );
    fflush( stdout );

#endif /* __DEBUG__ */


    //    info=1;

  }
  else{

    index = dim_row *i_column +i_row;

  }


  return index;

} 

int __my_matrix_arg_inverse( const int index, int* i_row_p, int* i_column_p, const int dim_row, const int dim_column ){

  /* dummies */
  int dim;
  int info=0;

  
  dim = dim_row *dim_column;

  if( index < 0 || index > dim-1 ){ //WARNING: no further checks are required

    fprintf( stderr, "ERROR: illegal scalar index.\n");

    info=1;

  }
  else{

    *i_row_p    = index %dim_row;

    *i_column_p = index /dim_row; //WARNING: as an integer operation

  }


  return info;

} 

int       __my_matrix_check( const matrix mat ){

  if( !(mat.matrix_dim_row) || !(mat.matrix_dim_column) || !(mat.matrix) || !(mat.column_p) ){
    fprintf( stderr, "MATRIX CHECK ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;
  
}

int       __my_matrix_consinstency( const matrix mat1, const matrix mat2 ){

  if( (MATRIX_CHECK( mat1 )) || (MATRIX_CHECK( mat2 )) || (mat1.matrix_dim_row != mat2.matrix_dim_row) || 
      (mat1.matrix_dim_column != mat2.matrix_dim_column) ){

    fprintf( stderr, "MATRIX CONSISTENCY ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;

}

int       __my_matrix_square_check( const matrix mat ){

  if( (MATRIX_CHECK( mat )) || ( mat.matrix_dim_row != mat.matrix_dim_column ) ){
    fprintf( stderr, "MATRIX SQUARE CHECK ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;
  
}

int    __my_matrix_set_to_zero( matrix_p mat_p ){

  int i, dim;


  if( (MATRIX_CHECK( *mat_p)) ) return 1;

  dim = mat_p->matrix_dim;

  for( i=0; i<dim; i++ ){
    //CMPLX_COPY( mat_p->matrix[i], CMPLX_ZERO );
    //mat_p->matrix[i] = CMPLX_ZERO; /* WARNING! That could be dangerous! */
    CMPLX_ZERO_VAR( mat_p->matrix[i] );
  }

  return 0;

}

int    __my_matrix_set_to_unit( matrix_p mat_p ){

  int i_row, i_column, dim_row, dim_column;


  if( (MATRIX_CHECK( *mat_p )) ) return 1;

  dim_row    = mat_p->matrix_dim_row;
  dim_column = mat_p->matrix_dim_column;

  if( dim_column != dim_row  ) return 1; /* only squared matrix */

  for( i_row=0; i_row<dim_row; i_row++ ){
    for( i_column=0; i_column<dim_column; i_column++ ){

      if( i_column != i_row  ){
	//CMPLX_COPY( mat_p->matrix[ MATRIX_ARG( i_row, i_column, dim_row, dim_column ) ], CMPLX_ZERO );
	//	mat_p->matrix[ MATRIX_ARG( i_row, i_column, dim_row, dim_column ) ] =  CMPLX_ZERO; /* WARNING! That could be dangerous! */
	CMPLX_ZERO_VAR( mat_p->matrix[ MATRIX_ARG( i_row, i_column, dim_row, dim_column ) ] );
      } else{
	//CMPLX_COPY( mat_p->matrix[ MATRIX_ARG( i_row, i_column, dim_row, dim_column ) ], CMPLX_UNIT );
	//	mat_p->matrix[ MATRIX_ARG( i_row, i_column, dim_row, dim_column ) ] =  CMPLX_UNIT; /* WARNING! That could be dangerous! */
	CMPLX_UNIT_VAR( mat_p->matrix[ MATRIX_ARG( i_row, i_column, dim_row, dim_column ) ] );
      }

    }
  }

  return 0;

}

int    __my_matrix_set_to_sigma_x( matrix_p mat_p ){

  if( (MATRIX_CHECK( *mat_p )) ) return 1;

  if( (mat_p->matrix_dim_row !=2 ) || (mat_p->matrix_dim_column != 2) ) return 1;

  mat_p->matrix[0] = CMPLX_ZERO;
  mat_p->matrix[1] = CMPLX_UNIT; 
  mat_p->matrix[2] = CMPLX_UNIT;
  mat_p->matrix[3] = CMPLX_ZERO; 
  /* WARNING! Those equalities  could be dangerous! CMPLX_COPY may be used instead */

  return 0;

}

int    __my_matrix_set_to_sigma_y( matrix_p mat_p ){

  if( (MATRIX_CHECK( *mat_p )) ) return 1;

  if( (mat_p->matrix_dim_row !=2 ) || (mat_p->matrix_dim_column != 2) ) return 1;


  mat_p->matrix[0] = CMPLX_ZERO;
  mat_p->matrix[1] = CMPLX_I; 
  mat_p->matrix[2] = CMPLX_INIT( 0.0e0, -1.0e0 ) ;
  mat_p->matrix[3] = CMPLX_ZERO; 
  /* WARNING! Those equalities  could be dangerous! CMPLX_COPY may be used instead */

  return 0;

}

int    __my_matrix_set_to_sigma_z( matrix_p mat_p ){

  if( (MATRIX_CHECK( *mat_p )) ) return 1;

  if( (mat_p->matrix_dim_row !=2 ) || (mat_p->matrix_dim_column != 2) ) return 1;

  mat_p->matrix[0] = CMPLX_UNIT;
  mat_p->matrix[1] = CMPLX_ZERO; 
  mat_p->matrix[2] = CMPLX_ZERO;
  mat_p->matrix[3] = CMPLX_INIT( -1.0e0, 0.0e0 ); 
  /* WARNING! Those equalities  could be dangerous! CMPLX_COPY may be used instead */

  return 0;

}

int     __my_matrix_array_zero( matrix_array_p mat_array_p){

  /* dummies */
  int i, dim;
  int info=0;


  dim = mat_array_p->dim;

  for( i=0; i<dim; i++ ){

    if( !info ){

      if( MATRIX_ZERO( mat_array_p->array[ i ] ) ) info=1;

    }
    else{

      break;

    }

  }

  return info;

}


/* imatrix */

int       __my_imatrix_check( const imatrix imat ){

  if( !(imat.imatrix_dim_row) || !(imat.imatrix_dim_column) || !(imat.imatrix) ){
    fprintf( stderr, "IMATRIX CHECK ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;
  
}

int       __my_imatrix_consinstency( const imatrix imat1, const imatrix imat2 ){

  if( (IMATRIX_CHECK( imat1 )) || (IMATRIX_CHECK( imat2 )) || (imat1.imatrix_dim_row != imat2.imatrix_dim_row) || 
      (imat1.imatrix_dim_column != imat2.imatrix_dim_column) ){

    fprintf( stderr, "MATRIX CONSISTENCY ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;

}

int __my_imatrix_allocate( const int dim_row, const int dim_column, imatrix_p imat_p ){


  int         i;


  imat_p->imatrix_dim_row     = dim_row;
  imat_p->imatrix_dim_column  = dim_column;

  imat_p->imatrix      = (int**) calloc( (size_t) imat_p->imatrix_dim_row, sizeof(int*) );

  if( (POINTER_CHECK( imat_p->imatrix, __my_imatrix_allocate  ) ) ) return 1;

  for( i=0; i<imat_p->imatrix_dim_row; i++ ){

    imat_p->imatrix[ i ] = (int*) calloc( (size_t) imat_p->imatrix_dim_column, sizeof(int) );

    if( (POINTER_CHECK( imat_p->imatrix[ i ], __my_imatrix_allocate  ) ) ){

      return 1;

      break;

    }

  }

  if( (IMATRIX_CHECK( *imat_p )) ) return 1;


  return 0;

}

int       __my_imatrix_free( imatrix_p imat_p){

  int i;

  if( (IMATRIX_CHECK( *imat_p )) ) return 1;

  for( i=0; i<imat_p->imatrix_dim_row; i++ ){

    free( imat_p->imatrix[ i ] );

  }

  free( imat_p->imatrix );


  imat_p->imatrix_dim_column =0;
  imat_p->imatrix_dim_row    =0;


  return 0;

}

int       __my_imatrix_copy( imatrix_p imat1_p, const imatrix imat2 ){

  int i, j;

  if( IMATRIX_CONSINSTENCY( *imat1_p, imat2 ) ) return 1;

  for( i=0; i<imat1_p->imatrix_dim_row ; i++){

    for( j=0; j<imat1_p->imatrix_dim_column ; j++){

      imat1_p->imatrix[i][j] = imat2.imatrix[i][j];

    }

  }

  return 0; 

}

int       __my_imatrix_compare( const imatrix imat1, const imatrix imat2 ){

  /* dummies */
  int i, j;
  int info=0;


  if( IMATRIX_CONSINSTENCY( imat1, imat2 ) ) info=1;


  for( i=0; i<imat1.imatrix_dim_row ; i++){

    for( j=0; j<imat1.imatrix_dim_column ; j++){

      if( imat1.imatrix[ i ][ j ] != imat2.imatrix[ i ][ j ] ){

	fprintf( stderr, "ERROR: the %d components are not equal\n", i );
	fflush( stderr );
	
	info=1;
	
	break;
	
      } 
      
    }

  }


  return info; 

}


int __my_imatrix_read( FILE* file_p, imatrix_p imat_p ){

  /* dummies */
  int i, j;
  int info=0;


  if( (FILE_CHECK( file_p, __my_imatrix_print ) ) || (IMATRIX_CHECK( *imat_p ) ) ) info=1;

  if( fread( ( void* ) &imat_p->imatrix_dim_row, sizeof( int ), 1, file_p ) < 1 ) info=1;

  if( fread( ( void* ) &imat_p->imatrix_dim_column, sizeof( int ), 1, file_p ) < 1 ) info=1;

  for( i=0; i<imat_p->imatrix_dim_row; i++ ){

    for( j=0; j<imat_p->imatrix_dim_column; j++ ){

      if( fread( ( void* ) &imat_p->imatrix[i][j], sizeof( int ), 1, file_p ) < 1 ) info=1;
  
    }

  }  


  return info;

}

int __my_imatrix_print( FILE* file_p, const imatrix imat){

  /* dummies */
  int i, j;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_imatrix_print ) ) || ( IMATRIX_CHECK( imat ) ) ) info=1;

  if( fwrite( ( const void* ) &imat.imatrix_dim_row, sizeof( int ), 1, file_p ) < 1 ) info=1;

  if( fwrite( ( const void* ) &imat.imatrix_dim_column, sizeof( int ), 1, file_p ) < 1 ) info=1;
  
  for( i=0; i<imat.imatrix_dim_row; i++ ){

    for( j=0; j<imat.imatrix_dim_column; j++ ){

      if( fwrite( ( const void* ) &imat.imatrix[ i ][ j ], sizeof( int ), 1, file_p ) < 1 ) info=1;
  
    }

  }  


  return info;

}

int __my_imatrix_print_plus( FILE* file_p, const imatrix imat){

  int dim_row, dim_column, i, j;


  if( (FILE_CHECK( file_p, __my_imatrix_print  )) || (IMATRIX_CHECK( imat )) ) return 1;

  dim_row    = imat.imatrix_dim_row;
  dim_column = imat.imatrix_dim_column;

  fprintf( file_p, "# imatrix dim_row    = %d\n", dim_row    );
  fprintf( file_p, "# imatrix dim_column = %d\n", dim_column );


  for( i=0; i<imat.imatrix_dim_row; i++){

    fprintf( file_p, "# " );

    for( j=0; j<( imat.imatrix_dim_column -1 ); j++){

      fprintf( file_p, "%d ", imat.imatrix[ i ][ j ] );
      
    }

    fprintf( file_p, "%d \n", imat.imatrix[ i ][ j ] );

  }
  fflush( file_p );


  return 0;

}


int    __my_imatrix_set_to_zero( imatrix_p mat_p ){

  int i, j;
  int imatrix_dim_row;
  int imatrix_dim_column;


  if( (IMATRIX_CHECK( *mat_p ) ) ) return 1;

  imatrix_dim_row    = mat_p->imatrix_dim_row;
  imatrix_dim_column = mat_p->imatrix_dim_column;


  for( i=0; i<imatrix_dim_row; i++ ){

    for( j=0; j<imatrix_dim_column; j++ ){

      mat_p->imatrix[i][j] = 0; /* WARNING! That could be dangerous! */

    }

  }


  return 0;

}

/* rmatrix */

int       __my_rmatrix_check( const rmatrix mat ){

  if( !(mat.rmatrix_dim_row) || !(mat.rmatrix_dim_column) || !(mat.rmatrix) ){
    fprintf( stderr, "RMATRIX CHECK ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;
  
}

int       __my_rmatrix_consinstency( const rmatrix mat1, const rmatrix mat2 ){

  if( (RMATRIX_CHECK( mat1 )) || (RMATRIX_CHECK( mat2 )) || (mat1.rmatrix_dim_row != mat2.rmatrix_dim_row) || 
      (mat1.rmatrix_dim_column != mat2.rmatrix_dim_column) ){

    fprintf( stderr, "RMATRIX CONSISTENCY ERROR.\n" );
    fflush( stderr );
    return 1;
  } else
    return 0;

}

int __my_rmatrix_allocate( const int dim_row, const int dim_column, rmatrix_p mat_p ){


  int         dim, i_column;
  double**    rcolumn_p;


  mat_p->rmatrix_dim_row     = dim_row;
  mat_p->rmatrix_dim_column  = dim_column;

  dim                = dim_row *dim_column;
  mat_p->rmatrix_dim  = dim;
  mat_p->rmatrix      = (double*) calloc( (size_t) dim, sizeof(double) );

  /* allocate a vector for every column */
  mat_p->rcolumn_p = (double**)  calloc( (size_t) dim_column, sizeof(double*) );

  if( (POINTER_CHECK( mat_p->rcolumn_p, __my_rmatrix_allocate  )) ) return 1;

  rcolumn_p =  mat_p->rcolumn_p;
  for( i_column=0; i_column<dim_column; i_column++, rcolumn_p++ ){

    *rcolumn_p = &mat_p->rmatrix[ i_column *dim_column ]; /* column beginning */
    if( (POINTER_CHECK( *rcolumn_p, __my_rmatrix_allocate  )) ) return 1;

  }

  if( (RMATRIX_CHECK( *mat_p )) ) return 1;

  return 0;

}

int       __my_rmatrix_free( rmatrix_p mat_p){

  double**  rcolumn_p;

  if( (RMATRIX_CHECK( *mat_p )) ) return 1;

  mat_p->rmatrix_dim_row    = 0;
  mat_p->rmatrix_dim_column = 0;
  mat_p->rmatrix_dim        = 0;

  free( mat_p->rmatrix );

  mat_p->rmatrix = NULL;

  /* deallocate columns */
  rcolumn_p   = mat_p->rcolumn_p;

  if( (POINTER_CHECK( rcolumn_p, __my_rmatrix_free  )) ) return 1;
  free( rcolumn_p );

  return 0;

}

int       __my_rmatrix_copy( rmatrix_p mat1_p, const rmatrix mat2 ){

  int i, dim;

  if( RMATRIX_CONSINSTENCY( *mat1_p, mat2 ) ) return 1;

  dim = mat2.rmatrix_dim;

  for( i=0; i<dim; i++){
    mat1_p->rmatrix[i] = mat2.rmatrix[i];
  }

  return 0; 

}

int       __my_rmatrix_compare( const rmatrix mat1, const rmatrix mat2 ){

  /* dummies */
  int i, dim;
  int info=0;


  if( RMATRIX_CONSINSTENCY( mat1, mat2 ) ) info=1;

  dim = mat2.rmatrix_dim;

  for( i=0; i<dim; i++){

    if( fabs( mat1.rmatrix[ i ] -mat2.rmatrix[ i ] ) < DBL_EPSILON ){

      fprintf( stderr, "ERROR: the %d components are not equal\n", i );
      fflush( stderr );

      info=1;

      break;

    } 

  }


  return info; 

}


int __my_rmatrix_read( FILE* file_p, rmatrix_p mat_p ){

  /* dummies */
  int i;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_rmatrix_print ) ) || ( RMATRIX_CHECK( *mat_p ) ) ) info=1;

  if( fread( ( void* ) &mat_p->rmatrix_dim_row, sizeof( int ), 1, file_p ) < 1 ) info=1;

  if( fread( ( void* ) &mat_p->rmatrix_dim_column, sizeof( int ), 1, file_p ) < 1 ) info=1;

  for( i=0; i<( mat_p->rmatrix_dim_row *mat_p->rmatrix_dim_column ); i++ ){

    if( fread( ( void* ) &mat_p->rmatrix[i], sizeof( double ), 1, file_p ) < 1 ) info=1;

  }

  
  return info;

}

int __my_rmatrix_print( FILE* file_p, const rmatrix mat){

  /* dummies */
  int i;
  int info=0;


  if( ( FILE_CHECK( file_p, __my_rmatrix_print ) ) || ( RMATRIX_CHECK( mat ) ) ) info=1;

  if( fwrite( ( const void* ) &mat.rmatrix_dim_row, sizeof( int ), 1, file_p ) < 1 ) info=1;

  if( fwrite( ( const void* ) &mat.rmatrix_dim_column, sizeof( int ), 1, file_p ) < 1 ) info=1;

  for( i=0; i<( mat.rmatrix_dim_row *mat.rmatrix_dim_column ); i++ ){

    if( fwrite( ( const void* ) &mat.rmatrix[ i ], sizeof( double ), 1, file_p ) < 1 ) info=1;

  }


  return info;

}

int __my_rmatrix_print_plus( FILE* file_p, const rmatrix mat){

  int dim_row, dim_column, i_row, i_column;


  if( (FILE_CHECK( file_p, __my_rmatrix_print  )) || (RMATRIX_CHECK( mat )) ) return 1;

  dim_row    = mat.rmatrix_dim_row;
  dim_column = mat.rmatrix_dim_column;

  fprintf( file_p, "# rmatrix dim_row    = %d\n", dim_row    );
  fprintf( file_p, "# rmatrix dim_column = %d\n", dim_column );

  for( i_row=0; i_row<dim_row; i_row++){

    fprintf( file_p, "# " );

    for( i_column=0; i_column<dim_column; i_column++){
      fprintf( file_p, DOUBLE_FORMAT, mat.rmatrix[ MATRIX_ARG( i_row, i_column, dim_row, dim_column ) ] );
      fprintf( file_p, "  ");
    }
    fprintf( file_p, "\n");
  }
  fflush( file_p );

  return 0;

}

int    __my_rmatrix_set_to_zero( rmatrix_p mat_p ){

  int i, dim;


  if( (RMATRIX_CHECK( *mat_p)) ) return 1;

  dim = mat_p->rmatrix_dim;

  for( i=0; i<dim; i++ ){
    mat_p->rmatrix[i] = 0.0e0; /* WARNING! That could be dangerous! */
  }

  return 0;

}
