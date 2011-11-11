
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
#include "PolyCEID_indexing.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* private variables */

#ifdef __DEBUG__

unsigned short int indexing_allocate_flag=0;

#endif /* __DEBUG__ */


int  N_levels_single_loc;
int  N_levels_many_loc;
int  N_atoms_loc;
int  N_coor_loc;
int  N_coor_red_loc;
int  max_rho_index_loc;
int  sqrt_max_rho_index_loc;
int  CEID_order_loc;
int  rho_index_border_length_loc;
int  rho_index_next_to_border_length_loc;


typedef struct{

  imatrix rho_index_table;                       // [sqrt_max_rho_index, N_coor]

  imatrix rho_index_table_change1_plus1;         // [sqrt_max_rho_index, N_coor]

  imatrix rho_index_table_change1_minus1;        // [sqrt_max_rho_index, N_coor]

  imatrix rho_index_table_change1_plus2;         // [sqrt_max_rho_index, N_coor]

  imatrix rho_index_table_change1_minus2;        // [sqrt_max_rho_index, N_coor]

  imatrix rho_index_table_change2_plus1_plus1;   // [sqrt_max_rho_index, N_coor *N_coor]

  imatrix rho_index_table_change2_plus1_minus1;  // [sqrt_max_rho_index, N_coor *N_coor]

  imatrix rho_index_table_change2_minus1_minus1; // [sqrt_max_rho_index, N_coor *N_coor]


} rho_index_tables;

rho_index_tables tables;


//------------------------------------------
//                  COORDINATES
//------------------------------------------

/* coordinate_index */

int coordinate_index( int i, int j ){

  /* dummies */
  int index=-1; //WARNING: if error occurs, an invalid index will be returned
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    //    index = MATRIX_ARG( i, j, N_coor_loc, N_coor_loc );
    index = coordinate_index_table.imatrix[i][j];

  }


  return index;

}

//------------------------------------------

/* reduced_coordinate_index */

int reduced_coordinate_index( int i, int j ){

  /* dummies */
  int index=-1; //WARNING: if error occurs, an invalid index will be returned
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    //    index = MATRIX_ARG( i, j, N_coor_red_loc, N_coor_red_loc );
    index = reduced_coordinate_index_table.imatrix[i][j];

  }


  return index;

}


//------------------------------------------
//                  ATOMS
//------------------------------------------

/* atom_index */

int atom_index( int i, int j ){

  /* dummies */
  int index=-1; //WARNING: if error occurs, an invalid index will be returned
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    //    index = MATRIX_ARG( i, j, N_atoms_loc, N_atoms_loc );
    index = atom_index_table.imatrix[i][j];

  }


  return index;

}

//------------------------------------------
//                 ELECTRONS
//------------------------------------------

/* electron_many_index */

int electron_many_index( int i, int j ){

  /* dummies */
  int index=-1; //WARNING: if error occurs, an invalid index will be returned
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    //    index = MATRIX_ARG( i, j, N_levels_many_loc, N_levels_many_loc );
    index = electron_many_index_table.imatrix[i][j];

  }


  return index;

}

//------------------------------------------

/* electron_single_index */

int electron_single_index( int i, int j ){

  /* dummies */
  int index=-1; //WARNING: if error occurs, an invalid index will be returned
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    //    index = MATRIX_ARG( i, j, N_levels_single_loc, N_levels_single_loc );
    index = electron_single_index_table.imatrix[i][j];

  }


  return index;

}

//------------------------------------------
//                  RHO
//------------------------------------------

/* rho_index */

/*

WARNING: very time-consuming!

2st step: from two ivectors of dimension N_coor_red
to a single scalar index

*/

int rho_index( ivector index1, ivector index2 ){

  /* dummies */
  int dummy_superindex1;
  int dummy_superindex2;
  int index=-1;
  int info=0;


#ifdef __DEBUG__

  if( !info ){

    if( IVECTOR_CONSINSTENCY( index1, index2 ) ){

      fprintf( stderr, "ERROR: indexing inconsistency.\n");
      fflush( stderr );

      info=1;
      
    }

  }

#endif /* __DEBUG__ */


  if( !info ){

    dummy_superindex1 = SEARCH_RHO_INDEX_TABLE( index1 ); // WARNING: time-consuming operation

    dummy_superindex2 = SEARCH_RHO_INDEX_TABLE( index2 ); // WARNING: time-consuming operation

    if( dummy_superindex1 < 0 || dummy_superindex1 > ( sqrt_max_rho_index_loc -1 ) || 
	dummy_superindex2 < 0 || dummy_superindex2 > ( sqrt_max_rho_index_loc -1 ) ){ // WARNING: it's necessary

      index = max_rho_index_loc;

    }
    else{

      index = RHO_INDEX_AUX( dummy_superindex1, dummy_superindex2 );

    }

  } /* end info conditional */


  return index;


}

//------------------------------------------

/* rho_index_aux */

int rho_index_aux( int index1, int index2 ){

  /* dummies */
  int index=-1;
  int info=0;


  if( !info ){

    //index = MATRIX_ARG( index1, index2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
    index = rho_index_aux_table.imatrix[index1][index2];
 
  }


  return index;


}

//------------------------------------------
//------------------------------------------

/* rho_index_change1_plus1_aux */

int rho_index_change1_plus1_aux( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_index_loc-1 ){

      index_changed = sqrt_max_rho_index_loc;

    }
    else{

      index_changed = tables.rho_index_table_change1_plus1.imatrix[ index][ i ];

    }

  }

  return index_changed;

}

//------------------------------------------

/* rho_index_change1_minus1_aux */

int rho_index_change1_minus1_aux( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_index_loc-1 ){

      index_changed = sqrt_max_rho_index_loc;

    }
    else{

      index_changed = tables.rho_index_table_change1_minus1.imatrix[ index][ i ];
    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change1_plus2_aux */

int rho_index_change1_plus2_aux( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_index_loc-1 ){

      index_changed = sqrt_max_rho_index_loc;

    }
    else{
      
      index_changed = tables.rho_index_table_change1_plus2.imatrix[ index][ i ];
    
    }

  }


  return index_changed;

}

//------------------------------------------

/*  rho_index_change1_minus2_aux */

int rho_index_change1_minus2_aux( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_index_loc-1 ){
      
      index_changed = sqrt_max_rho_index_loc;
      
    }
    else{

      index_changed = tables.rho_index_table_change1_minus2.imatrix[ index][ i ];

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_plus1_plus1_aux */

int rho_index_change2_plus1_plus1_aux( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

  if( j < 0 || j > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: j value is not legal [%d]\n", j );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_index_loc-1 ){
      
      index_changed = sqrt_max_rho_index_loc;
      
    }
    else{

      index_changed = tables.rho_index_table_change2_plus1_plus1.imatrix[ index ][ REDUCED_COORDINATE_INDEX( i, j) ];

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_plus1_minus1_aux */

int rho_index_change2_plus1_minus1_aux( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault 
  int info=0;

#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

  if( j < 0 || j > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: j value is not legal [%d]\n", j );

    info=1;

  }


#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_index_loc-1 ){

      index_changed = sqrt_max_rho_index_loc;

    }
    else{

      index_changed = tables.rho_index_table_change2_plus1_minus1.imatrix[ index ][ REDUCED_COORDINATE_INDEX( i, j) ];

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_minus1_minus1_aux */

int rho_index_change2_minus1_minus1_aux( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

  if( j < 0 || j > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: j value is not legal [%d]\n", j );

    info=1;

  }


#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_index_loc-1 ){

      index_changed = sqrt_max_rho_index_loc;

    }
    else{

      index_changed = tables.rho_index_table_change2_minus1_minus1.imatrix[ index ][ REDUCED_COORDINATE_INDEX( i, j) ];

    }

  }


  return index_changed;

}

//------------------------------------------
//------------------------------------------

/* rho_index_change1_plus1_first */

int rho_index_change1_plus1_first( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_INDEX_CHANGE1_PLUS1_AUX( dummy_superindex1, i );      
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change1_minus1_first */

int rho_index_change1_minus1_first( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_INDEX_CHANGE1_MINUS1_AUX( dummy_superindex1, i );

	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_index_loc-1 ){
	  
	  index_changed = max_rho_index_loc;
	  
	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change1_plus2_first */

int rho_index_change1_plus2_first( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_INDEX_CHANGE1_PLUS2_AUX( dummy_superindex1, i );
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;


}

//------------------------------------------

/*  rho_index_change1_minus2_first */

int rho_index_change1_minus2_first( int index, int i ){


  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;
      
    }
    else{


      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_INDEX_CHANGE1_MINUS2_AUX( dummy_superindex1, i );
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{
	  
	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_plus1_plus1_first */

int rho_index_change2_plus1_plus1_first( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_INDEX_CHANGE2_PLUS1_PLUS1_AUX( dummy_superindex1, i, j );
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_plus1_minus1_first */

int rho_index_change2_plus1_minus1_first( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_INDEX_CHANGE2_PLUS1_MINUS1_AUX( dummy_superindex1, i, j );
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_minus1_minus1_first */

int rho_index_change2_minus1_minus1_first( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_INDEX_CHANGE2_MINUS1_MINUS1_AUX( dummy_superindex1, i, j );

	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{
	
	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------
//------------------------------------------

/* rho_index_change1_plus1_second */

int rho_index_change1_plus1_second( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_INDEX_CHANGE1_PLUS1_AUX( dummy_superindex2, i );      
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change1_minus1_second */

int rho_index_change1_minus1_second( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_INDEX_CHANGE1_MINUS1_AUX( dummy_superindex2, i );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change1_plus2_second */

int rho_index_change1_plus2_second( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_INDEX_CHANGE1_PLUS2_AUX( dummy_superindex2, i );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{
	  
	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;


}

//------------------------------------------

/*  rho_index_change1_minus2_second */

int rho_index_change1_minus2_second( int index, int i ){


  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;
      
    }
    else{


      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_INDEX_CHANGE1_MINUS2_AUX( dummy_superindex2, i );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_plus1_plus1_second */

int rho_index_change2_plus1_plus1_second( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_INDEX_CHANGE2_PLUS1_PLUS1_AUX( dummy_superindex2, i, j );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_plus1_minus1_second */

int rho_index_change2_plus1_minus1_second( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_INDEX_CHANGE2_PLUS1_MINUS1_AUX( dummy_superindex2, i, j );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_index_change2_minus1_minus1_second */

int rho_index_change2_minus1_minus1_second( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_index_loc-1 ){

      index_changed = max_rho_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_INDEX_CHANGE2_MINUS1_MINUS1_AUX( dummy_superindex2, i, j );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_index_loc-1 ){

	  index_changed = max_rho_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );
	  index_changed = rho_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* utilities */

//------------------------------------------

//------------------------------------------
//                  GENERAL
//------------------------------------------

/* initialise indexing  */

int indexing_allocate( constants constants, unsigned short int flag ){

  /* dummies */
  int i,j;
  int index, index_aux;
  int info=0;


#ifdef __DEBUG__

  if( indexing_allocate_flag ){

    info=1;

    fprintf( stderr, "ERROR: indexing allocation already done.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    N_levels_single_loc                 = constants.N_levels_single;
    N_levels_many_loc                   = constants.N_levels_many;
    N_atoms_loc                         = constants.N_atoms;
    N_coor_loc                          = constants.N_coor;
    N_coor_red_loc                      = constants.N_coor_red;
    CEID_order_loc                      = constants.CEID_order;
    max_rho_index_loc                   = constants.max_rho_index;
    sqrt_max_rho_index_loc              = constants.sqrt_max_rho_index;
    rho_index_border_length_loc         = constants.rho_index_border_length;
    rho_index_next_to_border_length_loc = constants.rho_index_next_to_border_length;


    if( !flag ){

      /* rho_index_aux_table */
      if( IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, sqrt_max_rho_index_loc, rho_index_aux_table ) ) info=1;
      
      for( i=0; i<sqrt_max_rho_index_loc; i++ ){
	for( j=0; j<sqrt_max_rho_index_loc; j++ ){
	  
	  rho_index_aux_table.imatrix[i][j] = MATRIX_ARG( i, j, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc );

	}
      }
      
      /*
	fprintf( stdout, "rho_index_aux_table\n" );
	if(IMATRIX_PRINT_PLUS( stdout, rho_index_aux_table ) ) info=1;
	fprintf( stdout, "\n" );
      */
      

      /* coordinate_index_table */
      if( IMATRIX_ALLOCATE( N_coor_loc, N_coor_loc, coordinate_index_table ) ) info=1;
      
      for( i=0; i<N_coor_loc; i++ ){
	for( j=0; j<N_coor_loc; j++ ){

	  coordinate_index_table.imatrix[i][j] = MATRIX_ARG( i, j, N_coor_loc, N_coor_loc );

	}
      }
      
      /*
	fprintf( stdout, "coordinate_index_table\n" );
	if(IMATRIX_PRINT_PLUS( stdout, coordinate_index_table ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /* reduced_coordinate_index_table */
      if( IMATRIX_ALLOCATE( N_coor_red_loc, N_coor_red_loc, reduced_coordinate_index_table ) ) info=1;
      
      for( i=0; i<N_coor_red_loc; i++ ){
	for( j=0; j<N_coor_red_loc; j++ ){

	  reduced_coordinate_index_table.imatrix[i][j] = MATRIX_ARG( i, j, N_coor_red_loc, N_coor_red_loc );

	}
      }
      
      /*
	fprintf( stdout, "reduced_coordinate_index_table\n" );
	if(IMATRIX_PRINT_PLUS( stdout, reduced_coordinate_index_table ) ) info=1;
	fprintf( stdout, "\n" );
      */
      
      /* atoms_index_table */
      if( IMATRIX_ALLOCATE( N_atoms_loc, N_atoms_loc, atom_index_table ) ) info=1;
      
      for( i=0; i<N_atoms_loc; i++ ){
	for( j=0; j<N_atoms_loc; j++ ){
	  
	  atom_index_table.imatrix[i][j] = MATRIX_ARG( i, j, N_atoms_loc, N_atoms_loc );

	}
      }
      
      /*
	fprintf( stdout, "atom_index_table\n" );
	if(IMATRIX_PRINT_PLUS( stdout, atom_index_table ) ) info=1;
	fprintf( stdout, "\n" );
      */

      /* electron_many_index_table */
      if( IMATRIX_ALLOCATE( N_levels_many_loc, N_levels_many_loc, electron_many_index_table ) ) info=1;

      for( i=0; i<N_levels_many_loc; i++ ){
	for( j=0; j<N_levels_many_loc; j++ ){
	  
	  electron_many_index_table.imatrix[i][j] = MATRIX_ARG( i, j, N_levels_many_loc, N_levels_many_loc );

	}
      }
      
      /*
	fprintf( stdout, "electron_levels_many_index_table\n" );
	if(IMATRIX_PRINT_PLUS( stdout, electron_many_index_table ) ) info=1;
	fprintf( stdout, "\n" );
      */
      
    } /* end flag conditional */



    /* electron_single_index_table */
    if( IMATRIX_ALLOCATE( N_levels_single_loc, N_levels_single_loc, electron_single_index_table ) ) info=1;
      
    for( i=0; i<N_levels_single_loc; i++ ){
      for( j=0; j<N_levels_single_loc; j++ ){
	
	electron_single_index_table.imatrix[i][j] = MATRIX_ARG( i, j, N_levels_single_loc, N_levels_single_loc );
	
      }
    }

    /*
      fprintf( stdout, "electron_levels_single_index_table\n" );
      if(IMATRIX_PRINT_PLUS( stdout, electron_single_index_table ) ) info=1;
      fprintf( stdout, "\n" );
    */


    if( !flag ){

      /* dummy_rho_index1 */
      if( IVECTOR_ALLOCATE( N_coor_red_loc, dummy_rho_index1 ) ) info=1;
      
      /* dummy_rho_index2 */
      if( IVECTOR_ALLOCATE( N_coor_red_loc, dummy_rho_index2 ) ) info=1;
      
      /* rho_index_table */
      if(IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, N_coor_red_loc, tables.rho_index_table ) ) info=1;
      
      /* rho_index_table_change1_plus1 */
      if(IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, N_coor_red_loc, tables.rho_index_table_change1_plus1 ) ) info=1;

      /* rho_index_table_change1_minus1 */
      if(IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, N_coor_red_loc, tables.rho_index_table_change1_minus1 ) ) info=1;

      /* rho_index_table_change1_plus2 */
      if(IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, N_coor_red_loc, tables.rho_index_table_change1_plus2 ) ) info=1;

      /* rho_index_table_change1_minus2 */
      if(IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, N_coor_red_loc, tables.rho_index_table_change1_minus2 ) ) info=1;

      /* rho_index_table_change2_plus1_plus1 */
      if(IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, N_coor_red_loc *N_coor_red_loc, tables.rho_index_table_change2_plus1_plus1 ) ) info=1;

      /* rho_index_table_change2_plus1_minus1 */
      if(IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, N_coor_red_loc *N_coor_red_loc, tables.rho_index_table_change2_plus1_minus1 ) ) info=1;

      /* rho_index_table_change2_minus1_minus1 */
      if(IMATRIX_ALLOCATE( sqrt_max_rho_index_loc, N_coor_red_loc *N_coor_red_loc, tables.rho_index_table_change2_minus1_minus1 ) ) info=1;

      /* diagonal_rho_indeces */
      if( IVECTOR_ALLOCATE( sqrt_max_rho_index_loc, diagonal_rho_indices ) ) info=1;
    
      /* rho_index_conjugate */
      if( IVECTOR_ALLOCATE( max_rho_index_loc, rho_index_conjugate ) ) info=1;

      /* rho_index_border */
      if(IVECTOR_ALLOCATE( rho_index_border_length_loc, rho_index_border ) ) info=1;

      if( rho_index_next_to_border_length_loc > 0 ){

	/* rho_index_next_to_border */
	if(IVECTOR_ALLOCATE( rho_index_next_to_border_length_loc, rho_index_next_to_border ) ) info=1;

      }
     

#ifdef __DEBUG_PLUS__

      fprintf( stdout, "DOING: compute_rho_index_table\n" );

#endif /* __DEBUG_PLUS__ */

      if( !info ){

	if( COMPUTE_RHO_INDEX_TABLE( 1 ) ) info=1; //WARNING: it restarts from scratch

      }

      
#ifdef __DEBUG_PLUS__
    
      fprintf( stdout, "\n" );
      fprintf( stdout, "DOING: compute_rho_index_other_tables\n" );

#endif /* __DEBUG_PLUS__ */
      
    
      if( !info ){
	
	if( COMPUTE_RHO_INDEX_OTHER_TABLES() ) info=1;

      }


#ifdef __DEBUG_PLUS__

      fprintf( stdout, "\n" );
      fprintf( stdout, "DOING: diagonal_rho_indeces\n" );

#endif /* __DEBUG_PLUS__ */
      

      if( !info ){
	
	for( i=0; i<sqrt_max_rho_index_loc; i++ ){

	  diagonal_rho_indices.ivector[ i ] = i *( sqrt_max_rho_index_loc +1 );

	} 

      }

#ifdef __DEBUG_PLUS__

      fprintf( stdout, "\n" );
      fprintf( stdout, "DOING: rho_index_conjugate \n" );

#endif /* __DEBUG_PLUS__ */

    
      if( !info ){
      
	for( index=0; index<max_rho_index_loc; index++ ){

	  if( RHO_INDEX_INVERSE( index, dummy_rho_index1, dummy_rho_index2 ) ) info=1;
	
	  index_aux = RHO_INDEX( dummy_rho_index2, dummy_rho_index1 ); // Notice the inverted order!

	  rho_index_conjugate.ivector[ index ] = index_aux;

	} 

	
      }


#ifdef __DEBUG_PLUS__

      fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */

    }

  }    


#ifdef __DEBUG__

  if( !info ){

    indexing_allocate_flag=1;
    
  }
  else{

    indexing_allocate_flag=0;
      
  }

#endif /* __DEBUG__ */
  

  return info;

}

//------------------------------------------

/* free indexing  */

int indexing_free( unsigned short int flag ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;

    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){
    
    if( !flag ){
 
      if( rho_index_next_to_border_length_loc > 0 ){
	
	/* rho_index_next_to_border */
	if(IVECTOR_FREE( rho_index_next_to_border ) ) info=1;
	
      }
      
      /* rho_index_border */
      if(IVECTOR_FREE( rho_index_border ) ) info=1;
      
      /* rho_index_conjugate */
      if( IVECTOR_FREE( rho_index_conjugate ) ) info=1;
      
      /* diagonal_rho_indeces*/
      if( IVECTOR_FREE( diagonal_rho_indices ) ) info=1;
      
      /* rho_index_table_change2_minus1_minus1 */
      if(IMATRIX_FREE( tables.rho_index_table_change2_minus1_minus1 ) ) info=1;
    
      /* rho_index_table_change2_plus1_minus1 */
      if(IMATRIX_FREE( tables.rho_index_table_change2_plus1_minus1 ) ) info=1;
    
      /* rho_index_table_change2_plus1_plus1 */
      if(IMATRIX_FREE( tables.rho_index_table_change2_plus1_plus1 ) ) info=1;
    
      /* rho_index_table_change1_minus2 */
      if(IMATRIX_FREE( tables.rho_index_table_change1_minus2 ) ) info=1;
      
      /* rho_index_table_change1_plus2 */
      if(IMATRIX_FREE( tables.rho_index_table_change1_plus2 ) ) info=1;
    
      /* rho_index_table_change1_minus1 */
      if(IMATRIX_FREE( tables.rho_index_table_change1_minus1 ) ) info=1;
      
      /* rho_index_table_change1_plus1 */
      if(IMATRIX_FREE( tables.rho_index_table_change1_plus1 ) ) info=1;

      /* rho_index_table */
      if(IMATRIX_FREE( tables.rho_index_table) ) info=1;

      /* dummy_rho_index2 */
      if( IVECTOR_FREE( dummy_rho_index2 ) ) info=1;

      /* dummy_rho_index1 */
      if( IVECTOR_FREE( dummy_rho_index1 ) ) info=1;

    }

    /* electron_single_index_table */
    if( IMATRIX_FREE( electron_single_index_table ) )    info=1;

    if( !flag ){

      /* electron_many_index_table */
      if( IMATRIX_FREE( electron_many_index_table ) )      info=1;
      
      /* atoms_index_table */
      if( IMATRIX_FREE( atom_index_table ) )              info=1;

      /* reduced_coordinate_index_table */
      if( IMATRIX_FREE( reduced_coordinate_index_table ) ) info=1;
      
      /* coordinate_index_table */
      if( IMATRIX_FREE( coordinate_index_table ) )         info=1;

      /* coordinate_index_table */
      if( IMATRIX_FREE( rho_index_aux_table ) ) info=1;

    }

  }


#ifdef __DEBUG__

  if( !info ){

    indexing_allocate_flag=0;

  }
  else{

    indexing_allocate_flag=1;

  }

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------
//                  COORDINATES
//------------------------------------------

/* coordinate_index_inverse */

int coordinate_index_inverse( const int index, int* i_p, int* j_p ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */


  if( MATRIX_ARG_INVERSE( index, *i_p, *j_p, N_coor_loc, N_coor_loc ) ) info=1;


  return info;

}

//------------------------------------------

/* reduced_coordinate_index_inverse */

int reduced_coordinate_index_inverse( const int index, int* i_p, int* j_p ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */


  if( MATRIX_ARG_INVERSE( index, *i_p, *j_p, N_coor_red_loc, N_coor_red_loc ) ) info=1;


  return info;

}

//------------------------------------------
//                  ATOMS
//------------------------------------------

/* atom_index_inverse */

int atom_index_inverse( const int index, int* i_p, int* j_p ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */


  if( MATRIX_ARG_INVERSE( index, *i_p, *j_p, N_atoms_loc, N_atoms_loc ) ) info=1;


  return info;

}

//------------------------------------------
//                  ELECTRONS
//------------------------------------------

/* electron_many_index_inverse */

int electron_many_index_inverse( const int index, int* i_p, int* j_p ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */


  if( MATRIX_ARG_INVERSE( index, *i_p, *j_p, N_levels_many_loc, N_levels_many_loc ) ) info=1;


  return info;

}

//------------------------------------------

/* electron_single_index_inverse */

int electron_single_index_inverse( const int index, int* i_p, int* j_p ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */


  if( MATRIX_ARG_INVERSE( index, *i_p, *j_p, N_levels_single_loc, N_levels_single_loc ) ) info=1;


  return info;

}

//------------------------------------------
//                  RHO
//------------------------------------------

/* compute_rho_dimensions */

int compute_rho_dimensions( constants_p constants_p, int CEID_order, int N_coor_red ){

  /* constants */
  int*    sqrt_max_rho_index_p;
  int*    max_rho_index_p;
  int*    rho_index_border_length_p;
  int*    rho_index_next_to_border_length_p;
  /* dummies */
  int    info=0;
  int    i;
  double dummy;


  sqrt_max_rho_index_p              = &(constants_p->sqrt_max_rho_index);
  max_rho_index_p                   = &(constants_p->max_rho_index);
  rho_index_border_length_p         = &(constants_p->rho_index_border_length);
  rho_index_next_to_border_length_p = &(constants_p->rho_index_next_to_border_length);


  // sqrt_max_rho_index

  if( CEID_order > 0 ){

    i=CEID_order;

    dummy = (double) (N_coor_red +i )/ (double) ( i );

    while( i > 1 ){

      i--;

      dummy *= (double) (N_coor_red +i )/ (double) ( i );

    }


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "sqrt_max_rho_index = %le [as a double]\n", dummy );

#endif /* __DEBUG_PLUS__ */ 


    *sqrt_max_rho_index_p = (int) (dummy+ EPS); //WARNING! it seems to be important


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "sqrt_max_rho_index = %d [as an integer]\n", *sqrt_max_rho_index_p );
    fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */ 


  }
  else if( CEID_order ==  0 ){

    *sqrt_max_rho_index_p = 1; // Ehrenfest

  }
  else{

    fprintf( stderr, "ERROR in compute_sqrt_max_rho_index: CEID_order must be non-zero\n" );
    fflush( stderr );

    *sqrt_max_rho_index_p = 0;

    info=1;

  }

  // max_rho_index
  *max_rho_index_p = ( *sqrt_max_rho_index_p ) *( *sqrt_max_rho_index_p );


  //  rho_index_border_length
  *rho_index_border_length_p = ( *sqrt_max_rho_index_p ) *N_coor_red /( CEID_order +N_coor_red );

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_index_border_length = %d\n", *rho_index_border_length_p );

#endif /* __DEBUG_PLUS__ */


  //  rho_index_next_to_border_length
  if( CEID_order +N_coor_red > 1){

    *rho_index_next_to_border_length_p = ( *rho_index_border_length_p ) *CEID_order /( CEID_order +N_coor_red -1 );
    
  }
  else{
    
    *rho_index_next_to_border_length_p = 0;
    
  }
  

#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "rho_index_next_to_border_length = %d\n", *rho_index_next_to_border_length_p );
  fprintf( stdout, "\n" );
  fflush( stdout );
  
#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------

/* compute_rho_index_table */

int compute_rho_index_table( unsigned short int flag ){

  /* dummies */
  int i, j;
  static unsigned short int flag_int  =0;
  static int super_index              =0;
  static int border_index             =0;
  static int next_to_border_index     =0;
  static int counter                  =0;
  int sum;
  int info=0;


  /* resetting */
  if( flag ){

    flag_int             =0;
    super_index          =0;
    border_index         =0;
    next_to_border_index =0;
    counter              =0;

  }


  //  fprintf( stdout, "counter = %d\n", counter );

  for( i=0; i<CEID_order_loc+2; i++ ){

    //    fprintf( stdout, "i = %d\n", i );

    /*
      fprintf( stdout, "dummy_rho_index1\n" ) ;
      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
      fprintf( stdout, "\n" );
      fflush( stdout );
    */

    dummy_rho_index1.ivector[ counter ] = i;

    //    fprintf( stdout, "dummy_rho_index1[ %d ] = %d\n", counter, dummy_rho_index1.ivector[ counter ] );

    //    fprintf( stdout, "compute sum\n" );

    sum =0;
    for( j=0; j< N_coor_red_loc; j++ ){
      
      sum += dummy_rho_index1.ivector[ j ];
      
    }
    //    fprintf( stdout, "sum = %d\n", sum );


    if( sum > CEID_order_loc ){
      
      dummy_rho_index1.ivector[ counter ] = 0;

      counter--;

      //      fprintf( stdout, "sum> CEID_order: Breaking\n" );

      break;
      
    }
    else{ 

      if( counter <  N_coor_red_loc-1 ){

	counter++;
	
	compute_rho_index_table( flag_int ); // WARNING: recurrency

      }
      else if( counter == N_coor_red_loc-1 ){

#ifdef __DEBUG_PLUS__

	fprintf( stdout, "%d \t (", super_index );

#endif /* __DEBUG_PLUS__*/

	// write table
	for( j=0; j< N_coor_red_loc; j++ ){
	  
	  tables.rho_index_table.imatrix[ super_index ][ j ] = dummy_rho_index1.ivector[ j ];

#ifdef __DEBUG_PLUS__

	  fprintf( stdout, " %d ", tables.rho_index_table.imatrix[ super_index ][ j ] );
	
#endif /* __DEBUG_PLUS__*/

	}
#ifdef __DEBUG_PLUS__

	if( sum == CEID_order_loc ){

	  fprintf( stdout, ")* \n" );

	}
	else if( sum == CEID_order_loc -1 ){

	  fprintf( stdout, ")** \n" );

	}
	else{

	  fprintf( stdout, ") \n" );

	}

#endif /* __DEBUG_PLUS__*/

	if( sum == CEID_order_loc ){

	  rho_index_border.ivector[ border_index ] = super_index;

	  border_index++;

	}
	else if( sum == CEID_order_loc -1 ){

	  rho_index_next_to_border.ivector[ next_to_border_index ] = super_index;

	  next_to_border_index++;

	}

	super_index++;

      }
      else{
	
	info=1;

      }

    } /* end if */

  } /* end for loop */


  return info;

}
//------------------------------------------

/* search_rho_index_table */

int search_rho_index_table( ivector_p index_p ){

  /* dummies */
  int i;
  int sum;
  int counter;
  int super_index=-1;
  int info=0;
  

  if( index_p->ivector_dim != N_coor_red_loc ){

    fprintf( stderr, "ERROR in search_rho_index_table: index has got the wrong dimension\n" );

    info=1;

  }

  if( !info ){

    sum =0;
    for( i=0; i< N_coor_red_loc; i++ ){
      
      if( index_p->ivector[ i ] < 0 || index_p->ivector[ i ] > CEID_order_loc ){
	
	sum = CEID_order_loc+1;
	
	break;
	
      }
      
      sum += index_p->ivector[ i ];
      
    }
    
    if( sum > CEID_order_loc ){ 
      
      super_index = sqrt_max_rho_index_loc;
      
    }
    else{
    
      super_index = 0;
      
      counter=0;
      
      while( counter < N_coor_red_loc ){
	
	if( index_p->ivector[ counter ] != tables.rho_index_table.imatrix[ super_index ][ counter ] ){
	  
	  super_index++;
	  
	}
	else{
	  
	  counter++;
	
	}
	
      } /* end while loop */
      
    } /* end if loop */
    
  } /* info loop */
  

  return super_index;
  
}

//------------------------------------------

/* compute_rho_index_other_tables */

int compute_rho_index_other_tables( void ){

  /* dummies */
  int i, j, k;
  int dummy;
  int info=0;


  /* rho_index_table_change1_plus1 */

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc; j++ ){
      
      dummy_rho_index1.ivector[ j ] = tables.rho_index_table.imatrix[ i ][ j ];
      
    }
    
    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 ) );

    for( j=0; j<N_coor_red_loc; j++ ){

      dummy_rho_index1.ivector[ j ] += 1;
      
      dummy = SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 );

      //      fprintf( stdout, "\t super_index_changed = %d\n", dummy );

      tables.rho_index_table_change1_plus1.imatrix[ i ][ j ] = dummy;

      dummy_rho_index1.ivector[ j ] -= 1;

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_index_table_change1_plus1\n" );

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc-1; j++ ){

      fprintf( stdout, "%d, ", tables.rho_index_table_change1_plus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables.rho_index_table_change1_plus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_index_table_change1_minus1 */

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc; j++ ){
      
      dummy_rho_index1.ivector[ j ] = tables.rho_index_table.imatrix[ i ][ j ];
      
    }
    
    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 ) );

    for( j=0; j<N_coor_red_loc; j++ ){

      dummy_rho_index1.ivector[ j ] -= 1;
      
      dummy = SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 );

      //      fprintf( stdout, "\t super_index_changed = %d\n", dummy );

      tables.rho_index_table_change1_minus1.imatrix[ i ][ j ] = dummy;

      dummy_rho_index1.ivector[ j ] += 1;

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_index_table_change1_minus1\n" );

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc-1; j++ ){

      fprintf( stdout, "%d, ", tables.rho_index_table_change1_minus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables.rho_index_table_change1_minus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_index_table_change1_plus2 */

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc; j++ ){
      
      dummy_rho_index1.ivector[ j ] = tables.rho_index_table.imatrix[ i ][ j ];
      
    }
    
    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 ) );

    for( j=0; j<N_coor_red_loc; j++ ){

      dummy_rho_index1.ivector[ j ] += 2;
      
      dummy = SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 );

      //      fprintf( stdout, "\t super_index_changed = %d\n", dummy );

      tables.rho_index_table_change1_plus2.imatrix[ i ][ j ] = dummy;

      dummy_rho_index1.ivector[ j ] -= 2;

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_index_table_change1_plus2\n" );

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc-1; j++ ){

      fprintf( stdout, "%d, ", tables.rho_index_table_change1_plus2.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables.rho_index_table_change1_plus2.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_index_table_change1_minus2 */

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc; j++ ){
      
      dummy_rho_index1.ivector[ j ] = tables.rho_index_table.imatrix[ i ][ j ];
      
    }
    
    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 ) );

    for( j=0; j<N_coor_red_loc; j++ ){

      dummy_rho_index1.ivector[ j ] -= 2;
      
      dummy = SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 );

      //      fprintf( stdout, "\t super_index_changed = %d\n", dummy );

      tables.rho_index_table_change1_minus2.imatrix[ i ][ j ] = dummy;

      dummy_rho_index1.ivector[ j ] += 2;

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_index_table_change1_minus2\n" );

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc-1; j++ ){

      fprintf( stdout, "%d, ", tables.rho_index_table_change1_minus2.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables.rho_index_table_change1_minus2.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_index_table_change2_plus1_plus1 */

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc; j++ ){
      
	dummy_rho_index1.ivector[ j ] = tables.rho_index_table.imatrix[ i ][ j ];
      
    }

    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 ) );

    for( j=0; j<N_coor_red_loc; j++ ){

      for( k=0; k<N_coor_red_loc; k++ ){

	dummy_rho_index1.ivector[ j ] += 1;

	dummy_rho_index1.ivector[ k ] += 1;
	
	dummy = SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 );
	
	//	fprintf( stdout, "\t super_index_changed = %d\n", dummy );
	
	tables.rho_index_table_change2_plus1_plus1.imatrix[ i ][ j +k *N_coor_red_loc ] = dummy;

	dummy_rho_index1.ivector[ k ] -= 1;

	dummy_rho_index1.ivector[ j ] -= 1;

      }

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_index_table_change2_plus1_plus1\n" );

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc*N_coor_red_loc -1; j++ ){

      fprintf( stdout, "%d, ", tables.rho_index_table_change2_plus1_plus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables.rho_index_table_change2_plus1_plus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_index_table_change2_plus1_minus1 */

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc; j++ ){
      
	dummy_rho_index1.ivector[ j ] = tables.rho_index_table.imatrix[ i ][ j ];
      
    }

    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 ) );

    for( j=0; j<N_coor_red_loc; j++ ){

      for( k=0; k<N_coor_red_loc; k++ ){

	dummy_rho_index1.ivector[ j ] += 1;

	dummy_rho_index1.ivector[ k ] -= 1;
	
	dummy = SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 );
	
	//	fprintf( stdout, "\t super_index_changed = %d\n", dummy );
	
	tables.rho_index_table_change2_plus1_minus1.imatrix[ i ][ j +k *N_coor_red_loc ] = dummy;

	dummy_rho_index1.ivector[ k ] += 1;

	dummy_rho_index1.ivector[ j ] -= 1;

      }

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_index_table_change2_plus1_minus1\n" );

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc*N_coor_red_loc -1; j++ ){

      fprintf( stdout, "%d, ", tables.rho_index_table_change2_plus1_minus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables.rho_index_table_change2_plus1_minus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_index_table_change2_minus1_minus1 */

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc; j++ ){
      
	dummy_rho_index1.ivector[ j ] = tables.rho_index_table.imatrix[ i ][ j ];
      
    }

    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 ) );

    for( j=0; j<N_coor_red_loc; j++ ){

      for( k=0; k<N_coor_red_loc; k++ ){

	dummy_rho_index1.ivector[ j ] -= 1;

	dummy_rho_index1.ivector[ k ] -= 1;
	
	dummy = SEARCH_RHO_INDEX_TABLE( dummy_rho_index1 );
	
	//	fprintf( stdout, "\t super_index_changed = %d\n", dummy );
	
	tables.rho_index_table_change2_minus1_minus1.imatrix[ i ][ j +k *N_coor_red_loc ] = dummy;

	dummy_rho_index1.ivector[ k ] += 1;

	dummy_rho_index1.ivector[ j ] += 1;

      }

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_index_table_change2_minus1_minus1\n" );

  for( i=0; i<sqrt_max_rho_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc*N_coor_red_loc -1; j++ ){

      fprintf( stdout, "%d, ", tables.rho_index_table_change2_minus1_minus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables.rho_index_table_change2_minus1_minus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------

/* rho_index_inverse */

int rho_index_inverse( const int index, ivector_p index1_p, ivector_p index2_p ){

  /* dummies */
  int dummy_superindex1;
  int dummy_superindex2;
  int info=0;


#ifdef __DEBUG__

  if( !info ){

    if( IVECTOR_CONSINSTENCY( *index1_p, *index2_p ) ){

      fprintf( stderr, "ERROR: indexing inconsistency.\n");
      fflush( stderr );

      info=1;

    }

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_index_loc, sqrt_max_rho_index_loc ) ) info=1;

    if( !info ){

      if( RHO_INDEX_INVERSE_AUX( dummy_superindex1, *index1_p ) ) info=1;  

      if( RHO_INDEX_INVERSE_AUX( dummy_superindex2, *index2_p ) ) info=1;  

    }

  }


  return info;

}
//------------------------------------------

/* rho_index_inverse_aux */

int rho_index_inverse_aux( const int index, ivector_p index_p ){

  /* dummies */
  int i;
  int info=0;


  if( !info ){
    
    for( i=0; i<N_coor_red_loc; i++ ){
      
      index_p->ivector[ i ] = tables.rho_index_table.imatrix[ index ][ i ];
      
    }
    
  }


  return info;

}

//------------------------------------------

/* check indexing */

int check_indexing( void ){

  /* dummies */
  ivector index_1, index_2;
  ivector index_1_tmp, index_2_tmp;
  int index, index_aux;
  int i, j;
  int sum;
  int info=0;


#ifdef __DEBUG__

  if( !indexing_allocate_flag ){
    
    info=1;
    
    fprintf( stderr, "ERROR: indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

//------------------------------------------
//                  COORDINATES
//------------------------------------------

#ifdef __DEBUG__

  fprintf( stdout, "CHECKING COORDINATE_INDEX\n");
  fflush( stdout );
  
#endif /* __DEBUG__ */

  for( index=0; index<(N_coor_loc*N_coor_loc); index++ ){

    /*
      fprintf( stdout, "->index:     %d\n", index );
      fflush( stdout );
    */

    if( COORDINATE_INDEX_INVERSE( index, i, j ) ) info=1;

    if( !info ){

      /*
	fprintf( stdout, "i = %d \t j = %d\n", i, j );
	fflush( stdout );
      */

      index_aux =  COORDINATE_INDEX( i, j );

      /*      
	      fprintf( stdout, "->index_aux:     %d\n", index_aux );
	      fflush( stdout );
      */

      if( index_aux != index ){
	
	fprintf( stderr, "CONSISTENCY ERROR for coordinate_index %d\n", index );
	fflush( stderr );

	fprintf( stdout, "CONSISTENCY ERROR for coordinate index %d\n", index );
	fflush( stdout );

	info=1;

	break;
	
      }

    }    
    
  }

  //------------------------------------------
  //                  REDUCED_COORDINATES
  //------------------------------------------
  
#ifdef __DEBUG__

  fprintf( stdout, "CHECKING REDUCED_COORDINATE_INDEX\n");
  fflush( stdout );
  
#endif /* __DEBUG__ */

  for( index=0; index<(N_coor_red_loc*N_coor_red_loc); index++ ){

    /*
      fprintf( stdout, "->index:     %d\n", index );
      fflush( stdout );
    */

    if( REDUCED_COORDINATE_INDEX_INVERSE( index, i, j ) ) info=1;

    if( !info ){

      /*
	fprintf( stdout, "i = %d \t j = %d\n", i, j );
	fflush( stdout );
      */

      index_aux =  REDUCED_COORDINATE_INDEX( i, j );

      /*      
	      fprintf( stdout, "->index_aux:     %d\n", index_aux );
	      fflush( stdout );
      */

      if( index_aux != index ){
	
	fprintf( stderr, "CONSISTENCY ERROR for reduced_coordinate_index %d\n", index );
	fflush( stderr );

	fprintf( stdout, "CONSISTENCY ERROR for reduced_coordinate index %d\n", index );
	fflush( stdout );

	info=1;

	break;
	
      }

    }    
    
  }

//------------------------------------------
//                  RHO
//------------------------------------------

#ifdef __DEBUG__

  fprintf( stdout, "CHECKING RHO_INDEX\n");
  fflush( stdout );
  
#endif /* __DEBUG__ */

  
  fprintf( stdout, "\n" );
  fprintf( stdout, "----------\n" );
  fprintf( stdout, "RHO_INDEX_TABLE\n" );
  
  for( index=0; index<sqrt_max_rho_index_loc; index++ ){
      
    fprintf( stdout, "%d \t( ", index );
    
    sum=0;
    for( j=0; j<N_coor_red_loc; j++ ){
      
      fprintf( stdout, "%d ", tables.rho_index_table.imatrix[ index ][ j ] );

      sum += tables.rho_index_table.imatrix[ index ][ j ];

    }


    if( sum == CEID_order_loc ){
      
      fprintf( stdout, ")* \n" );
      
    }
    else if( sum == CEID_order_loc -1 ){
      
      fprintf( stdout, ")** \n" );
      
	}
    else{
      
      fprintf( stdout, ") \n" );

    }

  }

  fprintf( stdout, "\n" );

  fprintf( stdout, "rho_index_border\n" );
  if( IVECTOR_PRINT_PLUS( stdout, rho_index_border ) ) info=1;

  if( rho_index_next_to_border_length_loc > 0 ){

    fprintf( stdout, "rho_index_next_to_border\n" );
    if( IVECTOR_PRINT_PLUS( stdout, rho_index_next_to_border ) ) info=1;

  }

  if( !info ){

    if( IVECTOR_ALLOCATE( N_coor_red_loc, index_1 ) ) info=1;
    if( IVECTOR_ALLOCATE( N_coor_red_loc, index_2 ) ) info=1;

  }

  fprintf( stdout, "\n" );
  fprintf( stdout, "----------\n" );
  fprintf( stdout, "RHO_INDICES\n" );
  
  for( index=0; index<max_rho_index_loc; index++ ){
    
    fprintf( stdout, "%d \t( ", index );

    if( RHO_INDEX_INVERSE( index, index_1, index_2) ) info=1;

    for( j=0; j<N_coor_red_loc; j++ ){

      fprintf( stdout, "%d ", index_1.ivector[ j ] );

    }

    fprintf( stdout, ")\t( " );

    for( j=0; j<N_coor_red_loc; j++ ){

      fprintf( stdout, "%d ", index_2.ivector[ j ] );

    }

    fprintf( stdout, ")\t\t" );

    index_aux = RHO_INDEX( index_2, index_1 ); 

    if( index_aux != rho_index_conjugate.ivector[ index ] ){

      fprintf( stderr, "ERROR: the conjugate index of %d is not %d [%d]\n", index, index_aux, rho_index_conjugate.ivector[ index ] );

      break;

      info=1;

    }


    fprintf( stdout, "%d\n", index_aux );

  }

  fprintf( stdout, "----------\n");
  fprintf( stdout, "\n" );


  //  -------------------------------------------------
  //  -------------------------------------------------


  if( !info ){

    if( IVECTOR_ALLOCATE( N_coor_red_loc, index_1_tmp ) ) info=1;
    if( IVECTOR_ALLOCATE( N_coor_red_loc, index_2_tmp ) ) info=1;

  }

  if( !info ){

    if( IVECTOR_CONSINSTENCY( index_1, index_2 ) ) info=1;

  }

  for( index=0; index<max_rho_index_loc; index++ ){

    /*    
	  fprintf( stdout, "->index:     %d\n", index );
	  fflush( stdout );
    */

    if( RHO_INDEX_INVERSE( index, index_1, index_2) ) info=1;

    /*
      fprintf( stdout, "print index_1\n" ) ;
      fflush( stdout );
      if( IVECTOR_PRINT_PLUS( stdout, index_1) ) info=1;
      
      fprintf( stdout, "print index_2\n" ) ;
      fflush( stdout );
      if( IVECTOR_PRINT_PLUS( stdout, index_2) ) info=1;
      
      fprintf( stdout, "-------------------------------------\n");
      fflush( stdout );
    */

    /* ---------------------------------- */
    // FIRST
    /* ---------------------------------- */

    // check: RHO_INDEX_CHANGE1_PLUS1_FIRST
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE1_PLUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*    
	    fprintf( stdout, "-------------------------------------\n");
	    fprintf( stdout, "---> i=%d\n\n", i );
	    fflush( stdout );
      */

	index_aux = RHO_INDEX_CHANGE1_PLUS1_FIRST( index, i );

	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/

      if( index_aux < max_rho_index_loc ){
	
	if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1_tmp\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2_tmp\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_INDEX_CHANGE1_MINUS1_FIRST( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stdout );
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE1_MINUS1_FIRST
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE1_MINUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_INDEX_CHANGE1_MINUS1_FIRST( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_index_loc ){
	
	if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_INDEX_CHANGE1_PLUS1_FIRST( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stdout );
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE1_PLUS2_FIRST
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE1_PLUS2_FIRST\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_INDEX_CHANGE1_PLUS2_FIRST( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_index_loc ){
	
	if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_INDEX_CHANGE1_MINUS2_FIRST( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stdout );
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE1_MINUS2_FIRST
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE1_MINUS2_FIRST\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_INDEX_CHANGE1_MINUS2_FIRST( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_index_loc ){
	
	if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_INDEX_CHANGE1_PLUS2_FIRST( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stdout );
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE2_PLUS1_PLUS1_FIRST
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE2_PLUS1_PLUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	
	index_aux = RHO_INDEX_CHANGE2_PLUS1_PLUS1_FIRST( index, i, j );
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
       
	if( i != j && index_aux < max_rho_index_loc ){
	
	  if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*	
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index_aux, i, j );
	  
	  /*
	    fprintf( stdout, "index_aux [check]: %d\n\n", index_aux );
	    fflush( stdout );
	  */

	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stdout );
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	

	index_aux = RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( index, i, j );
	
	
	/*
	  fprintf( stdout, "index_aux: %d [before]\n\n", index_aux );
	  fflush( stdout );
	*/
	
	if( i != j && index_aux < max_rho_index_loc ){
	
	  if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	    fflush( stdout );
	  */	  

	  index_aux = RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( index_aux, j, i ); //WARNING: note the index inversion
	  
	  /*
	    fprintf( stdout, "index_aux: %d [after]\n\n", index_aux );
	    fflush( stdout );
	  */

	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stdout );
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	

	index_aux = RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index, i, j );
	
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
	
	if( i != j && index_aux < max_rho_index_loc ){
	
	  if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_INDEX_CHANGE2_PLUS1_PLUS1_FIRST( index_aux, i, j );
	  
	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stdout );
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    /* ---------------------------------- */

    /* ---------------------------------- */
    // SECOND
    /* ---------------------------------- */

    // check: RHO_INDEX_CHANGE1_PLUS1_SECOND
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE1_PLUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*    
	    fprintf( stdout, "-------------------------------------\n");
	    fprintf( stdout, "---> i=%d\n\n", i );
	    fflush( stdout );
      */

	index_aux = RHO_INDEX_CHANGE1_PLUS1_SECOND( index, i );

	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/

      if( index_aux < max_rho_index_loc ){
	
	if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1_tmp\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2_tmp\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_INDEX_CHANGE1_MINUS1_SECOND( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stdout );
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE1_MINUS1_SECOND
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE1_MINUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_INDEX_CHANGE1_MINUS1_SECOND( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_index_loc ){
	
	if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_INDEX_CHANGE1_PLUS1_SECOND( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stdout );
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE1_PLUS2_SECOND
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE1_PLUS2_SECOND\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_INDEX_CHANGE1_PLUS2_SECOND( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_index_loc ){
	
	if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_INDEX_CHANGE1_MINUS2_SECOND( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stdout );
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE1_MINUS2_SECOND
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE1_MINUS2_SECOND\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_INDEX_CHANGE1_MINUS2_SECOND( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_index_loc ){
	
	if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_INDEX_CHANGE1_PLUS2_SECOND( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stdout );
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE2_PLUS1_PLUS1_SECOND
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE2_PLUS1_PLUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	
	index_aux = RHO_INDEX_CHANGE2_PLUS1_PLUS1_SECOND( index, i, j );
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
       
	if( i != j && index_aux < max_rho_index_loc ){
	
	  if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*	
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_INDEX_CHANGE2_MINUS1_MINUS1_SECOND( index_aux, i, j );
	  
	  /*
	    fprintf( stdout, "index_aux [check]: %d\n\n", index_aux );
	    fflush( stdout );
	  */

	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stdout );
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE2_PLUS1_MINUS1_SECOND
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE2_PLUS1_MINUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	

	index_aux = RHO_INDEX_CHANGE2_PLUS1_MINUS1_SECOND( index, i, j );
	
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
	
	if( i != j && index_aux < max_rho_index_loc ){
	
	  if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_INDEX_CHANGE2_PLUS1_MINUS1_SECOND( index_aux, j, i ); //WARNING: note the index inversion
	  
	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stdout );
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_INDEX_CHANGE2_MINUS1_MINUS1_SECOND
    //    fprintf( stdout, "----> RHO_INDEX_CHANGE2_MINUS1_MINUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	

	index_aux = RHO_INDEX_CHANGE2_MINUS1_MINUS1_SECOND( index, i, j );
	
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
	
	if( i != j && index_aux < max_rho_index_loc ){
	
	  if( RHO_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_INDEX_CHANGE2_PLUS1_PLUS1_SECOND( index_aux, i, j );
	  
	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stdout );
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */

  } /* index loop */

  /* ---------------------------------- */
  /* ---------------------------------- */

#ifdef __DEBUG__

  fprintf( stdout, "CHECKING DIAGONAL_RHO_INDICES\n" );
  fflush( stdout );
    
#endif /* __DEBUG__ */

  for( index=0; index<diagonal_rho_indices.ivector_dim; index++ ){

    /*
      fprintf( stdout, "->index:     %d\n", index );
      fflush( stdout );
    */

    index_aux = diagonal_rho_indices.ivector[ index ];

    /*
      fprintf( stdout, "-->index_aux:     %d\n\n", index_aux );
      fflush( stdout );
    */

    if( RHO_INDEX_INVERSE( index_aux, index_1, index_2 ) ) info=1;

    /*
      fprintf( stdout, "print index_1\n" ) ;
      fflush( stdout );
      if( IVECTOR_PRINT_PLUS( stdout, index_1) ) info=1;

      fprintf( stdout, "print index_2\n" ) ;
      fflush( stdout );
      if( IVECTOR_PRINT_PLUS( stdout, index_2) ) info=1;
    */

    if( IVECTOR_ARE_EQUAL( index_1, index_2 ) ){

      fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
      fflush( stderr );
      
      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
      fflush( stdout );
		 
      info=1;

      break;

    }

  }

  if( !info ){

    if( IVECTOR_FREE( index_1 ) ) info=1;
    if( IVECTOR_FREE( index_2 ) ) info=1;

    if( IVECTOR_FREE( index_1_tmp ) ) info=1;
    if( IVECTOR_FREE( index_2_tmp ) ) info=1;

  }

  return info;

}

//------------------------------------------

