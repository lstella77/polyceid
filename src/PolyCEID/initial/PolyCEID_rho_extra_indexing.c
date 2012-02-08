
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
#include "PolyCEID_rho_extra_indexing.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* private variables */

#ifdef __DEBUG__

unsigned short int rho_extra_indexing_allocate_flag=0;

#endif /* __DEBUG__ */

int  N_coor_loc_extra;
int  N_coor_red_loc_extra;
int  max_rho_extra_index_loc;
int  sqrt_max_rho_extra_index_loc;
int  CEID_order_loc_extra;


typedef struct{

  imatrix rho_extra_index_table;                       // [sqrt_max_rho_extra_index, N_coor]

  imatrix rho_extra_index_table_change1_plus1;         // [sqrt_max_rho_extra_index, N_coor]

  imatrix rho_extra_index_table_change1_minus1;        // [sqrt_max_rho_extra_index, N_coor]

  imatrix rho_extra_index_table_change1_plus2;         // [sqrt_max_rho_extra_index, N_coor]

  imatrix rho_extra_index_table_change1_minus2;        // [sqrt_max_rho_extra_index, N_coor]

  imatrix rho_extra_index_table_change2_plus1_plus1;   // [sqrt_max_rho_extra_index, N_coor *N_coor]

  imatrix rho_extra_index_table_change2_plus1_minus1;  // [sqrt_max_rho_extra_index, N_coor *N_coor]

  imatrix rho_extra_index_table_change2_minus1_minus1; // [sqrt_max_rho_extra_index, N_coor *N_coor]


} rho_extra_index_tables;

rho_extra_index_tables tables_extra;


//------------------------------------------
//                  RHO_EXTRA
//------------------------------------------

/* rho_extra_index */

/*

  2st step: from two ivectors of dimension N_coor_red_extra
  to a single scalar index

*/

int rho_extra_index( ivector index1, ivector index2 ){

  /* dummies */
  int dummy_superindex1;
  int dummy_superindex2;
  int index=-1;
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( !info ){

    if( IVECTOR_CONSINSTENCY( index1, index2 ) ){

      fprintf( stderr, "ERROR: rho_extra_indexing inconsistency.\n");
      fflush( stderr );

      info=1;
      
    }

  }

#endif /* __DEBUG__ */


  if( !info ){

    dummy_superindex1 = SEARCH_RHO_EXTRA_INDEX_TABLE( index1 ); // WARNING: time-consuming operation

    dummy_superindex2 = SEARCH_RHO_EXTRA_INDEX_TABLE( index2 ); // WARNING: time-consuming operation

    if( dummy_superindex1 < 0 || dummy_superindex1 > ( sqrt_max_rho_extra_index_loc -1 ) || 
	dummy_superindex2 < 0 || dummy_superindex2 > ( sqrt_max_rho_extra_index_loc -1 ) ){ // WARNING: it's necessary

      index = max_rho_extra_index_loc;

    }
    else{

      index = RHO_EXTRA_INDEX_AUX( dummy_superindex1, dummy_superindex2 );

    }

  } /* end info conditional */


  return index;


}

//------------------------------------------

/* rho_extra_index_aux */

int rho_extra_index_aux( int index1, int index2 ){

  /* dummies */
  int index=-1;
  int info=0;


  if( !info ){

    //index = MATRIX_ARG( index1, index2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
    index = rho_extra_index_aux_table.imatrix[index1][index2];

  }


  return index;


}

//------------------------------------------
//------------------------------------------

/* rho_extra_index_change1_plus1_aux */

int rho_extra_index_change1_plus1_aux( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_extra_index_loc-1 ){

      index_changed = sqrt_max_rho_extra_index_loc;

    }
    else{

      index_changed = tables_extra.rho_extra_index_table_change1_plus1.imatrix[ index][ i ];

    }

  }

  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change1_minus1_aux */

int rho_extra_index_change1_minus1_aux( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_extra_index_loc-1 ){

      index_changed = sqrt_max_rho_extra_index_loc;

    }
    else{

      index_changed = tables_extra.rho_extra_index_table_change1_minus1.imatrix[ index][ i ];
    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change1_plus2_aux */

int rho_extra_index_change1_plus2_aux( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_extra_index_loc-1 ){

      index_changed = sqrt_max_rho_extra_index_loc;

    }
    else{
      
      index_changed = tables_extra.rho_extra_index_table_change1_plus2.imatrix[ index][ i ];
    
    }

  }


  return index_changed;

}

//------------------------------------------

/*  rho_extra_index_change1_minus2_aux */

int rho_extra_index_change1_minus2_aux( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_extra_index_loc-1 ){
      
      index_changed = sqrt_max_rho_extra_index_loc;
      
    }
    else{

      index_changed = tables_extra.rho_extra_index_table_change1_minus2.imatrix[ index][ i ];

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_plus1_plus1_aux */

int rho_extra_index_change2_plus1_plus1_aux( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

  if( j < 0 || j > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: j value is not legal [%d]\n", j );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_extra_index_loc-1 ){
      
      index_changed = sqrt_max_rho_extra_index_loc;
      
    }
    else{

      index_changed = tables_extra.rho_extra_index_table_change2_plus1_plus1.imatrix[ index ][ REDUCED_COORDINATE_INDEX( i, j) ];

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_plus1_minus1_aux */

int rho_extra_index_change2_plus1_minus1_aux( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault 
  int info=0;

#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

  if( j < 0 || j > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: j value is not legal [%d]\n", j );

    info=1;

  }


#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_extra_index_loc-1 ){

      index_changed = sqrt_max_rho_extra_index_loc;

    }
    else{

      index_changed = tables_extra.rho_extra_index_table_change2_plus1_minus1.imatrix[ index ][ REDUCED_COORDINATE_INDEX( i, j) ];

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_minus1_minus1_aux */

int rho_extra_index_change2_minus1_minus1_aux( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

  if( j < 0 || j > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: j value is not legal [%d]\n", j );

    info=1;

  }


#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > sqrt_max_rho_extra_index_loc-1 ){

      index_changed = sqrt_max_rho_extra_index_loc;

    }
    else{

      index_changed = tables_extra.rho_extra_index_table_change2_minus1_minus1.imatrix[ index ][ REDUCED_COORDINATE_INDEX( i, j) ];

    }

  }


  return index_changed;

}

//------------------------------------------
//------------------------------------------

/* rho_extra_index_change1_plus1_first */

int rho_extra_index_change1_plus1_first( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( dummy_superindex1, i );      
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change1_minus1_first */

int rho_extra_index_change1_minus1_first( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_EXTRA_INDEX_CHANGE1_MINUS1_AUX( dummy_superindex1, i );

	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_extra_index_loc-1 ){
	  
	  index_changed = max_rho_extra_index_loc;
	  
	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change1_plus2_first */

int rho_extra_index_change1_plus2_first( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( dummy_superindex1, i );
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;


}

//------------------------------------------

/*  rho_extra_index_change1_minus2_first */

int rho_extra_index_change1_minus2_first( int index, int i ){


  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;
      
    }
    else{


      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_EXTRA_INDEX_CHANGE1_MINUS2_AUX( dummy_superindex1, i );
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{
	  
	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_plus1_plus1_first */

int rho_extra_index_change2_plus1_plus1_first( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( dummy_superindex1, i, j );
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_plus1_minus1_first */

int rho_extra_index_change2_plus1_minus1_first( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_AUX( dummy_superindex1, i, j );
	
	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_minus1_minus1_first */

int rho_extra_index_change2_minus1_minus1_first( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex1_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex1_changed = RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_AUX( dummy_superindex1, i, j );

	if( dummy_superindex1_changed < 0 || dummy_superindex1_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{
	
	  //index_changed = MATRIX_ARG( dummy_superindex1_changed, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1_changed][dummy_superindex2];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------
//------------------------------------------

/* rho_extra_index_change1_plus1_second */

int rho_extra_index_change1_plus1_second( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS1_AUX( dummy_superindex2, i );      
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change1_minus1_second */

int rho_extra_index_change1_minus1_second( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_EXTRA_INDEX_CHANGE1_MINUS1_AUX( dummy_superindex2, i );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change1_plus2_second */

int rho_extra_index_change1_plus2_second( int index, int i ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;

    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_EXTRA_INDEX_CHANGE1_PLUS2_AUX( dummy_superindex2, i );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{
	  
	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;


}

//------------------------------------------

/*  rho_extra_index_change1_minus2_second */

int rho_extra_index_change1_minus2_second( int index, int i ){


  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;
      
    }
    else{


      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_EXTRA_INDEX_CHANGE1_MINUS2_AUX( dummy_superindex2, i );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_plus1_plus1_second */

int rho_extra_index_change2_plus1_plus1_second( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_AUX( dummy_superindex2, i, j );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_plus1_minus1_second */

int rho_extra_index_change2_plus1_minus1_second( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_AUX( dummy_superindex2, i, j );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

	}

      }

    }

  }


  return index_changed;

}

//------------------------------------------

/* rho_extra_index_change2_minus1_minus1_second */

int rho_extra_index_change2_minus1_minus1_second( int index, int i, int j ){

  /* dummies */
  int index_changed=-1; //WARNING: default value would give a segmentation fault
  int dummy_superindex1, dummy_superindex2;
  int dummy_superindex2_changed;
  int info=0;


#ifdef __DEBUG__

  if( i < 0 || i > N_coor_loc_extra-1 ){

    fprintf( stderr, "ERROR: i value is not legal [%d]\n", i );

    info=1;

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( index < 0 || index > max_rho_extra_index_loc-1 ){

      index_changed = max_rho_extra_index_loc;
      
    }
    else{

      if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;
      
      if( !info ){
	
	dummy_superindex2_changed = RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_AUX( dummy_superindex2, i, j );
	
	if( dummy_superindex2_changed < 0 || dummy_superindex2_changed > sqrt_max_rho_extra_index_loc-1 ){

	  index_changed = max_rho_extra_index_loc;

	}
	else{

	  //index_changed = MATRIX_ARG( dummy_superindex1, dummy_superindex2_changed, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );
	  index_changed = rho_extra_index_aux_table.imatrix[dummy_superindex1][dummy_superindex2_changed];

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

/* initialise rho_extra_indexing  */

int rho_extra_indexing_allocate( constants constants ){

  /* dummies */
  int i, j;
  double dummy;
  int info=0;
  

#ifdef __DEBUG__

  if( rho_extra_indexing_allocate_flag ){
    
    info=1;

    fprintf( stderr, "ERROR: rho_extra_indexing allocation already done.\n");
    fflush( stderr );
    
  }

#endif /* __DEBUG__ */

  if( !info ){

    N_coor_loc_extra                    = constants.N_coor;
    N_coor_red_loc_extra                = constants.N_coor_red;
    CEID_order_loc_extra                = constants.CEID_order+2; //WARNING: notice the +2


    // sqrt_max_rho_extra_index_loc
    
    i=CEID_order_loc_extra;

    dummy = (double) (N_coor_red_loc_extra +i )/ (double) ( i );

    while( i > 1 ){
      
      i--;
      
      dummy *= (double) (N_coor_red_loc_extra +i )/ (double) ( i );

    }
    

#ifdef __DEBUG_PLUS__

    fprintf( stdout, "sqrt_max_rho_extra_index_loc = %le [as a double]\n", dummy );

#endif /* __DEBUG_PLUS__ */ 


    sqrt_max_rho_extra_index_loc = (int) (dummy+ EPS); //WARNING! it seems to be important


#ifdef __DEBUG_PLUS__
    
    fprintf( stdout, "sqrt_max_rho_extra_index_loc = %d [as an integer]\n", sqrt_max_rho_extra_index_loc );
    fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */ 

    
    // max_rho_extra_index_loc
    max_rho_extra_index_loc = sqrt_max_rho_extra_index_loc *sqrt_max_rho_extra_index_loc;


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "max_rho_extra_index_loc = %d [as an integer]\n", max_rho_extra_index_loc );
    fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */ 


    /* rho_extra_index_aux_table */
    if( IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc, rho_extra_index_aux_table ) ) info=1;

    for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){
      for( j=0; j<sqrt_max_rho_extra_index_loc; j++ ){

        rho_extra_index_aux_table.imatrix[i][j] = MATRIX_ARG( i, j, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc );

      }
    }

    /*   
      fprintf( stdout, "rho_extra_index_aux_table\n" );
      if(IMATRIX_PRINT_PLUS( stdout, rho_extra_index_aux_table ) ) info=1;
      fprintf( stdout, "\n" );
    */

    /* dummy_rho_extra_index1 */
    if( IVECTOR_ALLOCATE( N_coor_red_loc_extra, dummy_rho_extra_index1 ) ) info=1;

    /* dummy_rho_extra_index2 */
    if( IVECTOR_ALLOCATE( N_coor_red_loc_extra, dummy_rho_extra_index2 ) ) info=1;

    /* rho_index_aux_translate */
    if( IVECTOR_ALLOCATE( constants.sqrt_max_rho_index, rho_index_aux_translate ) ) info=1;

    /* rho_extra_index_table */
    if(IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, N_coor_red_loc_extra, tables_extra.rho_extra_index_table ) ) info=1;

    /* rho_extra_index_table_change1_plus1 */
    if(IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, N_coor_red_loc_extra, tables_extra.rho_extra_index_table_change1_plus1 ) ) info=1;

    /* rho_extra_index_table_change1_minus1 */
    if(IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, N_coor_red_loc_extra, tables_extra.rho_extra_index_table_change1_minus1 ) ) info=1;

    /* rho_extra_index_table_change1_plus2 */
    if(IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, N_coor_red_loc_extra, tables_extra.rho_extra_index_table_change1_plus2 ) ) info=1;

    /* rho_extra_index_table_change1_minus2 */
    if(IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, N_coor_red_loc_extra, tables_extra.rho_extra_index_table_change1_minus2 ) ) info=1;

    /* rho_extra_index_table_change2_plus1_plus1 */
    if(IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, N_coor_red_loc_extra *N_coor_red_loc_extra, tables_extra.rho_extra_index_table_change2_plus1_plus1 ) ) info=1;

    /* rho_extra_index_table_change2_plus1_minus1 */
    if(IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, N_coor_red_loc_extra *N_coor_red_loc_extra, tables_extra.rho_extra_index_table_change2_plus1_minus1 ) ) info=1;

    /* rho_extra_index_table_change2_minus1_minus1 */
    if(IMATRIX_ALLOCATE( sqrt_max_rho_extra_index_loc, N_coor_red_loc_extra *N_coor_red_loc_extra, tables_extra.rho_extra_index_table_change2_minus1_minus1 ) ) info=1;

  }


#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "DOING: compute_rho_extra_index_table\n" );

#endif /* __DEBUG_PLUS__ */

  if( !info ){

    if( COMPUTE_RHO_EXTRA_INDEX_TABLE( 1 ) ) info=1; //WARNING: it restarts from scratch

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "\n" );
  fprintf( stdout, "DOING: compute_rho_extra_index_other_tables\n" );

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    if( COMPUTE_RHO_EXTRA_INDEX_OTHER_TABLES() ) info=1;

  }


#ifdef __DEBUG__

  if( !info ){

    rho_extra_indexing_allocate_flag=1;

  }
  else{

    rho_extra_indexing_allocate_flag=0;

  }

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------

/* free rho_extra_indexing  */

int rho_extra_indexing_free( void ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;

    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */


  if( !info ){

    /* rho_extra_index_table_change2_minus1_minus1 */
    if(IMATRIX_FREE( tables_extra.rho_extra_index_table_change2_minus1_minus1 ) ) info=1;
    
    /* rho_extra_index_table_change2_plus1_minus1 */
    if(IMATRIX_FREE( tables_extra.rho_extra_index_table_change2_plus1_minus1 ) ) info=1;
    
    /* rho_extra_index_table_change2_plus1_plus1 */
    if(IMATRIX_FREE( tables_extra.rho_extra_index_table_change2_plus1_plus1 ) ) info=1;
    
    /* rho_extra_index_table_change1_minus2 */
    if(IMATRIX_FREE( tables_extra.rho_extra_index_table_change1_minus2 ) ) info=1;
    
    /* rho_extra_index_table_change1_plus2 */
    if(IMATRIX_FREE( tables_extra.rho_extra_index_table_change1_plus2 ) ) info=1;
    
    /* rho_extra_index_table_change1_minus1 */
    if(IMATRIX_FREE( tables_extra.rho_extra_index_table_change1_minus1 ) ) info=1;
    
    /* rho_extra_index_table_change1_plus1 */
    if(IMATRIX_FREE( tables_extra.rho_extra_index_table_change1_plus1 ) ) info=1;

    /* rho_extra_index_table */
    if(IMATRIX_FREE( tables_extra.rho_extra_index_table) ) info=1;

    /* rho_index_aux_translate */
    if( IVECTOR_FREE( rho_index_aux_translate ) ) info=1;

    /* dummy_rho_extra_index2 */
    if( IVECTOR_FREE( dummy_rho_extra_index2 ) ) info=1;

    /* dummy_rho_extra_index1 */
    if( IVECTOR_FREE( dummy_rho_extra_index1 ) ) info=1;

    /* rho_extra_index_aux_table */
    if( IMATRIX_FREE( rho_extra_index_aux_table ) ) info=1;

  }


#ifdef __DEBUG__

  if( !info ){

    rho_extra_indexing_allocate_flag=0;

  }
  else{

    rho_extra_indexing_allocate_flag=1;

  }

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------
//                  RHO_EXTRA
//------------------------------------------

/* compute_rho_extra_index_table */

int compute_rho_extra_index_table( unsigned short int flag ){

  /* dummies */
  int i, j;
  static unsigned short int flag_int  =0;
  static int super_index              =0;
  static int counter                  =0;
  int sum;
  int info=0;


  /* resetting */
  if( flag ){

    flag_int             =0;
    super_index          =0;
    counter              =0;

  }


  //  fprintf( stdout, "counter = %d\n", counter );

  for( i=0; i<CEID_order_loc_extra+2; i++ ){

    //    fprintf( stdout, "i = %d\n", i );

    dummy_rho_extra_index1.ivector[ counter ] = i;

    //    fprintf( stdout, "dummy_rho_extra_index1[ %d ] = %d\n", counter, dummy_rho_extra_index1.ivector[ counter ] );

    //    fprintf( stdout, "compute sum\n" );

    sum =0;
    for( j=0; j< N_coor_red_loc_extra; j++ ){
      
      sum += dummy_rho_extra_index1.ivector[ j ];
      
    }
    //    fprintf( stdout, "sum = %d\n", sum );


    if( sum > CEID_order_loc_extra ){
      
      dummy_rho_extra_index1.ivector[ counter ] = 0;

      counter--;

      //      fprintf( stdout, "sum> CEID_order_loc_extra: Breaking\n" );

      break;
      
    }
    else{ 

      if( counter <  N_coor_red_loc_extra-1 ){

	counter++;
	
	compute_rho_extra_index_table( flag_int ); // WARNING: recurrency

      }
      else if( counter == N_coor_red_loc_extra-1 ){

#ifdef __DEBUG_PLUS__

	fprintf( stdout, "%d \t (", super_index );

#endif /* __DEBUG_PLUS__*/

	// write table
	for( j=0; j< N_coor_red_loc_extra; j++ ){
	  
	  tables_extra.rho_extra_index_table.imatrix[ super_index ][ j ] = dummy_rho_extra_index1.ivector[ j ];

#ifdef __DEBUG_PLUS__

	  fprintf( stdout, " %d ", tables_extra.rho_extra_index_table.imatrix[ super_index ][ j ] );
	
#endif /* __DEBUG_PLUS__*/

	}
#ifdef __DEBUG_PLUS__

	if( sum == CEID_order_loc_extra -2 ){

	  fprintf( stdout, ")* \n" );

	}
	else if( sum == CEID_order_loc_extra -3 ){

	  fprintf( stdout, ")** \n" );

	}
	else{

	  fprintf( stdout, ") \n" );

	}

#endif /* __DEBUG_PLUS__*/


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

/* search_rho_extra_index_table */

int search_rho_extra_index_table( ivector_p index_p ){

  /* dummies */
  int i;
  int sum;
  int counter;
  int super_index=-1;
  int info=0;

  
  if( index_p->ivector_dim != N_coor_red_loc_extra ){

    fprintf( stderr, "ERROR in search_rho_extra_index_table: index has got the wrong dimension\n" );

    info=1;

  }

  if( !info ){

    sum =0;
    for( i=0; i< N_coor_red_loc_extra; i++ ){
      
      if( index_p->ivector[ i ] < 0 || index_p->ivector[ i ] > CEID_order_loc_extra ){
	
	sum = CEID_order_loc_extra+1;
	
	break;
	
      }
      
      sum += index_p->ivector[ i ];
      
    }
    
    if( sum > CEID_order_loc_extra ){ 
      
      super_index = sqrt_max_rho_extra_index_loc;
      
    }
    else{
    
      super_index = 0;
      
      counter=0;
      
      while( counter < N_coor_red_loc_extra ){
	
	if( index_p->ivector[ counter ] != tables_extra.rho_extra_index_table.imatrix[ super_index ][ counter ] ){
	  
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

/* compute_rho_extra_index_other_tables */

int compute_rho_extra_index_other_tables ( void ){

  /* dummies */
  int i, j, k;
  int dummy;
  int info=0;


  /* rho_extra_index_table_change1_plus1 */

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc_extra; j++ ){
      
      dummy_rho_extra_index1.ivector[ j ] = tables_extra.rho_extra_index_table.imatrix[ i ][ j ];
      
    }
    
    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 ) );

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      dummy_rho_extra_index1.ivector[ j ] += 1;
      
      dummy = SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 );

      //      fprintf( stdout, "\t super_index_changed = %d\n", dummy );

      tables_extra.rho_extra_index_table_change1_plus1.imatrix[ i ][ j ] = dummy;

      dummy_rho_extra_index1.ivector[ j ] -= 1;

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_extra_index_table_change1_plus1\n" );

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc_extra-1; j++ ){

      fprintf( stdout, "%d, ", tables_extra.rho_extra_index_table_change1_plus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables_extra.rho_extra_index_table_change1_plus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_extra_index_table_change1_minus1 */

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc_extra; j++ ){
      
      dummy_rho_extra_index1.ivector[ j ] = tables_extra.rho_extra_index_table.imatrix[ i ][ j ];
      
    }
    
    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 ) );

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      dummy_rho_extra_index1.ivector[ j ] -= 1;
      
      dummy = SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 );

      //      fprintf( stdout, "\t super_index_changed = %d\n", dummy );

      tables_extra.rho_extra_index_table_change1_minus1.imatrix[ i ][ j ] = dummy;

      dummy_rho_extra_index1.ivector[ j ] += 1;

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_extra_index_table_change1_minus1\n" );

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc_extra-1; j++ ){

      fprintf( stdout, "%d, ", tables_extra.rho_extra_index_table_change1_minus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables_extra.rho_extra_index_table_change1_minus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_extra_index_table_change1_plus2 */

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc_extra; j++ ){
      
      dummy_rho_extra_index1.ivector[ j ] = tables_extra.rho_extra_index_table.imatrix[ i ][ j ];
      
    }
    
    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 ) );

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      dummy_rho_extra_index1.ivector[ j ] += 2;
      
      dummy = SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 );

      //      fprintf( stdout, "\t super_index_changed = %d\n", dummy );

      tables_extra.rho_extra_index_table_change1_plus2.imatrix[ i ][ j ] = dummy;

      dummy_rho_extra_index1.ivector[ j ] -= 2;

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_extra_index_table_change1_plus2\n" );

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc_extra-1; j++ ){

      fprintf( stdout, "%d, ", tables_extra.rho_extra_index_table_change1_plus2.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables_extra.rho_extra_index_table_change1_plus2.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_extra_index_table_change1_minus2 */

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc_extra; j++ ){
      
      dummy_rho_extra_index1.ivector[ j ] = tables_extra.rho_extra_index_table.imatrix[ i ][ j ];
      
    }
    
    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 ) );

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      dummy_rho_extra_index1.ivector[ j ] -= 2;
      
      dummy = SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 );

      //      fprintf( stdout, "\t super_index_changed = %d\n", dummy );

      tables_extra.rho_extra_index_table_change1_minus2.imatrix[ i ][ j ] = dummy;

      dummy_rho_extra_index1.ivector[ j ] += 2;

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_extra_index_table_change1_minus2\n" );

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc_extra-1; j++ ){

      fprintf( stdout, "%d, ", tables_extra.rho_extra_index_table_change1_minus2.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables_extra.rho_extra_index_table_change1_minus2.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_extra_index_table_change2_plus1_plus1 */

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc_extra; j++ ){
      
	dummy_rho_extra_index1.ivector[ j ] = tables_extra.rho_extra_index_table.imatrix[ i ][ j ];
      
    }

    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 ) );

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      for( k=0; k<N_coor_red_loc_extra; k++ ){

	dummy_rho_extra_index1.ivector[ j ] += 1;

	dummy_rho_extra_index1.ivector[ k ] += 1;
	
	dummy = SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 );
	
	//	fprintf( stdout, "\t super_index_changed = %d\n", dummy );
	
	tables_extra.rho_extra_index_table_change2_plus1_plus1.imatrix[ i ][ j +k *N_coor_red_loc_extra ] = dummy;

	dummy_rho_extra_index1.ivector[ k ] -= 1;

	dummy_rho_extra_index1.ivector[ j ] -= 1;

      }

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_extra_index_table_change2_plus1_plus1\n" );

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc_extra*N_coor_red_loc_extra -1; j++ ){

      fprintf( stdout, "%d, ", tables_extra.rho_extra_index_table_change2_plus1_plus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables_extra.rho_extra_index_table_change2_plus1_plus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_extra_index_table_change2_plus1_minus1 */

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc_extra; j++ ){
      
	dummy_rho_extra_index1.ivector[ j ] = tables_extra.rho_extra_index_table.imatrix[ i ][ j ];
      
    }

    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 ) );

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      for( k=0; k<N_coor_red_loc_extra; k++ ){

	dummy_rho_extra_index1.ivector[ j ] += 1;

	dummy_rho_extra_index1.ivector[ k ] -= 1;
	
	dummy = SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 );
	
	//	fprintf( stdout, "\t super_index_changed = %d\n", dummy );
	
	tables_extra.rho_extra_index_table_change2_plus1_minus1.imatrix[ i ][ j +k *N_coor_red_loc_extra ] = dummy;

	dummy_rho_extra_index1.ivector[ k ] += 1;

	dummy_rho_extra_index1.ivector[ j ] -= 1;

      }

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_extra_index_table_change2_plus1_minus1\n" );

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc_extra*N_coor_red_loc_extra -1; j++ ){

      fprintf( stdout, "%d, ", tables_extra.rho_extra_index_table_change2_plus1_minus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables_extra.rho_extra_index_table_change2_plus1_minus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* rho_extra_index_table_change2_minus1_minus1 */

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    for( j=0; j<N_coor_red_loc_extra; j++ ){
      
	dummy_rho_extra_index1.ivector[ j ] = tables_extra.rho_extra_index_table.imatrix[ i ][ j ];
      
    }

    //    fprintf( stdout, "\t super_index = %d\t", SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 ) );

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      for( k=0; k<N_coor_red_loc_extra; k++ ){

	dummy_rho_extra_index1.ivector[ j ] -= 1;

	dummy_rho_extra_index1.ivector[ k ] -= 1;
	
	dummy = SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_extra_index1 );
	
	//	fprintf( stdout, "\t super_index_changed = %d\n", dummy );
	
	tables_extra.rho_extra_index_table_change2_minus1_minus1.imatrix[ i ][ j +k *N_coor_red_loc_extra ] = dummy;

	dummy_rho_extra_index1.ivector[ k ] += 1;

	dummy_rho_extra_index1.ivector[ j ] += 1;

      }

    }

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "rho_extra_index_table_change2_minus1_minus1\n" );

  for( i=0; i<sqrt_max_rho_extra_index_loc; i++ ){

    fprintf( stdout, "%d \t( ", i );


    for( j=0; j<N_coor_red_loc_extra*N_coor_red_loc_extra -1; j++ ){

      fprintf( stdout, "%d, ", tables_extra.rho_extra_index_table_change2_minus1_minus1.imatrix[ i ][ j ] );

    }
    fprintf( stdout, "%d ", tables_extra.rho_extra_index_table_change2_minus1_minus1.imatrix[ i ][ j ] );
    fprintf( stdout, ")\n");

  }

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------

/* rho_extra_index_inverse */

int rho_extra_index_inverse( const int index, ivector_p index1_p, ivector_p index2_p ){

  /* dummies */
  int dummy_superindex1;
  int dummy_superindex2;
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){
    
    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );

  }

  if( !info ){

    if( IVECTOR_CONSINSTENCY( *index1_p, *index2_p ) ){

      fprintf( stderr, "ERROR: rho_extra_indexing inconsistency.\n");
      fflush( stderr );

      info=1;

    }

  }

#endif /* __DEBUG__ */


  if( !info ){

    if( MATRIX_ARG_INVERSE( index, dummy_superindex1, dummy_superindex2, sqrt_max_rho_extra_index_loc, sqrt_max_rho_extra_index_loc ) ) info=1;

    if( !info ){

      if( RHO_EXTRA_INDEX_INVERSE_AUX( dummy_superindex1, *index1_p ) ) info=1;  

      if( RHO_EXTRA_INDEX_INVERSE_AUX( dummy_superindex2, *index2_p ) ) info=1;  

    }

  }


  return info;

}
//------------------------------------------

/* rho_extra_index_inverse_aux */

int rho_extra_index_inverse_aux( const int index, ivector_p index_p ){

  /* dummies */
  int i;
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){
    
    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */


  if( !info ){
    
    for( i=0; i<N_coor_red_loc_extra; i++ ){
      
      index_p->ivector[ i ] = tables_extra.rho_extra_index_table.imatrix[ index ][ i ];
      
    }
    
  }


  return info;

}

//------------------------------------------

/* rho_index_aux_to_rho_extra_index_aux */

int rho_index_aux_to_rho_extra_index_aux( int rho_index_aux ){

  /* dummies */
  int rho_extra_index_aux;
  int info=0;


  if( RHO_INDEX_INVERSE_AUX( rho_index_aux, dummy_rho_index1 ) ) info=1;

  if( !info ){

    rho_extra_index_aux = SEARCH_RHO_EXTRA_INDEX_TABLE( dummy_rho_index1 );

  }
  else{

    rho_extra_index_aux = sqrt_max_rho_extra_index_loc;

  }


  return rho_extra_index_aux;

}

//------------------------------------------

/* compute_rho_index_aux_translate*/

int compute_rho_index_aux_translate( const constants constants ){

  /* constants */
  int sqrt_max_rho_index;
  /* dummies */
  int i;
  int info=0;


  sqrt_max_rho_index = constants.sqrt_max_rho_index;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){

    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );
    
  }

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "\n" );
  fprintf( stdout, "DOING: compute_rho_index_aux_translate\n" );
  fprintf( stdout, "\n" );
  fprintf( stdout, "rho_index_aux --> rho_extra_index_aux\n" );

#endif /* __DEBUG_PLUS__ */


  for( i=0; i<sqrt_max_rho_index; i++ ){

    rho_index_aux_translate.ivector[ i ] = RHO_INDEX_AUX_TO_RHO_EXTRA_INDEX_AUX( i );


#ifdef __DEBUG_PLUS__

    fprintf( stdout, "%d --> %d\n", i, RHO_INDEX_AUX_TO_RHO_EXTRA_INDEX_AUX( i ) );

#endif /* __DEBUG_PLUS__ */

  }

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------

/* check rho_extra_indexing */

int check_rho_extra_indexing( void ){

  /* dummies */
  ivector index_1, index_2;
  ivector index_1_tmp, index_2_tmp;
  int index, index_aux;
  int i, j;
  int sum;
  int info=0;


#ifdef __DEBUG__

  if( !rho_extra_indexing_allocate_flag ){
    
    info=1;
    
    fprintf( stderr, "ERROR: rho_extra_indexing not allocated yet.\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

//------------------------------------------
//                  RHO_EXTRA
//------------------------------------------

#ifdef __DEBUG__

  fprintf( stdout, "CHECKING RHO_EXTRA_INDEX\n");
  fflush( stdout );
  
#endif /* __DEBUG__ */

  
  fprintf( stdout, "\n" );
  fprintf( stdout, "----------\n" );
  fprintf( stdout, "RHO_EXTRA_INDEX_TABLE\n" );
  
  for( index=0; index<sqrt_max_rho_extra_index_loc; index++ ){
      
    fprintf( stdout, "%d \t( ", index );
    
    sum=0;
    for( j=0; j<N_coor_red_loc_extra; j++ ){
      
      fprintf( stdout, "%d ", tables_extra.rho_extra_index_table.imatrix[ index ][ j ] );

      sum += tables_extra.rho_extra_index_table.imatrix[ index ][ j ];

    }

    if( sum == CEID_order_loc_extra -2 ){
      
      fprintf( stdout, ")* \n" );
      
    }
    else if( sum == CEID_order_loc_extra -3 ){
      
      fprintf( stdout, ")** \n" );
      
	}
    else{
      
      fprintf( stdout, ") \n" );

    }

  }

  fprintf( stdout, "\n" );


  if( !info ){

    if( IVECTOR_ALLOCATE( N_coor_red_loc_extra, index_1 ) ) info=1;
    if( IVECTOR_ALLOCATE( N_coor_red_loc_extra, index_2 ) ) info=1;

  }

  fprintf( stdout, "\n" );
  fprintf( stdout, "----------\n" );
  fprintf( stdout, "RHO_EXTRA_INDICES\n" );
  
  for( index=0; index<max_rho_extra_index_loc; index++ ){
    
    fprintf( stdout, "%d \t( ", index );

    if( RHO_EXTRA_INDEX_INVERSE( index, index_1, index_2) ) info=1;

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      fprintf( stdout, "%d ", index_1.ivector[ j ] );

    }

    fprintf( stdout, ")\t( " );

    for( j=0; j<N_coor_red_loc_extra; j++ ){

      fprintf( stdout, "%d ", index_2.ivector[ j ] );

    }

    fprintf( stdout, ")\t\t" );

    index_aux = RHO_EXTRA_INDEX( index_2, index_1 ); 

    fprintf( stdout, "%d\n", index_aux );

  }

  fprintf( stdout, "----------\n");
  fprintf( stdout, "\n" );


  //  -------------------------------------------------
  //  -------------------------------------------------


  if( !info ){

    if( IVECTOR_ALLOCATE( N_coor_red_loc_extra, index_1_tmp ) ) info=1;
    if( IVECTOR_ALLOCATE( N_coor_red_loc_extra, index_2_tmp ) ) info=1;

  }

  if( !info ){

    if( IVECTOR_CONSINSTENCY( index_1, index_2 ) ) info=1;

  }

  for( index=0; index<max_rho_extra_index_loc; index++ ){

    /*
      fprintf( stdout, "->index:     %d\n", index );
      fflush( stdout );
    */

    if( RHO_EXTRA_INDEX_INVERSE( index, index_1, index_2) ) info=1;

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

    // check: RHO_EXTRA_INDEX_CHANGE1_PLUS1_FIRST
    //    fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE1_PLUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*    
	    fprintf( stdout, "-------------------------------------\n");
	    fprintf( stdout, "---> i=%d\n\n", i );
	    fflush( stdout );
      */

	index_aux = RHO_EXTRA_INDEX_CHANGE1_PLUS1_FIRST( index, i );

	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/

      if( index_aux < max_rho_extra_index_loc ){
	
	if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1_tmp\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2_tmp\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE1_MINUS1_FIRST( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  /*
	      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	      fflush( stdout );
	  */
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE1_MINUS1_FIRST
    // fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE1_MINUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_EXTRA_INDEX_CHANGE1_MINUS1_FIRST( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_extra_index_loc ){
	
	if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE1_PLUS1_FIRST( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  /*  
	      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	      fflush( stdout );
	  */
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE1_PLUS2_FIRST
    // fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE1_PLUS2_FIRST\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_EXTRA_INDEX_CHANGE1_PLUS2_FIRST( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_extra_index_loc ){
	
	if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE1_MINUS2_FIRST( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  /*  
	      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	      fflush( stdout );
	  */
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE1_MINUS2_FIRST
    // fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE1_MINUS2_FIRST\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_EXTRA_INDEX_CHANGE1_MINUS2_FIRST( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_extra_index_loc ){
	
	if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE1_PLUS2_FIRST( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  /*  
	      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	      fflush( stdout );
	  */
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_FIRST
    //    fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc_extra; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_FIRST( index, i, j );
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
       
	if( i != j && index_aux < max_rho_extra_index_loc ){
	
	  if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*	
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index_aux, i, j );
	  
	  /*
	    fprintf( stdout, "index_aux [check]: %d\n\n", index_aux );
	    fflush( stdout );
	  */

	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    /*	
		fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
		fflush( stdout );
	    */
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_FIRST
    //    fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc_extra; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	

	index_aux = RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( index, i, j );
	
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
	
	if( i != j && index_aux < max_rho_extra_index_loc ){
	
	  if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( index_aux, j, i ); //WARNING: note the index inversion
	  
	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    /*
		fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
		fflush( stdout );
	    */
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_FIRST
    // fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_FIRST\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc_extra; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	

	index_aux = RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( index, i, j );
	
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
	
	if( i != j && index_aux < max_rho_extra_index_loc ){
	
	  if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_FIRST( index_aux, i, j );
	  
	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    /*  
		fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
		fflush( stdout );
	    */
	    
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

    // check: RHO_EXTRA_INDEX_CHANGE1_PLUS1_SECOND
    //    fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE1_PLUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*    
	    fprintf( stdout, "-------------------------------------\n");
	    fprintf( stdout, "---> i=%d\n\n", i );
	    fflush( stdout );
      */

	index_aux = RHO_EXTRA_INDEX_CHANGE1_PLUS1_SECOND( index, i );

	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/

      if( index_aux < max_rho_extra_index_loc ){
	
	if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1_tmp\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2_tmp\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE1_MINUS1_SECOND( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  /*
	      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	      fflush( stdout );
	  */
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE1_MINUS1_SECOND
    // fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE1_MINUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_EXTRA_INDEX_CHANGE1_MINUS1_SECOND( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_extra_index_loc ){
	
	if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE1_PLUS1_SECOND( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  /*  
	      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	      fflush( stdout );
	  */
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE1_PLUS2_SECOND
    // fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE1_PLUS2_SECOND\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_EXTRA_INDEX_CHANGE1_PLUS2_SECOND( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_extra_index_loc ){
	
	if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE1_MINUS2_SECOND( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  /*  
	      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	      fflush( stdout );
	  */
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE1_MINUS2_SECOND
    // fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE1_MINUS2_SECOND\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      
      index_aux = RHO_EXTRA_INDEX_CHANGE1_MINUS2_SECOND( index, i );


      /*
	fprintf( stdout, "index_aux: %d\n\n", index_aux );
	fflush( stdout );
      */

      if( index_aux < max_rho_extra_index_loc ){
	
	if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	/*
	  fprintf( stdout, "print index_1\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	  
	  fprintf( stdout, "print index_2\n" ) ;
	  fflush( stdout );
	  if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE1_PLUS2_SECOND( index_aux, i );

	if( index_aux != index ){

	  fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	  fflush( stderr );
	  
	  /*  
	      fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
	      fflush( stdout );
	  */
	  
	  info=1;
	  
	  break;
	
	}

      }

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_SECOND
    //    fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc_extra; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	
	index_aux = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_SECOND( index, i, j );
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
       
	if( i != j && index_aux < max_rho_extra_index_loc ){
	
	  if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*	
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_SECOND( index_aux, i, j );
	  
	  /*
	    fprintf( stdout, "index_aux [check]: %d\n\n", index_aux );
	    fflush( stdout );
	  */

	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    /*	
		fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
		fflush( stdout );
	    */
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_SECOND
    //    fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc_extra; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	

	index_aux = RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_SECOND( index, i, j );
	
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
	
	if( i != j && index_aux < max_rho_extra_index_loc ){
	
	  if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_EXTRA_INDEX_CHANGE2_PLUS1_MINUS1_SECOND( index_aux, j, i ); //WARNING: note the index inversion
	  
	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    /*
		fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
		fflush( stdout );
	    */
	    
	    info=1;
	  
	    break;
	    
	  }

	}

      } /* j loop */
      
      if( info ) break;

    } /* i loop */

    /* ---------------------------------- */
    // check: RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_SECOND
    // fprintf( stdout, "----> RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_SECOND\n" );

    for( i=0; i<N_coor_red_loc_extra; i++ ){

      /*
	fprintf( stdout, "-------------------------------------\n");
	fprintf( stdout, "---> i=%d\n\n", i );
	fflush( stdout );
      */
      
      for( j=0; j<N_coor_red_loc_extra; j++ ){

	/*
	  fprintf( stdout, "-------------------------------------\n");
	  fprintf( stdout, "---> j=%d\n\n", j );
	  fflush( stdout );
	*/
	

	index_aux = RHO_EXTRA_INDEX_CHANGE2_MINUS1_MINUS1_SECOND( index, i, j );
	
	
	/*
	  fprintf( stdout, "index_aux: %d\n\n", index_aux );
	  fflush( stdout );
	*/
	
	if( i != j && index_aux < max_rho_extra_index_loc ){
	
	  if( RHO_EXTRA_INDEX_INVERSE( index_aux, index_1_tmp, index_2_tmp ) ) info=1;
	
	  /*
	    fprintf( stdout, "print index_1\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_1_tmp ) ) info=1;
	    
	    fprintf( stdout, "print index_2\n" ) ;
	    fflush( stdout );
	    if( IVECTOR_PRINT_PLUS( stdout, index_2_tmp ) ) info=1;
	  */
	  
	  index_aux = RHO_EXTRA_INDEX_CHANGE2_PLUS1_PLUS1_SECOND( index_aux, i, j );
	  
	  if( index_aux != index ){

	    fprintf( stderr, "CONSISTENCY ERROR for index %d\n", index );
	    fflush( stderr );
	  
	    /*  
		fprintf( stdout, "CONSISTENCY ERROR for index %d\n", index );
		fflush( stdout );
	    */
	    
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

  if( !info ){

    if( IVECTOR_FREE( index_1 ) ) info=1;
    if( IVECTOR_FREE( index_2 ) ) info=1;

    if( IVECTOR_FREE( index_1_tmp ) ) info=1;
    if( IVECTOR_FREE( index_2_tmp ) ) info=1;

  }

  return info;

}
