
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
#include "PolyCEID_CI_table.h"


/*********************
  FUNCTIONS & MACROS
*********************/


int CI_table_no_spin( constants_p constants_p ){

  /* constants */
  unsigned short int flag_no_spin_flip;
  int total_levels;
  int CAS_levels;
  int total_electrons;
  int CAS_electrons;
  imatrix_p single_level_occupation_p;
  /* dummy */
  imatrix table_tmp;
  imatrix table_tmp2;
  int frozen_electrons=0;
  int CAS_dim=0;
  int dim_space=0;
  int i, j;
  int info=0;

  
  flag_no_spin_flip           =  constants_p->flag_no_spin_flip;
  total_levels                =  constants_p->N_levels_single;
  CAS_levels                  =  constants_p->N_levels_single_CAS;
  total_electrons             =  constants_p->N_electrons;
  CAS_electrons               =  constants_p->N_electrons_CAS;
  single_level_occupation_p   = &(constants_p->single_level_occupation);


  if( CAS_electrons > total_electrons || total_electrons < 1 ){

    fprintf( stderr, "ERROR: invalid total_electrons [%d]\n", total_electrons ); 
    fflush( stderr );

    info=1;
 
  }

  if( total_electrons %2 || CAS_electrons %2 ){

    fprintf( stderr, "ERROR: only closed shells\n" );
    fflush( stderr );

    info=1;

  }


  if( !info ){

    CAS_dim = CHOOSE( CAS_levels, CAS_electrons /2 );
    
    //    fprintf( stdout, "CAS_dim = %d\n", CAS_dim );
    
    dim_space = CAS_dim *CAS_dim;
    
    //    fprintf( stdout, "dim_space = %d\n\n", dim_space );
    
    if( IMATRIX_ALLOCATE( CAS_dim, CAS_levels, table_tmp ) ) info=1;
  

    if( IMATRIX_ALLOCATE( dim_space, SPIN_DEG *total_levels, table_tmp2 ) ) info=1;
  

    if( CI_TABLE_CHOOSE( 0,0, CAS_levels, CAS_electrons /2, table_tmp ) ) info=1;
  
    /*  
	fprintf( stdout, "table_tmp:\n" );
	if( IMATRIX_PRINT_PLUS( stdout, table_tmp ) ) info=1;
	fprintf( stdout, "\n" ); 
    */

  }

  if( !info ){

    frozen_electrons = total_electrons -CAS_electrons;
    
    //    fprintf( stdout, "frozen_electrons = %d\n", frozen_electrons );


    if( frozen_electrons %2 ){
      
      fprintf( stderr, "ERROR: only an even number of frozen electrons are allowed\n" );
      fflush( stderr );
      
      info=1;
      
    }

    if( CAS_levels + ( frozen_electrons /2 ) > total_levels ){
      
      fprintf( stderr, "ERROR: invalid frozen_electrons [%d]\n", frozen_electrons  );
      fflush( stderr );
      
      info=1;
      
    }


  }

  if( !info ){

    if( IMATRIX_ZERO( table_tmp2 ) ) info=1;

    /* set frozen electrons */
    for( i=0; i<dim_space; i++ ){

      for( j=0; j<frozen_electrons; j++ ){


	table_tmp2.imatrix[ i ][ j ] = 1;


      } /* end loop j*/

    } /* end loop i */
 
  }

  if( !info ){

    /* make table*/
    for( i=0; i<dim_space; i++ ){

      for( j=0; j< CAS_levels; j++ ){


	table_tmp2.imatrix[ i ][ frozen_electrons +2 *j ] = 
	  table_tmp.imatrix[ i /CAS_dim ][ j ]; /* spin up */

	table_tmp2.imatrix[ i ][ frozen_electrons +2 *j +1 ] = 
	  table_tmp.imatrix[ i %CAS_dim ][ j ]; /* spin down */



      } /* end loop j*/

    } /* end loop i */

  }
  

#ifdef __DEBUG__
  
  fprintf( stdout, "table_tmp2:\n" );
  if( IMATRIX_PRINT_PLUS( stdout, table_tmp2 ) ) info=1;
  fprintf( stdout, "\n" ); 
  
#endif /* __DEBUG__ */


  if( flag_no_spin_flip ){

    if( PURGE_SPIN_FLIP( table_tmp2, *single_level_occupation_p ) ) info=1;

  }
  else{

    if( IMATRIX_COPY( *single_level_occupation_p, table_tmp2 ) ) info=1;

  }


  /* energy pruning */
  if( ENERGY_PRUNING( *constants_p ) ) info=1; 

#ifdef __DEBUG__
  
  fprintf( stdout, "single_level_occupation:\n" );
  if( IMATRIX_PRINT_PLUS( stdout, *single_level_occupation_p ) ) info=1;
  fprintf( stdout, "\n" ); 
  fflush( stdout );

#endif /* __DEBUG__ */
  

  if( IMATRIX_FREE( table_tmp2 ) ) info=1;
  
  if( IMATRIX_FREE( table_tmp ) ) info=1;
  

  return info;

}

//------------------------------------------

int CI_table_spin( constants_p constants_p ){

  /* constants */
  unsigned short int flag_no_spin_flip;
  int total_levels;
  int CAS_levels;
  int total_electrons;
  int CAS_electrons;
  imatrix_p single_level_occupation_p;
  /* dummy */
  imatrix table_tmp;
  imatrix table_tmp2;
  int frozen_electrons=0;
  int dim_space=0;
  int i, j;
  int info=0;

  
  flag_no_spin_flip           =  constants_p->flag_no_spin_flip;
  total_levels                =  constants_p->N_levels_single;
  CAS_levels                  =  constants_p->N_levels_single_CAS;
  total_electrons             =  constants_p->N_electrons;
  CAS_electrons               =  constants_p->N_electrons_CAS;
  single_level_occupation_p   = &(constants_p->single_level_occupation);


  if( CAS_electrons > total_electrons || total_electrons < 1 ){

    fprintf( stderr, "ERROR: invalid total_electrons [%d]\n", total_electrons ); 
    fflush( stderr );

    info=1;
 
  }


  if( !info ){

    dim_space = CHOOSE( SPIN_DEG *CAS_levels, CAS_electrons );
    
    //    fprintf( stdout, "dim_space = %d\n\n", dim_space );
    

    if( IMATRIX_ALLOCATE( dim_space, SPIN_DEG *CAS_levels,   table_tmp  ) ) info=1;
  
    if( IMATRIX_ALLOCATE( dim_space, SPIN_DEG *total_levels, table_tmp2 ) ) info=1;
  

    if( CI_TABLE_CHOOSE( 0,0, SPIN_DEG *CAS_levels, CAS_electrons, table_tmp ) ) info=1;


#ifdef __DEBUG_PLUS__
  
    fprintf( stdout, "table_tmp:\n" );
    if( IMATRIX_PRINT_PLUS( stdout, table_tmp ) ) info=1;
    fprintf( stdout, "\n" ); 

#endif /* __DEBUG_PLUS__ */


  }

  if( !info ){

    frozen_electrons = total_electrons -CAS_electrons;
    
    //    fprintf( stdout, "frozen_electrons = %d\n", frozen_electrons );


    if( SPIN_DEG *CAS_levels +frozen_electrons > SPIN_DEG *total_levels ){
      
      fprintf( stderr, "ERROR: invalid frozen_electrons [%d]\n", frozen_electrons  );
      fflush( stderr );
      
      info=1;
      
    }


  }

  if( !info ){

    if( IMATRIX_ZERO( table_tmp2 ) ) info=1;

    /* set frozen electrons */
    for( i=0; i<dim_space; i++ ){

      for( j=0; j<frozen_electrons; j++ ){


	table_tmp2.imatrix[ i ][ j ] = 1;


      } /* end loop j*/

    } /* end loop i */
 
  }

  if( !info ){

    /* make table*/
    for( i=0; i<dim_space; i++ ){

      for( j=0; j< SPIN_DEG *CAS_levels; j++ ){


	table_tmp2.imatrix[ i ][ frozen_electrons +j ] = 
	  table_tmp.imatrix[ i ][ j ];


      } /* end loop j*/

    } /* end loop i */

  }
  

#ifdef __DEBUG__
  
  fprintf( stdout, "table_tmp2:\n" );
  if( IMATRIX_PRINT_PLUS( stdout, table_tmp2 ) ) info=1;
  fprintf( stdout, "\n" ); 
  
#endif /* __DEBUG__ */

  /*
    fprintf( stdout, "single_level_occupation [before]:\n" );
    if( IMATRIX_PRINT_PLUS( stdout, *single_level_occupation_p ) ) info=1;
    fprintf( stdout, "\n" ); 
    fflush( stdout );
  */

  if( flag_no_spin_flip ){

    if( PURGE_SPIN_FLIP( table_tmp2, *single_level_occupation_p ) ) info=1;

  }
  else{

    if( IMATRIX_COPY( *single_level_occupation_p, table_tmp2 ) ) info=1;

  }


#ifdef __DEBUG__
  
  fprintf( stdout, "single_level_occupation [after purge_spin_flip]:\n" );
  if( IMATRIX_PRINT_PLUS( stdout, *single_level_occupation_p ) ) info=1;
  fprintf( stdout, "\n" ); 
  fflush( stdout );

#endif /* __DEBUG__ */

  
  /* energy pruning */
  if( ENERGY_PRUNING( *constants_p ) ) info=1; 


#ifdef __DEBUG__
  
  fprintf( stdout, "single_level_occupation [after energy_pruning]:\n" );
  if( IMATRIX_PRINT_PLUS( stdout, *single_level_occupation_p ) ) info=1;
  fprintf( stdout, "\n" ); 
  fflush( stdout );

#endif /* __DEBUG__ */
  

  if( IMATRIX_FREE( table_tmp2 ) ) info=1;
  
  if( IMATRIX_FREE( table_tmp ) ) info=1;
  

  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------


int CI_table_choose( const int i, const int j, const int CAS_levels, const int CAS_electrons, imatrix_p table_p  ){

  /* dummies */
  int length;
  int index;
  int k;
  int coeff;
  static int counter=0;
  static int info=0;


  //  fprintf( stdout, "counter = %d\n", counter );

  if( 0 == counter ){

    if( IMATRIX_ZERO( *table_p ) ) info=1;

  }

  counter++;

  /*
    fprintf( stdout, "i = %d, j = %d\n", i, j );

    fprintf( stdout, "CAS_levels = %d\n", CAS_levels );
    
    fprintf( stdout, "CAS_electrons = %d\n\n", CAS_electrons );
  */

  /* Checks */
  if( info || CAS_electrons>CAS_levels ){

    fprintf( stderr, "ERROR: invalid condition [1]\n" );
    fflush( stderr );

    info=1;

  }

  if( i<0 || i>table_p->imatrix_dim_row  ){

    fprintf( stderr, "ERROR: i\n" );
    fflush( stderr );

    info=1;

  }

  if( j<0 || j>table_p->imatrix_dim_column  ){

    fprintf( stderr, "ERROR: j\n" );
    fflush( stderr );

    info=1;

  }


  /* routine */
  if( CAS_electrons>1 ){

    length=0;

    //    fprintf( stdout, "#-------------------------\n" );

    for( k=0; k<( CAS_levels -CAS_electrons +1 ); k++ ){

      /*
	fprintf( stdout, "--> k = %d\n", k );

	fprintf( stdout, "length = %d\n", length );
      */

      coeff = CHOOSE( CAS_levels-k-1, CAS_electrons-1 );

      //      fprintf( stdout, "coeff = %d\n", coeff );


      //      fprintf( stdout, "i+length = %d, j+k = %d\n", i+length, j+k );
      for( index=0; index<coeff; index++ ){

	table_p->imatrix[ i+length+index ][ j+k ]=1;

      } /* end index loop */

      /*
	fprintf( stdout, "CI_table: internal\n" );
	if( IMATRIX_PRINT_PLUS( stdout, *table_p ) ) info=1;
	fprintf( stdout, "\n" ); 

	fprintf( stdout, "CAS_levels-1 = %d, CAS_electrons-1 = %d\n", CAS_levels-1, CAS_electrons-1);
	fflush( stdout );
      */

      if( CI_TABLE_CHOOSE( i+length, j+k+1, CAS_levels-k-1, CAS_electrons-1, *table_p  ) ) info=1; /* WARNING: recursion */
      
      length += coeff;

    } /* end k loop */

    //    fprintf( stdout, "#-------------------------\n" );

  }
  else if( 1==CAS_electrons ){


    for( index=0; index<CAS_levels; index++ ){
      
      table_p->imatrix[ i+index ][ j+index ]=1;

    } /* end index loop */


  }
  else{

    fprintf( stderr, "ERROR: invalid condition [2]\n" );
    fflush( stderr );

    info=1;

  }

  return info;

}

//------------------------------------------

int purge_spin_flip( const imatrix imat_in, imatrix_p imat_out_p ){

  /* dummies */
  ivector ivec1;
  ivector ivec2;
  ivector purged_ivector;
  int  dim_row_in;
  int  dim_row_out;
  int  dim_column;
  int  i, j;
  int  h, k;
  int  counter=0;
  int  info=0;


  dim_row_in  = imat_in.imatrix_dim_row;

  dim_row_out = imat_out_p->imatrix_dim_row;

  dim_column = imat_in.imatrix_dim_column;


  if( dim_row_out > dim_row_in || imat_out_p->imatrix_dim_column != dim_column ){

    fprintf( stderr, "ERROR: incompatible matrices\n" );
    fflush( stderr );

    info=1;

  }


  if( !info ){

    if( IVECTOR_ALLOCATE( dim_column, ivec1 ) ) info=1;

    if( IVECTOR_ALLOCATE( dim_column, ivec2 ) ) info=1;

    if( IVECTOR_ALLOCATE( dim_row_in, purged_ivector ) ) info=1;

  }


  for( i=0; i<dim_row_in; i++ ){


    /* ivector copy */
    for( k=0; k<dim_column; k++ ){

      ivec1.ivector[ k ] = imat_in.imatrix[ i ][ k ];

    }

    /*
      fprintf( stdout, "ivec1:\n" );
      IVECTOR_PRINT_PLUS( stdout, ivec1 );
      fprintf( stdout, "\n" );
    */
    
    for( j=(i+1); j<dim_row_in; j++ ){

      if( counter >= dim_row_out ){

        fprintf( stderr, "ERROR: the counter has exceeded the number of rows [%d]\n", counter);
      
        info=1;

        break;

      }


      /* ivector copy */
      for( k=0; k<(dim_column)/SPIN_DEG; k++ ){

	for( h=0; h<SPIN_DEG; h++ ){

	  ivec2.ivector[ k *SPIN_DEG + h  ] = imat_in.imatrix[ j ][ k *SPIN_DEG +SPIN_DEG -1 -h  ];

	} /* end h loop */

      } /* end k loop */

      /*
        fprintf( stdout, "ivec2:\n" );
        IVECTOR_PRINT_PLUS( stdout, ivec2 );
        fprintf( stdout, "\n" );
      */

      if( !IVECTOR_COMPARE_MUTE( ivec1, ivec2 ) ){

	purged_ivector.ivector[ j ] = 1;

	counter++;

	break;

      }


    } /* end j loop */

  } /* end i loop */


#ifdef __DEBUG__

  fprintf( stdout, "purged_ivector:\n" );
  if( IVECTOR_PRINT_PLUS( stdout, purged_ivector ) ) info=1;
  fprintf( stdout, "\n" ); 

#endif /* __DEBUG__ */


  if( counter +dim_row_out < dim_row_in ){

    fprintf( stderr, "ERROR: the counter is smaller than the number of rows [%d]\n", counter);

    info=1;
    
  }


  if( !info ){

    /* final_copy */
    counter=0;

    for( i=0; i<dim_row_in; i++ ){

      if( !(purged_ivector.ivector[ i ] ) ){

	for( k=0; k<dim_column; k++ ){

	  imat_out_p->imatrix[ counter ][ k ]= imat_in.imatrix[ i ][ k ];
	  
	} /* end k loop */

	counter++;

      }

    } /* end i loop */

  }


  if( IVECTOR_FREE( purged_ivector ) ) info=1;

  if( IVECTOR_FREE( ivec2 ) ) info=1;

  if( IVECTOR_FREE( ivec1 ) ) info=1;


  return info;

}

//------------------------------------------

int energy_pruning( constants_p constants_p ){

  /* constanst */
  int        N_levels_single;
  int*       N_levels_many_p;
  imatrix    single_level_occupation;
  rvector    initial_positions;
  double     cutoff_energy;
  /* state */
  double     sum;
  double     E0=0.0e0;
  rvector    energy_rvector;
  ivector    sorting_ivector;
  imatrix    pruned_imatrix;
  rvector    eigenvalues;
  matrix     eigenvectors;
  matrix     dummy_matrix_single;
  int        pp[ NP_CONSTANTS ];
  int        level1;
  int        counter;
  int        i, k;
  int        info=0;


  N_levels_single               =   constants_p->N_levels_single;
  N_levels_many_p               = &(constants_p->N_levels_many);
  single_level_occupation       =   constants_p->single_level_occupation;
  initial_positions             =   constants_p->initial_atoms.positions;
  cutoff_energy                 =   constants_p->cutoff_energy;


  //  fprintf( stdout, "N_levels_many = %d\n", *N_levels_many_p );


  if( !info ){

    if( RVECTOR_ALLOCATE( *N_levels_many_p, energy_rvector ) )                     info=1;

    if( IVECTOR_ALLOCATE( *N_levels_many_p, sorting_ivector ) )                    info=1;

    if( RVECTOR_ALLOCATE( N_levels_single, eigenvalues ) )                         info=1;

    if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, eigenvectors ) )        info=1;

    if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, dummy_matrix_single ) ) info=1;

  }


  /* initialise hamiltonian */
  if( INITIALISE_HAMILTONIAN_SINGLE( *constants_p ) ) info=1;


  /* allocate indexing */
  if( INDEXING_ALLOCATE( *constants_p, 1 ) ) info=1; 


  /* H_matrix_single*/
  if( H_MATRIX_SINGLE_CHECK( *constants_p, initial_positions, dummy_matrix_single ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "H_matrix_single\n" );
  if( MATRIX_PRINT_PLUS( stdout, dummy_matrix_single ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* diagonalisation */
  if( DIAGONALISATION( dummy_matrix_single, eigenvectors, eigenvalues ) ) info=1;


#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "eigenvectors\n" );
  if( MATRIX_PRINT_PLUS( stdout, eigenvectors ) ) info=1;
  fprintf( stdout, "\n" );
  
  fprintf( stdout, "eigenvalues\n" );
  if( RVECTOR_PRINT_PLUS( stdout, eigenvalues ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* compute many-body energies */
  for( i=0; i<*N_levels_many_p; i++ ){

    /* set dummy to zero*/
    sum = 0.0e0;

    for( k=0; k<N_levels_single*SPIN_DEG; k++ ){

      //      fprintf( stdout, "single_level %d, single_level_occupation %d\n", k, single_level_occupation.imatrix[i][ k ] );

      if( single_level_occupation.imatrix[i][ k ] ){

	level1 = k /SPIN_DEG;

	//	fprintf( stdout, "OK, level1 = %d\n", level1 );

	/* sum update*/
	sum += eigenvalues.rvector[ level1 ];
	  
      } /* end conditional */

    } /* end k loop */

    //    fprintf( stdout, "energy[ %d ] = %le\n", i, sum );

    energy_rvector.rvector[ i ] = sum; 

  }


  // sorting by energy
  if( RVECTOR_QSORT( energy_rvector, sorting_ivector ) ) info=1;



  /* ground state energy */

  E0 = energy_rvector.rvector[ sorting_ivector.ivector[ 0 ] ]; 
  
  /*
    fprintf( stdout, "# %d E0 = %le\n", sorting_ivector.ivector[ 0 ], E0 );
    fflush( stdout );
  */

  /* set counter to zero */
  counter=0;

  /* energy pruning loop */
  for( i=0; i<*N_levels_many_p; i++ ){


    /* writing */
    if(  energy_rvector.rvector[ sorting_ivector.ivector[ i ] ] <= E0 +cutoff_energy ){

      counter++;

    }
    else{

      break;

    }

  } /* end i loop */


  //  fprintf( stdout, "counter = %d\n", counter );


  if( counter < *N_levels_many_p ){

    /* update N_levels_many */
    *N_levels_many_p = counter;

    /* set flag_pruning*/
    constants_p->flag_pruning=1;

  }


  if(  *N_levels_many_p < 1 ){

    fprintf( stderr, "ERROR: cutoff_energy too small [%le]\n", cutoff_energy );
    fflush( stderr );

    info=1;

  }


  if( !info ){

    if( IMATRIX_ALLOCATE( *N_levels_many_p, SPIN_DEG *N_levels_single, pruned_imatrix ) ) info=1;

  }

      
  // reallocate single_level_occupation;
  if( !info ){
    
    /* final_copy */
    for( i=0; i<*N_levels_many_p; i++ ){

      for( k=0; k<(SPIN_DEG *N_levels_single); k++ ){

	pruned_imatrix.imatrix[ i ][ k ] = constants_p->single_level_occupation.imatrix[ sorting_ivector.ivector[ i ] ][ k ];
	      
      } /* end k loop */

    } /* end i loop */

  }      


      
#ifdef __DEBUG_PLUS__
      
  fprintf( stdout, "single_level_occupation\n" );
  if( IMATRIX_PRINT_PLUS( stdout, constants_p->single_level_occupation ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "pruned_imatrix\n" );
  if( IMATRIX_PRINT_PLUS( stdout, pruned_imatrix ) ) info=1;
  fprintf( stdout, "\n" );
      
#endif /* __DEBUG_PLUS__ */
  

  if( constants_p->flag_pruning ){    
	
    /* reallocate single_level_occupation */
    if( IMATRIX_FREE( constants_p->single_level_occupation) ) info=1;
    
    if( IMATRIX_ALLOCATE( *N_levels_many_p, SPIN_DEG *N_levels_single, constants_p->single_level_occupation ) ) info=1;
	
  }

  /* final copy */
  if( IMATRIX_COPY( constants_p->single_level_occupation, pruned_imatrix ) ) info=1;
	
	
  if( constants_p->flag_pruning ){    

    // reallocate initial_many_body_state
    if( RVECTOR_FREE( constants_p->initial_many_body_state ) ) info=1;
  
    if( RVECTOR_ALLOCATE( *N_levels_many_p, constants_p->initial_many_body_state ) ) info=1;


    // reallocate excited_many_body_state
    if( RVECTOR_FREE( constants_p->excited_many_body_state ) ) info=1;

    if( RVECTOR_ALLOCATE( *N_levels_many_p, constants_p->excited_many_body_state ) ) info=1;


  }


  if( !info ){

    if( IMATRIX_FREE( pruned_imatrix ) ) info=1;

  }


  if( constants_p->flag_pruning ){    

    // reallocate constants

    /* single_level_occupation_reduced */
    if( IMATRIX_FREE( constants_p->single_level_occupation_reduced ) ) info=1;
    
    if( IMATRIX_ALLOCATE( *N_levels_many_p, N_levels_single, constants_p->single_level_occupation_reduced ) ) info=1;
    

    /* symmetry_multiplicity  */
    if( IMATRIX_FREE( constants_p->symmetry_multiplicity ) ) info=1;

    if( IMATRIX_ALLOCATE( *N_levels_many_p, 2*N_SYMMETRIES +1, constants_p->symmetry_multiplicity ) ) info=1;


  }

  /* constants allocation */
  pp[0] = constants_p->N_atoms;
  pp[1] = constants_p->N_levels_single;
  pp[2] = constants_p->N_levels_many;
  pp[3] = constants_p->N_chain;
  
  if( CONSTANTS_ALLOCATE( NP_CONSTANTS, pp, *constants_p, 1 ) ) info=1; //WARNING: initial allocation only many_body_hybridisation_table and symmetry_coefficients


  /* indexing free */
  if( INDEXING_FREE( 1 ) ) info=1;


  if( !info ){

    if( MATRIX_FREE( dummy_matrix_single ) ) info=1;

    if( MATRIX_FREE( eigenvectors ) )        info=1;

    if( RVECTOR_FREE( eigenvalues ) )        info=1;
 
    if( IVECTOR_FREE( sorting_ivector ) )    info=1;

    if( RVECTOR_FREE( energy_rvector ) )     info=1;

  }


  return info;

}

//------------------------------------------
