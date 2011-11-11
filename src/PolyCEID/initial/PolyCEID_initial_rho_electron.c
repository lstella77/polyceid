
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
#include "PolyCEID_initial_rho_electron.h"


/* private varaiables */
#ifdef __DEBUG__

unsigned short int compute_initial_rho_electron_flag=0;
unsigned short int compute_excited_rho_electron_flag=0;
unsigned short int workspace_initial_rho_electron_flag=0;

#endif /* __DEBUG__ */

matrix  H_tmp_eigenvectors;
rvector H_tmp_eigenvalues;


/*********************
  FUNCTIONS & MACROS
*********************/

/* compute_initial_rho_electron */

int compute_initial_rho_electron( const constants constants, state_p state_p ){

  /* constants */
  int              N_levels_many;
  rvector          initial_many_body_state;
  /* state */
  matrix_p         initial_rho_electron_p;
  /* dummies */
  int              i,j;
  int              index;
  int              info=0;


  N_levels_many               =  constants.N_levels_many;
  initial_many_body_state     =  constants.initial_many_body_state;

  initial_rho_electron_p      = &(state_p->initial_rho_electron);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_initial_rho_electron\n");

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    for( i=0; i<N_levels_many; i++ ){
     
      for( j=0; j<N_levels_many; j++ ){
 
	index   = ELECTRON_MANY_INDEX( i, j );
	  
	initial_rho_electron_p->matrix[ index ] = CMPLX( ( initial_many_body_state.rvector[ i ] ) *( initial_many_body_state.rvector[ j ] ) );

      } /* end j loop */

    } /* end i loop */

  }
	
    
#ifdef __DEBUG_PLUS__

  fprintf( stdout, "\n" );
  fprintf( stdout, "initial_rho_electron\n" );
  if( MATRIX_PRINT_PLUS( stdout, *initial_rho_electron_p )) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__*/


  return info;

}

//------------------------------------------

/* compute_excited_rho_electron */

int compute_excited_rho_electron( const constants constants, state_p state_p ){

  /* constants */
  int              N_levels_many;
  rvector          excited_many_body_state;
  /* state */
  matrix_p         initial_rho_electron_p;
  /* dummies */
  int              i,j;
  int              index;
  int              info=0;


  N_levels_many               =  constants.N_levels_many;
  excited_many_body_state     =  constants.excited_many_body_state;

  initial_rho_electron_p      = &(state_p->initial_rho_electron);


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: compute_excited_rho_electron\n");

#endif /* __DEBUG_PLUS__ */


  if( !info ){

    for( i=0; i<N_levels_many; i++ ){
     
      for( j=0; j<N_levels_many; j++ ){
 
	index   = ELECTRON_MANY_INDEX( i, j );
	  
	initial_rho_electron_p->matrix[ index ] = CMPLX( ( excited_many_body_state.rvector[ i ] ) *( excited_many_body_state.rvector[ j ] ) );

      } /* end j loop */

    } /* end i loop */

  }
	
    
#ifdef __DEBUG_PLUS__

  fprintf( stdout, "\n" );
  fprintf( stdout, "excited_rho_electron\n" );
  if( MATRIX_PRINT_PLUS( stdout, *initial_rho_electron_p )) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__*/


  return info;

}

//------------------------------------------

int allocate_symmetry_multiplicity( constants_p constants_p ){

  /* constants */
  unsigned short int  flag_no_spin_flip;
  int                 N_levels_single;
  int                 N_levels_many;
  imatrix             single_level_occupation;
  imatrix             single_level_occupation_reduced;
  imatrix_p           symmetry_multiplicity_p;
  /* dummies */
  int                 i;
  int                 k;
  int                 mult;
  int                 j;
  int                 mult_sign;
  int                 excitation_order;
  int                 info=0;


  flag_no_spin_flip               =  constants_p->flag_no_spin_flip;
  N_levels_single                 =  constants_p->N_levels_single;
  N_levels_many                   =  constants_p->N_levels_many;
  single_level_occupation         =  constants_p->single_level_occupation;
  single_level_occupation_reduced =  constants_p->single_level_occupation_reduced;
  symmetry_multiplicity_p         = &(constants_p->symmetry_multiplicity);


#ifdef __DEBUG__

  fprintf( stdout, "DOING: allocate_symmetry_multiplicity\n" );

#endif /* __DEBUG__ */

  
  for( i=0; i<N_levels_many; i++ ){

    // spin_multiplicity
    
    if( flag_no_spin_flip ){

      /* check spin symmetry */
      mult=0;

      for( j=0; j<N_levels_single; j++ ){

        for( k=0; k<SPIN_DEG; k++ ){

	  if( single_level_occupation.imatrix[i][ j*SPIN_DEG +k ] != single_level_occupation.imatrix[i][ j*SPIN_DEG +SPIN_DEG -1 -k ] ){

	    mult=1;

	    break;

	  }

        } /* end k loop */

        if( mult ) break;

      } /* end j loop */

      if( mult ){

        mult = 2;

      }
      else{

        mult = 1;
      
      }

      symmetry_multiplicity_p->imatrix[ i ][ 0 ] = mult;

    }
    else{

      symmetry_multiplicity_p->imatrix[ i ][ 0 ] = 1;

    }

    // spin_sign
    mult_sign=0;

    for( j=0; j<N_levels_single; j++ ){

      if( 2 == single_level_occupation_reduced.imatrix[ i ][ j ] ){

	mult_sign++;

      }

    } /* end j loop */

    if( mult_sign %2  ){

      // is odd
      mult_sign = -1;

    }
    else{

      // is even
      mult_sign = +1;
      
    }

    symmetry_multiplicity_p->imatrix[ i ][ 1 ] = mult_sign;
  

    /* excitation order */
    excitation_order=0;

    for( j=0; j<(SPIN_DEG *N_levels_single); j++ ){

      if( single_level_occupation.imatrix[ i ][ j ] != 
	  single_level_occupation.imatrix[ 0 ][ j ] ){

        excitation_order++;

      }

    } /* end j loop */
                                                           
    if( excitation_order %2 ){

      // is odd
      fprintf( stderr, "ERROR: illegal excitation_order [ %d/2 ]\n", excitation_order );
      fflush( stderr );

      info=1;

    }
    else{

      // is even
      excitation_order /= 2;

      symmetry_multiplicity_p->imatrix[ i ][ 2 *N_SYMMETRIES ] = excitation_order;
  
    }


  } /* end i loop */


  // global phase
  if( symmetry_multiplicity_p->imatrix[ 0 ][ 1 ] < 0 ){

    for( i=0; i<N_levels_many; i++ ){

      symmetry_multiplicity_p->imatrix[ i ][ 1 ] *= -1; //WARNING: change the sign

    } /* end i loop */

  }


  return info;

}

//------------------------------------------

int allocate_many_body_hybridisation_table( constants_p constants_p ){

  /* constants */
  int        N_levels_single;
  int        N_levels_many;
  imatrix    single_level_occupation;
  imatrix_p  many_body_hybridisation_table_p; 
  vector_p   symmetry_coefficient_p;
  imatrix    symmetry_multiplicity;
  /* dummies */
  ivector    ivec1, ivec2;
  int        flag_rec;
  int        counter;
  int        counter_success;
  int        sign;
  int        global_sign;
  int        i, j, k, l;
  int        orbital1, orbital2;
  int        orbital1_rec, orbital2_rec;
  int        spin1, spin2;
  int        spin1_rec, spin2_rec;
  int        info=0;


  N_levels_single                 =  constants_p->N_levels_single;
  N_levels_many                   =  constants_p->N_levels_many;
  single_level_occupation         =  constants_p->single_level_occupation;
  many_body_hybridisation_table_p = &(constants_p->many_body_hybridisation_table);
  symmetry_coefficient_p          = &(constants_p->symmetry_coefficient);
  symmetry_multiplicity           =  constants_p->symmetry_multiplicity;


  /* allocation */
  if( IVECTOR_ALLOCATE( SPIN_DEG *N_levels_single, ivec1 ) ) info=1;

  if( IVECTOR_ALLOCATE( SPIN_DEG *N_levels_single, ivec2 ) ) info=1;


#ifdef __DEBUG__

  fprintf( stdout, "DOING: allocate_many_body_hybridisation_table\n" );

#endif /* __DEBUG__ */


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "single_level_occupation [check]\n" );
  if( IMATRIX_PRINT_PLUS( stdout, single_level_occupation ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */ 


  /* set counter to zero */
  counter=0;

  for( i=0; i<N_levels_many; i++ ){

    for( j=(i+1); j<N_levels_many; j++ ){


    // fprintf( stdout, "--> indeces: (%d,%d)\n\n", i, j );

      /* default values */
      orbital2 = orbital1 = orbital2_rec = orbital1_rec = N_levels_single;
      spin2    = spin1    = spin2_rec    = spin1_rec    = SPIN_DEG;
      sign     = global_sign = 0;

      flag_rec=0;
      counter_success=0;


      // original

      /* initial copy */
      for( k=0; k<( SPIN_DEG *N_levels_single ); k++ ){
	  
	  ivec1.ivector[ k ] = single_level_occupation.imatrix[ i ][ k ];

	  ivec2.ivector[ k ] = single_level_occupation.imatrix[ j ][ k ];

      } /* end k loop */

      
      if( ALLOCATE_MANY_BODY_HYBRIDISATION_TABLE_AUX( *constants_p, ivec1, ivec2, orbital1, spin1, orbital2, spin2, sign ) ) info=1;

      // fprintf( stdout, "%d, %d, %d, %d, %d\n", orbital1, spin1, orbital2, spin2, sign  );

      if( orbital1 != N_levels_single && spin1 != SPIN_DEG &&
	  orbital2 != N_levels_single && spin2 != SPIN_DEG &&
          spin1 == spin2 ){

	orbital1_rec = orbital1;
	spin1_rec    = spin1;
	
	orbital2_rec = orbital2;
	spin2_rec    = spin2;

	flag_rec=1;

      }

      if( sign != 0 && spin1 == spin2 ){

	counter_success++;

	global_sign += sign;

        /*	
          fprintf( stdout, "counter_success = %d\n", counter_success );
	  fprintf( stdout, "global_sign     = %d\n", global_sign );
	*/

      }


      /***************
       * spin_flip j *
       ***************/
      if( symmetry_multiplicity.imatrix[j][0] > 1 ){
	

      // fprintf( stdout, "spin_flip [%d]\n", j );

	/* initial copy */
	for( k=0; k<N_levels_single; k++ ){
	  
	  for( l=0; l<SPIN_DEG; l++ ){
	  
	    ivec2.ivector[ k*SPIN_DEG +l ] = single_level_occupation.imatrix[ j ][ k*SPIN_DEG +SPIN_DEG -1 -l ];

	  } /* end l loop */
	  
	} /* end k loop */


	if( ALLOCATE_MANY_BODY_HYBRIDISATION_TABLE_AUX( *constants_p, ivec1, ivec2, orbital1, spin1, orbital2, spin2, sign ) ) info=1;

	// fprintf( stdout, "%d, %d, %d, %d, %d\n", orbital1, spin1, orbital2, spin2, sign  );

	if( orbital1 != N_levels_single && spin1 != SPIN_DEG &&
	    orbital2 != N_levels_single && spin2 != SPIN_DEG &&
            spin1 == spin2 ){
	  
	  if( flag_rec ){

	    if( orbital1_rec != orbital1 || orbital2_rec != orbital2 ){
	      
	      fprintf( stderr, "ERROR: spin-orbital inconsistency [%d,%d]\n", i, j );
	      fflush( stderr );
	      
	      info=1;
	      
	      break;
	      
	    }

	  }
	  else{

	    orbital1_rec = orbital1;
	    spin1_rec    = spin1;
	    
	    orbital2_rec = orbital2;
	    spin2_rec    = spin2;
	    
	    flag_rec=1;
	      
	  }

	}

	if( sign != 0 && spin1 == spin2 ){

	  counter_success++;

	  global_sign += ( symmetry_multiplicity.imatrix[ j ][ 1 ] ) *sign;

          /*	
            fprintf( stdout, "counter_success = %d\n", counter_success );
	    fprintf( stdout, "global_sign     = %d\n", global_sign );
	  */

	}

      } /* end conditional flip j*/


      /***************
       * spin_flip i *
       ***************/
      if( symmetry_multiplicity.imatrix[i][0] > 1 ){
	

	// fprintf( stdout, "spin_flip [%d]\n", i );

	/* initial copy */
	for( k=0; k<N_levels_single; k++ ){
	  
	  for( l=0; l<SPIN_DEG; l++ ){
	  
	    ivec1.ivector[ k*SPIN_DEG +l ]  = single_level_occupation.imatrix[ i ][ k*SPIN_DEG +SPIN_DEG -1 -l ];

	    ivec2.ivector[ SPIN_DEG *k +l ] = single_level_occupation.imatrix[ j ][ SPIN_DEG *k +l ];

	  } /* end l loop */
	  
	} /* end k loop */

	if( ALLOCATE_MANY_BODY_HYBRIDISATION_TABLE_AUX( *constants_p, ivec1, ivec2, orbital1, spin1, orbital2, spin2, sign ) ) info=1;

	// fprintf( stdout, "%d, %d, %d, %d, %d\n", orbital1, spin1, orbital2, spin2, sign  );

	if( orbital1 != N_levels_single && spin1 != SPIN_DEG &&
	    orbital2 != N_levels_single && spin2 != SPIN_DEG &&
            spin1 == spin2 ){
	  
	  if( flag_rec ){

	    if( orbital1_rec != orbital1 || orbital2_rec != orbital2 ){
	      
	      fprintf( stderr, "ERROR: spin-orbital inconsistency [%d,%d]\n", i, j );
	      fflush( stderr );
	      
	      info=1;
	      
	      break;
	      
	    }

	  }
	  else{
	    
	    orbital1_rec = orbital1;
	    spin1_rec    = spin1;
	    
	    orbital2_rec = orbital2;
	    spin2_rec    = spin2;
	    
	    flag_rec=1;
	    
	  }

	}

	if( sign != 0 && spin1 == spin2 ){

	  counter_success++;

	  global_sign += ( symmetry_multiplicity.imatrix[ i ][ 1 ] ) *sign;

	  /* 
	    fprintf( stdout, "counter_success = %d\n", counter_success );
	    fprintf( stdout, "global_sign     = %d\n", global_sign );
	  */

	}


	/*********************
	 * spin_flip i and j *
	 *********************/
	if( symmetry_multiplicity.imatrix[j][0] > 1 ){
	

	  // fprintf( stdout, "spin_flip [ both %d and %d ]\n", i, j );

	  /* initial copy */
	  for( k=0; k<N_levels_single; k++ ){
	  
	    for( l=0; l<SPIN_DEG; l++ ){
		
	      ivec2.ivector[ k*SPIN_DEG +l ] = single_level_occupation.imatrix[ j ][ k*SPIN_DEG +SPIN_DEG -1 -l ];

	    } /* end l loop */
	      
	  } /* end k loop */
	  
	  if( ALLOCATE_MANY_BODY_HYBRIDISATION_TABLE_AUX( *constants_p, ivec1, ivec2, orbital1, spin1, orbital2, spin2, sign ) ) info=1;

	  // fprintf( stdout, "%d, %d, %d, %d, %d\n", orbital1, spin1, orbital2, spin2, sign  );

	  if( orbital1 != N_levels_single && spin1 != SPIN_DEG &&
	      orbital2 != N_levels_single && spin2 != SPIN_DEG && 
              spin1 == spin2 ){
	    
	    if( flag_rec ){

	      if( orbital1_rec != orbital1 || orbital2_rec != orbital2 ){
	      
		fprintf( stderr, "ERROR: spin-orbital inconsistency [%d,%d]\n", i, j );
		fflush( stderr );
		
		info=1;
	      
		break;
	      
	      }

	    }
	    else{

	      orbital1_rec = orbital1;
	      spin1_rec    = spin1;
	      
	      orbital2_rec = orbital2;
	      spin2_rec    = spin2;
	      
	      flag_rec=1;
	      
	    }
	    
	  }

	  if( sign != 0 && spin1 == spin2 ){

	    counter_success++;
	    
	    global_sign += ( symmetry_multiplicity.imatrix[ i ][ 1 ] ) *( symmetry_multiplicity.imatrix[ j ][ 1 ] ) *sign;

	    /* 
	      fprintf( stdout, "counter_success = %d\n", counter_success );
	      fprintf( stdout, "global_sign     = %d\n", global_sign );
	    */

	  }
	  
	} /* end conditional flip j*/

      } /* end conditional flip i*/


      /* write on many_body_hybridisation_table */
      many_body_hybridisation_table_p->imatrix[counter][0] = i;
      many_body_hybridisation_table_p->imatrix[counter][1] = j;
      
      many_body_hybridisation_table_p->imatrix[counter][2] = orbital1_rec;
      many_body_hybridisation_table_p->imatrix[counter][3] = spin1_rec;
      
      many_body_hybridisation_table_p->imatrix[counter][4] = orbital2_rec;
      many_body_hybridisation_table_p->imatrix[counter][5] = spin2_rec;


#ifdef __DEBUG_PLUS__

      fprintf( stdout, "global_sign     = %d\n", global_sign );
      fprintf( stdout, "counter_success = %d\n", counter_success );
      fprintf( stdout, "\n");

#endif /* __DEBUG_PLUS__ */


      symmetry_coefficient_p->vector[ counter ] = CMPLX( (double) ( global_sign ) /sqrt( (double) symmetry_multiplicity.imatrix[ i ][ 0 ] *
											 (double) symmetry_multiplicity.imatrix[ j ][ 0 ] ) );


      /* increase counter */	
      counter++;

      if( info ) break;

    } /* end j loop */
    
    if( info ) break;
    
  } /* end i loop */


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "many_body_hybridisation_table [check]\n" );
  if( IMATRIX_PRINT_PLUS( stdout, *many_body_hybridisation_table_p ) ) info=1;
  fprintf( stdout, "\n" );
    
  fprintf( stdout, "symmetry_coefficients [check]\n" );
  if( VECTOR_PRINT_PLUS( stdout, *symmetry_coefficient_p ) ) info=1;
  fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */ 



  /* deallocation */
  if( IVECTOR_FREE( ivec2 ) ) info=1;

  if( IVECTOR_FREE( ivec1 ) ) info=1;


  return info;

}

//------------------------------------------

int allocate_many_body_hybridisation_table_aux( const constants constants, const ivector ivec1, const ivector ivec2, int* orbital1, int* spin1, int* orbital2, int* spin2, int* sign ){

  /* constants */
  int N_levels_single;
  /* dummies */
  int k;
  int flag1, flag2;
  int flag_fail=0;
  int sign_counter=0;
  int info=0;


  N_levels_single = constants.N_levels_single;


  /* default values */
  flag2 = flag1 = 0;
  *sign         = 0;


#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "original occupation1\n" );
  for( k=0; k<(N_levels_single*SPIN_DEG); k++ ){
	
    fprintf( stdout, "%d ", ivec1.ivector[k] );
	
  }
  fprintf( stdout, "\n" );


  fprintf( stdout, "original occupation2\n" );
  for( k=0; k<(N_levels_single*SPIN_DEG); k++ ){
	
    fprintf( stdout, "%d ", ivec2.ivector[ k ] );
	
  }
  fprintf( stdout, "\n\n" );

#endif /* __DEBUG_PLUS__ */


  for( k=0; k<(N_levels_single*SPIN_DEG); k++ ){
	
    if( ivec1.ivector[k] > ivec2.ivector[ k ] ){
	  
      if( !flag1 ){
	
	*orbital1 = k /SPIN_DEG;
	*spin1    = k %SPIN_DEG;
	
	flag1=1;
	
      }
      else{

	*orbital2 = *orbital1 = N_levels_single;
	*spin2    = *spin1    = SPIN_DEG; 

	flag_fail=1;

	break;
	      
      }

    }
    else if( ivec1.ivector[k] < ivec2.ivector[ k ] ){

      if( !flag2 ){
	    
	*orbital2 = k /SPIN_DEG;
	*spin2    = k %SPIN_DEG;
	
	flag2=1;

      }
      else{

	*orbital2 = *orbital1 = N_levels_single;
	*spin2    = *spin1    = SPIN_DEG; 

	flag_fail=1;

	break;

      }

    }
    else{

      if( !flag1 ){

	sign_counter += ivec1.ivector[ k ];
	  
      }

      if( !flag2 ){
	    
	sign_counter += ivec2.ivector[ k ];
	      
      }
      
    } /* end conditional */
    
  } /* end k loop */

  /* check */
  if( !flag1 && !flag2 ){
	
    fprintf( stderr, "ERROR: two many-body states are equal\n");
    
    info=1;
	  
  }
  else{

    // fermionic sign;
    if( !flag_fail ){

      if( sign_counter%2 ){

	*sign = -1;
	
      }
      else{
	
	*sign = +1;
	
      }

    }
    else{

      *orbital2 = *orbital1 = N_levels_single;
      *spin2    = *spin1    = SPIN_DEG; 
      *sign     = 0;
      
    }

    // fprintf( stdout, "sign = %d\n\n", *sign );
    
  }


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int allocate_compute_initial_rho_electron_workspace( const constants constants ){

  /* constants */
  int N_levels_many;
  /* dummies */
  int info=0;


  N_levels_many = constants.N_levels_many;


#ifdef __DEBUG__

  if( workspace_initial_rho_electron_flag ){

    info=1;

    fprintf( stderr, "ERROR: workspace_initial_rho_electron already allocate\n");
    fflush( stderr );

  }


#endif /* __DEBUG__ */

  if( !info ){

    /* allocation */
    if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, H_tmp_eigenvectors ) ) info=1;

    if( RVECTOR_ALLOCATE( N_levels_many, H_tmp_eigenvalues ) ) info=1;

  }


#ifdef __DEBUG__

  if( !info ){

    workspace_initial_rho_electron_flag=1;

  }
  else{

    workspace_initial_rho_electron_flag=0;

  }

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------

int free_compute_initial_rho_electron_workspace( void ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !workspace_initial_rho_electron_flag ){

    info=1;

    fprintf( stderr, "ERROR: workspace_initial_rho_electron not allocated yet\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */


  if( !info ){

    /* deallocation */
    if( RVECTOR_FREE( H_tmp_eigenvalues ) ) info=1;

    if( MATRIX_FREE( H_tmp_eigenvectors ) ) info=1;

  }


#ifdef __DEBUG__

  if( !info ){

    workspace_initial_rho_electron_flag=0;

  }
  else{

    workspace_initial_rho_electron_flag=1;

  }

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------
//------------------------------------------

#ifdef __DEBUG__

/* reset compute_initial_rho_electron_flag */

int reset_compute_initial_rho_electron_flag( void ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !compute_initial_rho_electron_flag ){

    info=1;

    fprintf( stderr, "ERROR: initial_rho_electron not updated yet\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    compute_initial_rho_electron_flag=0;

  }
  else{

    compute_initial_rho_electron_flag=1;

  }


  return info;

}

//------------------------------------------

/* reset compute_excited_rho_electron_flag */

int reset_compute_excited_rho_electron_flag( void ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( !compute_excited_rho_electron_flag ){

    info=1;

    fprintf( stderr, "ERROR: excited_rho_electron not updated yet\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    compute_excited_rho_electron_flag=0;

  }
  else{

    compute_excited_rho_electron_flag=1;

  }


  return info;

}

//------------------------------------------
//------------------------------------------

/* force compute_initial_rho_electron_flag */

int force_compute_initial_rho_electron_flag( void ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( compute_initial_rho_electron_flag ){

    info=1;

    fprintf( stderr, "ERROR: initial_rho_electron already updated\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    compute_initial_rho_electron_flag=1;

  }
  else{

    compute_initial_rho_electron_flag=0;

  }


  return info;

}

//------------------------------------------

/* force compute_excited_rho_electron_flag */

int force_compute_excited_rho_electron_flag( void ){

  /* dummies */
  int info=0;


#ifdef __DEBUG__

  if( compute_excited_rho_electron_flag ){

    info=1;

    fprintf( stderr, "ERROR: excited_rho_electron already updated\n");
    fflush( stderr );

  }

#endif /* __DEBUG__ */

  if( !info ){

    compute_excited_rho_electron_flag=1;

  }
  else{

    compute_excited_rho_electron_flag=0;

  }


  return info;

}

//------------------------------------------
//------------------------------------------

/* test_compute_initial_rho_electron_flag */

int test_compute_initial_rho_electron_flag( void ){

  /* dummies */
  int info=0;


  if( compute_initial_rho_electron_flag ) info=1;


  return info;

}

//------------------------------------------

/* test_compute_excited_rho_electron_flag */

int test_compute_excited_rho_electron_flag( void ){

  /* dummies */
  int info=0;


  if( compute_excited_rho_electron_flag ) info=1;


  return info;

}

//------------------------------------------

#endif /* __DEBUG__ */
