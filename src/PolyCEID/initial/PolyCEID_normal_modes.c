
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
#include "PolyCEID_normal_modes.h"


#define EPS_LOC   1000 *( DBL_EPSILON )


/*********************
  FUNCTIONS & MACROS
*********************/

//------------------------------------------

/* compute phonons */

int compute_phonons( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p ){

  /* constants  */
  unsigned short int flag_relaxation;
  /* state */
  rvector_p          positions_def_p;
  rvector_p          positions_tmp_p;
  /* dummies */
  //  int i;
  //  double com;
  int info=0;


  flag_relaxation = constants.flag_relaxation;

  positions_def_p = &(config_def_p->atoms.positions);
  positions_tmp_p = &(config_tmp_p->atoms.positions);


#ifdef __DEBUG__

  fprintf( stdout, "DOING: compute_phonons \n" );

#endif /* __DEBUG__ */

  if( RVECTOR_COPY( *positions_tmp_p, *positions_def_p ) ) info=1;

  // RVECTOR_PRINT_PLUS( stdout, *positions_tmp_p );	

  if( RELAXATION( constants, *state_p, *config_tmp_p ) ) info=1;

  // RVECTOR_PRINT_PLUS( stdout, *positions_tmp_p );	

  if( COMPUTE_ORIGINAL_HESSIAN( constants, *state_p, *config_tmp_p ) ) info=1;

#ifdef __CHECK_SYMMETRIES__

  if( CHECK_SYMMETRIES( constants, *state_p, *config_tmp_p ) ) info=1;

#endif /* __CHECK_SYMMETRIES__ */

  if( COMPUTE_NORMAL_MODES( constants, *state_p, *config_def_p ) ) info=1;

  if( flag_relaxation ){

#ifdef __DEBUG__

    fprintf( stdout, "#---------------------------------------------------#\n" );
    fprintf( stdout, "# WARNING, CHANGED INITIAL CONDITION: relaxed positions\n" );
    fprintf( stdout, "#---------------------------------------------------#\n" );
    fprintf( stdout, "\n" );

#endif /* __DEBUG__ */

    if( RVECTOR_COPY( *positions_def_p, *positions_tmp_p ) ) info=1;

  }


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

/* compute original_hessian */

/*
  WARNING: the hessian has been rescaled by taking in to account the atom masses! 
*/

// BUGFIX: this was originally TMP
int compute_original_hessian( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_coor;
  double          mass_ave;
  rvector         masses;
  /* state */
  matrix_p        initial_rho_electron_p;
  matrix_array_p  K_matrix_p;
  matrix_p        dummy_matrix1_p;
  matrix_p        dummy_matrix2_p;
  matrix_p        original_hessian_p;
  /* dummies */
  int index;
  int i, j;
  int info=0;
  complex trace;


  N_coor                   =  constants.N_coor;
  mass_ave                 =  config_p->atoms.mass_ave;
  masses                   =  config_p->atoms.masses;

  initial_rho_electron_p   = &(state_p->initial_rho_electron);
  K_matrix_p               = &(config_p->electrons.K_matrix);
  dummy_matrix1_p          = &state_p->dummy_matrix1;
  dummy_matrix2_p          = &state_p->dummy_matrix2;
  original_hessian_p       = &(state_p->phonons.original_hessian);


#ifdef __DEBUG__

  fprintf( stdout, "DOING: compute_original_hessian \n" );

#endif /* __DEBUG__ */


  if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;

  if( INITIALISE_EHRENFEST_FRAME( constants, *state_p, *config_p ) ) info=1;

  if( K_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;


#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "K_matrix \n" );
  if( MATRIX_ARRAY_PRINT_PLUS( stdout, *K_matrix_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  for( i=0; i<N_coor; i++ ){

    for( j=0; j<N_coor; j++ ){

      index=COORDINATE_INDEX( i, j );

      /*
	fprintf( stdout, "i = %d, j = %d, index = %d\n", i, j, index );
	
	fprintf( stdout, "K_matrix\n" );
	if( MATRIX_PRINT_PLUS( stdout, K_matrix_p->array[ index ]  ) ) info=1;
	fprintf( stdout, "\n" );
      */

      if( MATRIX_COPY( *dummy_matrix1_p, K_matrix_p->array[ index ] ) ) info=1; //WARNING: it might be better to copy and than using BLAS

      if( MATRIX_MATRIX_PRODUCT( *dummy_matrix1_p, *initial_rho_electron_p, *dummy_matrix2_p ) ) info=1;

      /*
	fprintf( stdout, "dummy_matrix2\n" );
	if( MATRIX_PRINT_PLUS( stdout,  *dummy_matrix2_p ) ) info=1;
	fprintf( stdout, "\n" );
      */

      trace = MATRIX_TRACE( *dummy_matrix2_p );


      if( CMPLX_NORM( trace ) > EPS_LOC ){              // WARNING: why regularisation here?

	original_hessian_p->matrix[ index ] = trace;

      }
      else{

	original_hessian_p->matrix[ index ] = CMPLX_ZERO;

      }

    } /* j loop */

  } /* i_loop */


#ifndef __NO_HESSIAN_CORRECTIONS__


#ifdef __DEBUG__

  /* corrections */
  fprintf( stdout, "original_hessian [before corrections]\n" );
  if( MATRIX_PRINT_PLUS( stdout, *original_hessian_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


  /* Hamiltonian */
  if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;

  if( H_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;

  if( F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;


  /* correction to the bare hessian */
  if( HESSIAN_CORRECTIONS( constants, *state_p, *config_p, *original_hessian_p ) ) info=1;


#endif /* __NO_HESSIAN_CORRECTIONS__ */


  /* mass rescaling */
  for( i=0; i<N_coor; i++ ){
    
    for( j=0; j<N_coor; j++ ){
      
      index=COORDINATE_INDEX( i, j );
      
      if( config_p->atoms.mask.ivector[ i ] && config_p->atoms.mask.ivector[ j ] ){

        original_hessian_p->matrix[ index ] = CMPLX_PRODUCT( original_hessian_p->matrix[ index ], 
                                              CMPLX( mass_ave /sqrt( ( masses.rvector[ i ] ) *( masses.rvector[ j ] ) ) ) ); // WARNING: overwriting

      }
					    
    } /* j loop */
    
  } /* i_loop */


#ifdef __DEBUG__

  fprintf( stdout, "original_hessian\n" );
  if( MATRIX_PRINT_PLUS( stdout,  *original_hessian_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------

/* compute normal modes */

/*
  WARNING: the hessian has been rescaled by taking in to account the atom masses! 
*/

int compute_normal_modes( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  int          N_coor;
  int          N_coor_red;
  ivector      relevant_modes;
  double       mass_ave;
  rvector      masses;
  /* state */
  rvector_p    canonical_frequencies_p;
  matrix_p     original_hessian_p;
  matrix_p     canonical_transform_p;
  matrix_p     symplectic_transform_R_p;
  matrix_p     conjugate_symplectic_transform_R_p;
  matrix_p     symplectic_transform_P_p;
  matrix_p     conjugate_symplectic_transform_P_p;
  rvector_p    omega_p;
  rvector_p    ar_p;
  rvector_p    ap_p;
  /* dummies */
  int          i, j; 
  int          index;
  double       average;
  int          i_mode;
  int          info=0;
#ifdef __DEBUG__
  int          i_debug;
  matrix       dummy_matrix1_coor;
  matrix       dummy_matrix2_coor;
  matrix       dummy_matrix3_coor;
#endif /* __DEBUG__ */


  flag_normal_mode_expansion         =  constants.flag_normal_mode_expansion;
  N_coor                             =  constants.N_coor;
  N_coor_red                         =  constants.N_coor_red;
  relevant_modes                     =  constants.relevant_modes;
  mass_ave                           =  config_p->atoms.mass_ave;
  masses                             =  config_p->atoms.masses;

  canonical_frequencies_p            = &(state_p->phonons.canonical_frequencies);
  original_hessian_p                 = &(state_p->phonons.original_hessian);
  canonical_transform_p              = &(state_p->phonons.canonical_transform);
  symplectic_transform_R_p           = &(state_p->phonons.symplectic_transform_R);
  conjugate_symplectic_transform_R_p = &(state_p->phonons.conjugate_symplectic_transform_R);
  symplectic_transform_P_p           = &(state_p->phonons.symplectic_transform_P);
  conjugate_symplectic_transform_P_p = &(state_p->phonons.conjugate_symplectic_transform_P);
  omega_p                            = &(state_p->phonons.omega);
  ar_p                               = &(state_p->phonons.ar);
  ap_p                               = &(state_p->phonons.ap);


#ifdef __DEBUG__

  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix1_coor ) ) info=1;
  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix2_coor ) ) info=1;
  if( MATRIX_ALLOCATE( N_coor, N_coor, dummy_matrix3_coor ) ) info=1;

#endif /* __DEBUG__ */


  if( !info ){

    /* diagonalisation */
    if( DIAGONALISATION( *original_hessian_p, *canonical_transform_p, *canonical_frequencies_p ) ) info=1; 
    
  }
  

#ifdef __DEBUG__

  fprintf( stdout, "canonical transform\n" );
  if( MATRIX_PRINT_PLUS( stdout, *canonical_transform_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "canonical frequencies\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *canonical_frequencies_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


#ifdef __DEBUG__
  if( !info ){

    /* diagonalisation check */
    if( MATRIX_ZERO( dummy_matrix1_coor ) ) info=1;
    
    for( i=0; i<N_coor; i++ ){

      dummy_matrix1_coor.matrix[ COORDINATE_INDEX( i, i ) ].z[ 0 ] = canonical_frequencies_p->rvector[ i ];
      
    }

    if( MATRIX_ADJOINT( *canonical_transform_p, dummy_matrix2_coor ) ) info=1;
    
    if( MATRIX_MATRIX_PRODUCT( dummy_matrix1_coor, dummy_matrix2_coor, dummy_matrix3_coor ) ) info=1;

    if( MATRIX_MATRIX_PRODUCT( *canonical_transform_p, dummy_matrix3_coor, dummy_matrix2_coor ) ) info=1;

    if( MATRIX_DIF( *original_hessian_p, dummy_matrix2_coor, dummy_matrix1_coor ) ) info=1;

    /*
      fprintf( stdout, "matrix difference\n" );
      if( MATRIX_PRINT_PLUS( stdout, dummy_matrix1_coor ) ) info=1;
      fprintf( stdout, "\n" );
    */

    /*
      fprintf( stdout, "norm = %le\n", MATRIX_NORM( dummy_matrix1_coor ) );
      fflush( stdout );
    */

    if( MATRIX_NORM( dummy_matrix1_coor ) > EPS ){

      fprintf( stderr, "ERROR: inconsistency in original_hessian diagonalisation.\n" ); 
      fflush( stderr );

      info=1;

    }

  }

#endif /* _DEBUG__ */


  if( !info ){

    /* normal modes non-negativity check */
    /* WARNING: the centre-of-mass mode is always zero */
    for( i=0; i<N_coor; i++ ){

      if( canonical_frequencies_p->rvector[ i ] < -EPS_LOC ){

	fprintf( stderr, "# ERROR: the %dth normal mode is negative.\n", i );
	fflush( stderr );

        info=1;

      }

      if( fabs( canonical_frequencies_p->rvector[ i ] ) < EPS_LOC ){

        if( flag_normal_mode_expansion ){

 	  config_p->atoms.masses_aux.rvector[ i ] = 0.0e0;      

	}

      }

    } /* end i loop */

  }

  /* metric constants */
  if( !info ){

    /* avaerage of the effective hessian diagonal entries */
    average = 0.0e0;

    if( !flag_normal_mode_expansion ){

      for( i=0; i<N_coor; i++ ){

        average += REAL( original_hessian_p->matrix[ COORDINATE_INDEX( i, i ) ] );

      } /* end i loop */

      average /= N_coor;

      if( average < EPS_LOC ){
      
        fprintf( stderr, "ERROR: the average of the effective hessian diagonal entries must be positive\n" );
        fflush( stderr );
      
        info =1;
      
      }


#ifdef __DEBUG_PLUS__

      fprintf( stdout, "Effective hessian average = %le\n", average );
      fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */

    }

    i_mode=0; // important!

    for( i=0; i<N_coor; i++ ){


      // fprintf( stdout, "---> i_mode = %d, i = %d, N_coor = %d \n", i_mode, i, N_coor );

      
      /* omegas */
      if( flag_normal_mode_expansion ){

        if( canonical_frequencies_p->rvector[ i ] > EPS ){

	  omega_p->rvector[ i ] = sqrt( ( canonical_frequencies_p->rvector[ i ] ) /mass_ave ); // WARNING: notice the use of mass_ave
	
        }
        else{

	  omega_p->rvector[ i ] = 0.0e0;

        }
      
      }
      else{

        //      omega_p->rvector[ i ] = sqrt( REAL( original_hessian_p->matrix[ COORDINATE_INDEX( i, i ) ] ) /mass_ave );
        //      This line has been modified because it does not satify the bosonic symmetry

        omega_p->rvector[ i ] = sqrt( average /mass_ave ); // WARNING: notice the use of mass_ave

      }

      /* relevant omega check */
      if( !info ){
	
	if( i == relevant_modes.ivector[ i_mode ] ){
	  
	  if( omega_p->rvector[ i ] < EPS ){
	    
	    fprintf( stderr, "# WARNING: omega[ %d ] is zero [%12.5f]\n", i, omega_p->rvector[ i ] );
	    fflush( stderr );
	    
	  }
      
	  if( i_mode < ( N_coor_red -1 ) ) i_mode++;
	  
	}
      
      }

      if( omega_p->rvector[ i ] > EPS ){

	ar_p->rvector[ i ]    = sqrt( HBAR /( omega_p->rvector[ i ] *mass_ave ) ); // WARNING: notice the use of mass_ave
	ap_p->rvector[ i ]    = HBAR /( ar_p->rvector[ i ] );
	
	/* 
	   warning: WRONG NORMALISATION!
	   ar_p->rvector[ i ]    = sqrt( 0.5e0 *HBAR /( omega_p->rvector[ i ] *mass_ave ) );
	   ap_p->rvector[ i ]    = 0.5e0 *HBAR /( ar_p->rvector[ i ] );
	*/
	
      }
      else{

	ar_p->rvector[ i ]    = 0.0e0;
	ap_p->rvector[ i ]    = 0.0e0;

      }

    } /* end i loop */


#ifdef __DEBUG__

    // Print the effective oscillator energies
    fprintf( stdout, "# -----------------------------\n" );
    fprintf( stdout, "# Effective oscillator energies\n" );
    for( i_debug=0; i_debug<N_coor; i_debug++ ){

      fprintf( stdout, "# omega[%d] = %le\n", i_debug, HBAR *omega_p->rvector[ i_debug ] );

    }        
    fprintf( stdout, "# -----------------------------\n" );
    fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


  }

  if( !info ){

    /* tentative phononic DOS */
    COMPUTE_PHONONIC_DOS( constants, *state_p );

  }

  /* symplectic matrix generation */
  if( flag_normal_mode_expansion ){

    if( !info ){

      for( i=0; i<N_coor; i++ ){
      
        for( j=0; j<N_coor; j++ ){
	
	  index=COORDINATE_INDEX( i, j );
	
          if( config_p->atoms.mask.ivector[ i ] ){

            symplectic_transform_R_p->matrix[ index ] = CMPLX_PRODUCT( CMPLX( sqrt( mass_ave /( masses.rvector[ i ] ) ) ), canonical_transform_p->matrix[ index ] );

	  }
	  else{

            symplectic_transform_R_p->matrix[ index ] = canonical_transform_p->matrix[ index ]; 
		  
          }

        } /* j loop */

      } /* i_loop */


#ifdef __DEBUG__

      fprintf( stdout, "symplectic transform R\n" );
      if( MATRIX_PRINT_PLUS( stdout, *symplectic_transform_R_p ) ) info=1;
      fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


      if( MATRIX_ADJOINT( *symplectic_transform_R_p, *conjugate_symplectic_transform_R_p ) ) info=1;
    
      if( MATRIX_INVERSE( *conjugate_symplectic_transform_R_p, *symplectic_transform_P_p ) ) info=1; //WARNING: this is very important

      if( MATRIX_ADJOINT( *symplectic_transform_P_p, *conjugate_symplectic_transform_P_p ) ) info=1;

    }

  }


#ifdef __DEBUG__

  if( MATRIX_FREE( dummy_matrix3_coor ) ) info=1;
  if( MATRIX_FREE( dummy_matrix2_coor ) ) info=1;
  if( MATRIX_FREE( dummy_matrix1_coor ) ) info=1;

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------
//------------------------------------------

/*  check_symmetries */

//BUGFIX: this was originally TMP
int check_symmetries( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int       N_coor;
  int       N_atoms;
  /* state */
  matrix_p  original_hessian_p;
  rvector_p positions_p;
  /* dummies */
  vector    dummy_vector_coor1;
  vector    dummy_vector_coor2;
  rvector   centroid;
  double    norm;
  //  complex   scalar;
  int       sdim;
  int       i, j;
  int       i_atom;
  int       i_comp;
  int       info=0;


  N_coor              =  constants.N_coor;
  N_atoms             =  constants.N_atoms;
  sdim                =  constants.spacial_dimension; 

  original_hessian_p  = &(state_p->phonons.original_hessian);
  positions_p         = &(config_p->atoms.positions);


  if( VECTOR_ALLOCATE( N_coor, dummy_vector_coor1 ) ) info=1;

  if( VECTOR_ALLOCATE( N_coor, dummy_vector_coor2 ) ) info=1;

  if( RVECTOR_ALLOCATE( sdim, centroid ) ) info=1;

#ifdef __DEBUG__

  fprintf( stdout, "checking hessian symmetries\n" );
  fprintf( stdout, "\n" );

  /*
    fprintf( stdout, "original_hessian\n" );
    if( MATRIX_PRINT_PLUS( stdout,  *original_hessian_p ) ) info=1;
    fprintf( stdout, "\n" );
  */

  //
  fprintf( stdout, "checking translation [1]\n" );
  fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


  if( VECTOR_ZERO( dummy_vector_coor1 ) ) info=1;

  for( i=0; i<N_coor; i++ ){

    if( i %sdim == 0 ){

      dummy_vector_coor1.vector[ i ] = CMPLX_UNIT;

    }

  }

#ifdef __DEBUG__

  fprintf( stdout, "original vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor1 ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


  if( MATRIX_VECTOR_PRODUCT( *original_hessian_p, dummy_vector_coor1, dummy_vector_coor2 ) ) info=1;


#ifdef __DEBUG__

  fprintf( stdout, "transformed vector\n" );
  if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor2 ) ) info=1;
  fprintf( stdout, "\n" );
  
#endif /* __DEBUG__ */

  if( !info ){

    norm = VECTOR_NORM( dummy_vector_coor2 );
    
    if( norm > EPS ){
      
      fprintf( stderr, "ERROR: the original_hessian is not invariant\n" );
      
      info=1;
      
    }

  }


  if( sdim > 1 ){

    if( !info ){

#ifdef __DEBUG__

      fprintf( stdout, "checking translation [2]\n" );
      fprintf( stdout, "\n" );
      
#endif /* __DEBUG__ */

      
      if( VECTOR_ZERO( dummy_vector_coor1 ) ) info=1;
      
      for( i=0; i<N_coor; i++ ){
	
	if( i %sdim == 1 ){

	  dummy_vector_coor1.vector[ i ] = CMPLX_UNIT;
	
	}
	
      }
      
#ifdef __DEBUG__
      
      fprintf( stdout, "original vector\n" );
      if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor1 ) ) info=1;
      fprintf( stdout, "\n" );
      
#endif /* __DEBUG__ */


      if( MATRIX_VECTOR_PRODUCT( *original_hessian_p, dummy_vector_coor1, dummy_vector_coor2 ) ) info=1;


#ifdef __DEBUG__
    
      fprintf( stdout, "transformed vector\n" );
      if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor2 ) ) info=1;
      fprintf( stdout, "\n" );
    
#endif /* __DEBUG__ */

    }

    if( !info ){

      norm = VECTOR_NORM( dummy_vector_coor2 );
      
      if( norm > EPS ){
	
	fprintf( stderr, "ERROR: the original_hessian is not invariant\n" );
	
	info=1;
      
      }
      
    }
    

    if( sdim > 2 ){

      if( !info ){

#ifdef __DEBUG__

	fprintf( stdout, "checking translation [3]\n" );
	fprintf( stdout, "\n" );
	
#endif /* __DEBUG__ */   

	
	if( VECTOR_ZERO( dummy_vector_coor1 ) ) info=1;
	
	for( i=0; i<N_coor; i++ ){
	  
	  if( i %sdim == 2 ){
	    
	    dummy_vector_coor1.vector[ i ] = CMPLX_UNIT;
	    
	  }
	  
	}
	
#ifdef __DEBUG__
	
	fprintf( stdout, "original vector\n" );
	if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor1 ) ) info=1;
	fprintf( stdout, "\n" );
	
#endif /* __DEBUG__ */
	
	
	if( MATRIX_VECTOR_PRODUCT( *original_hessian_p, dummy_vector_coor1, dummy_vector_coor2 ) ) info=1;


#ifdef __DEBUG__
    
	fprintf( stdout, "transformed vector\n" );
	if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor2 ) ) info=1;
	fprintf( stdout, "\n" );
      
    
#endif /* __DEBUG__ */

      }
     
      if( !info ){

	norm = VECTOR_NORM( dummy_vector_coor2 );
	
	if( norm > EPS ){
      
	  fprintf( stderr, "ERROR: the original_hessian is not invariant\n" );
	  
	  info=1;
      
	}
	
      }
      
    }    
    
  }


  /* rotations */
  if( sdim > 1 ){

    if( !info ){

#ifdef __DEBUG__
      
      fprintf( stdout, "checking rotation [1]\n" );
      fprintf( stdout, "\n" );
      
#endif /* __DEBUG__ */
  
      // centroid
      for( i=0; i<N_atoms; i++ ){
	
	for( j=0; j<sdim; j++ ){
	  
	  centroid.rvector[ j ] += positions_p->rvector[ j +i *sdim ];
	  
	}

      }
      
      for( j=0; j<sdim; j++ ){
	
	centroid.rvector[ j ] /= N_atoms;
      
      }

#ifdef __DEBUG__

      fprintf( stdout, "centroid\n" );
      if( RVECTOR_PRINT_PLUS( stdout, centroid ) ) info=1;
      fprintf( stdout, "\n" );
    
#endif /* __DEBUG__ */


      if( VECTOR_ZERO( dummy_vector_coor1 ) ) info=1;

      for( i=0; i<N_coor; i++ ){
	
	i_atom = i /sdim;
	i_comp = i %sdim;
      
	if( i_comp == 0 ){

	  dummy_vector_coor1.vector[ i ] = CMPLX(  positions_p->rvector[ 1 +sdim *i_atom ] -centroid.rvector[ 1 ] );

	}
	else if( i_comp == 1 ){
	
	  dummy_vector_coor1.vector[ i ] = CMPLX( -positions_p->rvector[ 0 +sdim *i_atom ] +centroid.rvector[ 0 ] );

	}
	
      }

#ifdef __DEBUG__

      fprintf( stdout, "original vector\n" );
      if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor1 ) ) info=1;
      fprintf( stdout, "\n" );
      
#endif /* __DEBUG__ */

    
      if( MATRIX_VECTOR_PRODUCT( *original_hessian_p, dummy_vector_coor1, dummy_vector_coor2 ) ) info=1;

#ifdef __DEBUG__

      fprintf( stdout, "transformed vector\n" );
      if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor2 ) ) info=1;
      fprintf( stdout, "\n" );

#endif /* __DEBUG__ */
      
    }

    if( !info ){

      norm = VECTOR_NORM( dummy_vector_coor2 );
      
      if( norm > EPS ){
      
	fprintf( stderr, "ERROR: the original_hessian is not invariant\n" );
	
	info=1;
	
      }
      
    }


    /*
      scalar = SCALAR_PRODUCT( dummy_vector_coor2, dummy_vector_coor1 );
    
      
      fprintf( stdout, "scalar_product\n" );
      if( CMPLX_PRINT_PLUS( stdout, scalar ) ) info=1;
      fprintf( stdout, "\n\n" );
    
      
      for( i=0; i<N_coor; i++ ){

      if( CMPLX_NORM( dummy_vector_coor1.vector[i] ) > EPS ){
      
      scalar = CMPLX_DIVISION( dummy_vector_coor2.vector[i], dummy_vector_coor1.vector[i] );

      if( CMPLX_PRINT_PLUS( stdout, scalar ) ) info=1;
      fprintf( stdout, "\n" );

      }
      else{
      
      if( CMPLX_PRINT_PLUS( stdout, CMPLX_ZERO ) ) info=1;
      fprintf( stdout, "\n" );

      }

      }
      fprintf( stdout, "\n" );

    */

    if( sdim > 2 ){

      if( !info ){

#ifdef __DEBUG__
	
	fprintf( stdout, "checking rotation [2]\n" );
	fprintf( stdout, "\n" );
	
#endif /* __DEBUG__ */

	if( VECTOR_ZERO( dummy_vector_coor1 ) ) info=1;
	
	for( i=0; i<N_coor; i++ ){
	  
	  i_atom = i /sdim;
	  i_comp = i %sdim;
      
	  if( i_comp == 1 ){
	
	    dummy_vector_coor1.vector[ i ] = CMPLX(  positions_p->rvector[ 2 +sdim *i_atom ] -centroid.rvector[ 2 ]);
	
	  }
	  else if( i_comp == 2 ){
	    
	    dummy_vector_coor1.vector[ i ] = CMPLX( -positions_p->rvector[ 1 +sdim *i_atom ] +centroid.rvector[ 1 ] );
	    
	  }
	  
	}
    
#ifdef __DEBUG__
    
	fprintf( stdout, "original vector\n" );
	if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor1 ) ) info=1;
	fprintf( stdout, "\n" );
	
#endif /* __DEBUG__ */

      
	if( MATRIX_VECTOR_PRODUCT( *original_hessian_p, dummy_vector_coor1, dummy_vector_coor2 ) ) info=1;
    
#ifdef __DEBUG__
	
	fprintf( stdout, "transformed vector\n" );
	if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor2 ) ) info=1;
      fprintf( stdout, "\n" );
      
#endif /* __DEBUG__ */

      }

      if( !info ){

	norm = VECTOR_NORM( dummy_vector_coor2 );
	
	if( norm > EPS ){
	  
	  fprintf( stderr, "ERROR: the original_hessian is not invariant\n" );
	  
	  info=1;
      
	}
	
      }
      
      
      if( !info ){

#ifdef __DEBUG__
	
	fprintf( stdout, "checking rotation [3]\n" );
	fprintf( stdout, "\n" );
	
#endif /* __DEBUG__ */
	
	if( VECTOR_ZERO( dummy_vector_coor1 ) ) info=1;
	
	for( i=0; i<N_coor; i++ ){
	  
	  i_atom = i /sdim;
	  i_comp = i %sdim;
	  
	  if( i_comp == 2 ){
	    
	    dummy_vector_coor1.vector[ i ] = CMPLX(  positions_p->rvector[ 0 +sdim *i_atom ] -centroid.rvector[ 0 ] );
	    
	  }
	  else if( i_comp == 0 ){
	  
	    dummy_vector_coor1.vector[ i ] = CMPLX( -positions_p->rvector[ 2 +sdim *i_atom ] +centroid.rvector[ 2 ] );
	
	  }
	
	}
    
#ifdef __DEBUG__
	
	fprintf( stdout, "original vector\n" );
	if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor1 ) ) info=1;
	fprintf( stdout, "\n" );
	
#endif /* __DEBUG__ */
      
      
	if( MATRIX_VECTOR_PRODUCT( *original_hessian_p, dummy_vector_coor1, dummy_vector_coor2 ) ) info=1;
    
#ifdef __DEBUG__
	
	fprintf( stdout, "transformed vector\n" );
	if( VECTOR_PRINT_PLUS( stdout, dummy_vector_coor2 ) ) info=1;
	fprintf( stdout, "\n" );
      
#endif /* __DEBUG__ */

      }

      if( !info ){
	
	norm = VECTOR_NORM( dummy_vector_coor2 );
    
	if( norm > EPS ){
      
	  fprintf( stderr, "ERROR: the original_hessian is not invariant\n" );
      
	  info=1;
      
	}

      }

    }

  }

  //
  if( RVECTOR_FREE( centroid ) )          info=1;
  
  if( VECTOR_FREE( dummy_vector_coor2 ) ) info=1;
  
  if( VECTOR_FREE( dummy_vector_coor1 ) ) info=1;
  

  return info;

}

//------------------------------------------

/* hessian corrections */

/*
  WARNING: the hessian has been rescaled by taking in to account the atom masses! 
*/

// BUGFIX: this was originally TMP
int hessian_corrections( const constants constants, state_p state_p, config_p config_p, matrix_p hessian_p ){

  /* constants */
  int              N_levels_many;
  int              N_coor;
  rvector          initial_many_body_state;
  /* state */
  matrix_p         initial_rho_electron_p;
  matrix_p         H_matrix_p;
  matrix_array_p   F_matrix_p;
  rvector_p        dummy_rvector_p;
  vector_p         dummy_vector1_p;
  vector_p         dummy_vector2_p;
  matrix_p         dummy_matrix1_p;
  matrix_p         dummy_matrix2_p;
  /* dummies */
  int              i_coor;
  int              i_level;
  int              i, j;
  complex          dummy;
  int              index;
  double           eigenenergy_loc;
  double           delta_energy;
  vector           eigenstate_loc;
  vector           F_loc;
  matrix           hessian_corrections;
#ifdef __DEBUG_PLUS__
  matrix           hessian_corrections_eigenvectors;
  rvector          hessian_corrections_eigenvalues;
#endif /* __DEBUG_PLUS__ */
  int info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: hessian corrections\n" );

#endif /* __DEBUG_PLUS__ */


  N_levels_many           =  constants.N_levels_many;
  N_coor                  =  constants.N_coor;
  initial_many_body_state =  constants.initial_many_body_state;

  initial_rho_electron_p  = &(state_p->initial_rho_electron);
  H_matrix_p              = &(config_p->electrons.H_matrix);
  F_matrix_p              = &(config_p->electrons.F_matrix);
  dummy_rvector_p         = &(state_p->dummy_rvector);
  dummy_vector1_p         = &(state_p->dummy_vector1);
  dummy_vector2_p         = &(state_p->dummy_vector2);
  dummy_matrix1_p         = &(state_p->dummy_matrix1);
  dummy_matrix2_p         = &(state_p->dummy_matrix2);
 

  /* allocations */
  if( VECTOR_ALLOCATE( N_levels_many, eigenstate_loc ) )       info=1;

  if( VECTOR_ALLOCATE( N_coor, F_loc ) )                       info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, hessian_corrections ) ) info=1;

#ifdef __DEBUG_PLUS__

  if( MATRIX_ALLOCATE( N_coor, N_coor, hessian_corrections_eigenvectors ) ) info=1;

  if( RVECTOR_ALLOCATE( N_coor, hessian_corrections_eigenvalues ) ) info=1;

#endif /* __DEBUG_PLUS__ */

  /* Hamiltonian */
  //  if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;

  //  if( H_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;

  //  if( F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "H_matrix\n" );
  if( MATRIX_PRINT_PLUS( stdout, *H_matrix_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* Ehrenfest Hamiltonian diagonalisation */
  if( DIAGONALISATION( *H_matrix_p, *dummy_matrix1_p, *dummy_rvector_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "eigenvectors\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "eigenvalues\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *dummy_rvector_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* eigenstate_loc */
  for( i=0; i<N_levels_many; i++ ){
     
    eigenstate_loc.vector[ i ] = CMPLX( initial_many_body_state.rvector[ i ] );
  
  }


#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "eigenstate_loc\n" );
  if( VECTOR_PRINT_PLUS( stdout, eigenstate_loc ) ) info=1;
  fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */

  
  /* eigenenergy_loc */
  if( MATRIX_MATRIX_PRODUCT( *H_matrix_p, *initial_rho_electron_p, *dummy_matrix2_p ) ) info=1; 

  eigenenergy_loc = REAL( MATRIX_TRACE( *dummy_matrix2_p ) );



#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "eigenenergy_loc\n" );
  fprintf( stdout, "%le\n", eigenenergy_loc );
  fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */


  /* set to zero */
  if( MATRIX_ZERO( hessian_corrections ) ) info=1;

  for( i_level=0; i_level<N_levels_many; i_level++ ){
   
    /* create F_loc */
    for( i_coor=0; i_coor<N_coor; i_coor++ ){
      
      /* virtual eigenstate copy */
      for( i=0; i<N_levels_many; i++ ){
	  
	dummy_vector1_p->vector[ i ] = dummy_matrix1_p->matrix[ ELECTRON_MANY_INDEX( i, i_level ) ];
	  
      }
	
      if( MATRIX_VECTOR_PRODUCT( F_matrix_p->array[ i_coor ], *dummy_vector1_p, *dummy_vector2_p ) ) info=1;
	
      F_loc.vector[ i_coor ] = SCALAR_PRODUCT( eigenstate_loc, *dummy_vector2_p );
	
    } /* end F_loc */
      

#ifdef __DEBUG_PLUS__

    fprintf( stdout, "F_loc\n" );
    if( VECTOR_PRINT_PLUS( stdout, F_loc ) ) info=1;
    fprintf( stdout, "\n" );
    
#endif /* __DEBUG_PLUS__ */


    /* degenerancy check */
    delta_energy = eigenenergy_loc - dummy_rvector_p->rvector[ i_level ];
    
    if( fabs( delta_energy ) > EPS ){
      
      /* update hessian_corrections */
      for( i=0; i<N_coor; i++ ){
	
	for( j=0; j<N_coor; j++ ){
	  
	  index = COORDINATE_INDEX( i, j );
	    
	  dummy = CMPLX_PRODUCT( F_loc.vector[ i ], CONJ( F_loc.vector[ j ] ) );  

          dummy = CMPLX_PRODUCT( CMPLX(2.0e0), dummy ); //WARNING: this is due to the definition of the Hessian matrix
	    
	  dummy = CMPLX_DIVISION( dummy, CMPLX( delta_energy ) ); // WARNING: overwriting
	    
	  hessian_corrections.matrix[ index ] = CMPLX_SUM( hessian_corrections.matrix[ index ], dummy ); // WARNING: overwriting

	} /* j loop */
	  
      } /* i loop */
	
    } /* delta_energy conditional */

  } /* i_level loop */
    

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "hessian corrections\n" );
  if( MATRIX_PRINT_PLUS( stdout, hessian_corrections) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG_PLUS__

  /* hessian corrections diagonalisation */
  if( DIAGONALISATION( hessian_corrections, hessian_corrections_eigenvectors, hessian_corrections_eigenvalues ) ) info=1;

  fprintf( stdout, "eigenvectors\n" );
  if( MATRIX_PRINT_PLUS( stdout, hessian_corrections_eigenvectors ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "eigenvalues\n" );
  if( RVECTOR_PRINT_PLUS( stdout, hessian_corrections_eigenvalues ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG_PLUS__

  /* corrections */
  fprintf( stdout, "hessian [before corrections]\n" );
  if( MATRIX_PRINT_PLUS( stdout, *hessian_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  if( MATRIX_SUM( *hessian_p, hessian_corrections, *hessian_p ) ) info=1; // WARNING: overwriting


#ifdef __DEBUG_PLUS__

  /* corrections */
  fprintf( stdout, "hessian [after corrections]\n" );
  if( MATRIX_PRINT_PLUS( stdout, *hessian_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* deallocations */

#ifdef __DEBUG_PLUS__

  if( RVECTOR_FREE( hessian_corrections_eigenvalues ) ) info=1;

  if( MATRIX_FREE( hessian_corrections_eigenvectors ) ) info=1;

#endif /* __DEBUG_PLUS__ */

  if( MATRIX_FREE( hessian_corrections ) ) info=1;

  if( VECTOR_FREE( F_loc ) ) info=1;

  if( VECTOR_FREE( eigenstate_loc ) ) info=1;


  return info;

}

//------------------------------------------

/* hessian corrections_aux */

/*
  WARNING: the hessian has been rescaled by taking in to account the atom masses! 
*/

//BUGFIX: this was orginally DEF
int hessian_corrections_aux( const constants constants, state_p state_p, config_p config_p, vector electronic_state, matrix_p hessian_p ){

  /* constants */
  int              N_levels_many;
  int              N_coor;
  /* state */
  matrix_p         H_matrix_p;
  matrix_array_p   F_matrix_p;
  rvector_p        dummy_rvector_p;
  vector_p         dummy_vector1_p;
  vector_p         dummy_vector2_p;
  matrix_p         dummy_matrix1_p;
  /* dummies */
  int              i_coor;
  int              i_level;
  int              i, j;
  complex          dummy;
  int              index;
  double           eigenenergy_loc;
  double           delta_energy;
  vector           F_loc;
  matrix           hessian_corrections;
#ifdef __DEBUG_PLUS__
  matrix           hessian_corrections_eigenvectors;
  rvector          hessian_corrections_eigenvalues;
#endif /* __DEBUG_PLUS__ */
  int info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: hessian corrections\n" );

#endif /* __DEBUG_PLUS__ */


  N_levels_many           =  constants.N_levels_many;
  N_coor                  =  constants.N_coor;

  H_matrix_p              = &(config_p->electrons.H_matrix);
  F_matrix_p              = &(config_p->electrons.F_matrix);
  dummy_rvector_p         = &(state_p->dummy_rvector);
  dummy_vector1_p         = &(state_p->dummy_vector1);
  dummy_vector2_p         = &(state_p->dummy_vector2);
  dummy_matrix1_p         = &(state_p->dummy_matrix1);
 

  /* allocations */
  if( VECTOR_ALLOCATE( N_coor, F_loc ) )                       info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, hessian_corrections ) ) info=1;

#ifdef __DEBUG_PLUS__

  if( MATRIX_ALLOCATE( N_coor, N_coor, hessian_corrections_eigenvectors ) ) info=1;

  if( RVECTOR_ALLOCATE( N_coor, hessian_corrections_eigenvalues ) ) info=1;

#endif /* __DEBUG_PLUS__ */

  /* Hamiltonian */
  //  if( DISTANCES_UPDATE( constants, *state_p, *config_p ) ) info=1;

  //  if( H_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;

  //  if( F_MATRIX_MANY_UPDATE( constants, *state_p, *config_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "H_matrix\n" );
  if( MATRIX_PRINT_PLUS( stdout, *H_matrix_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* Ehrenfest Hamiltonian diagonalisation */
  if( DIAGONALISATION( *H_matrix_p, *dummy_matrix1_p, *dummy_rvector_p ) ) info=1;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "eigenvectors\n" );
  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "eigenvalues\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *dummy_rvector_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* eigenenergy_loc */
  if( MATRIX_VECTOR_PRODUCT( *H_matrix_p, electronic_state, *dummy_vector1_p ) ) info=1; 

  eigenenergy_loc = REAL( SCALAR_PRODUCT( electronic_state, *dummy_vector1_p ) );



#ifdef __DEBUG_PLUS__
  
  fprintf( stdout, "eigenenergy_loc\n" );
  fprintf( stdout, "%le\n", eigenenergy_loc );
  fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */


  /* set to zero */
  if( MATRIX_ZERO( hessian_corrections ) ) info=1;

  for( i_level=0; i_level<N_levels_many; i_level++ ){
   
    /* create F_loc */
    for( i_coor=0; i_coor<N_coor; i_coor++ ){
      
      /* virtual eigenstate copy */
      for( i=0; i<N_levels_many; i++ ){
	  
	dummy_vector1_p->vector[ i ] = dummy_matrix1_p->matrix[ ELECTRON_MANY_INDEX( i, i_level ) ];
	  
      }
	
      if( MATRIX_VECTOR_PRODUCT( F_matrix_p->array[ i_coor ], *dummy_vector1_p, *dummy_vector2_p ) ) info=1;
	
      F_loc.vector[ i_coor ] = SCALAR_PRODUCT( electronic_state, *dummy_vector2_p );
	
    } /* end F_loc */
      

#ifdef __DEBUG_PLUS__

    fprintf( stdout, "F_loc\n" );
    if( VECTOR_PRINT_PLUS( stdout, F_loc ) ) info=1;
    fprintf( stdout, "\n" );
    
#endif /* __DEBUG_PLUS__ */


    /* degenerancy check */
    delta_energy = eigenenergy_loc - dummy_rvector_p->rvector[ i_level ];
    
    if( fabs( delta_energy ) > EPS ){
      
      /* update hessian_corrections */
      for( i=0; i<N_coor; i++ ){
	
	for( j=0; j<N_coor; j++ ){
	  
	  index = COORDINATE_INDEX( i, j );
	    
	  dummy = CMPLX_PRODUCT( F_loc.vector[ i ], CONJ( F_loc.vector[ j ] ) );  

          dummy = CMPLX_PRODUCT( CMPLX(2.0e0), dummy ); //WARNING: this is due to the definition of the Hessian matrix
	    
	  dummy = CMPLX_DIVISION( dummy, CMPLX( delta_energy ) ); // WARNING: overwriting
	    
	  hessian_corrections.matrix[ index ] = CMPLX_SUM( hessian_corrections.matrix[ index ], dummy ); // WARNING: overwriting

	} /* j loop */
	  
      } /* i loop */
	
    } /* delta_energy conditional */

  } /* i_level loop */
    

#ifdef __DEBUG_PLUS__

  fprintf( stdout, "hessian corrections\n" );
  if( MATRIX_PRINT_PLUS( stdout, hessian_corrections) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG_PLUS__

  /* hessian corrections diagonalisation */
  if( DIAGONALISATION( hessian_corrections, hessian_corrections_eigenvectors, hessian_corrections_eigenvalues ) ) info=1;

  fprintf( stdout, "eigenvectors\n" );
  if( MATRIX_PRINT_PLUS( stdout, hessian_corrections_eigenvectors ) ) info=1;
  fprintf( stdout, "\n" );

  fprintf( stdout, "eigenvalues\n" );
  if( RVECTOR_PRINT_PLUS( stdout, hessian_corrections_eigenvalues ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG_PLUS__

  /* corrections */
  fprintf( stdout, "hessian [before corrections]\n" );
  if( MATRIX_PRINT_PLUS( stdout, *hessian_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  if( MATRIX_SUM( *hessian_p, hessian_corrections, *hessian_p ) ) info=1; // WARNING: overwriting


#ifdef __DEBUG_PLUS__

  /* corrections */
  fprintf( stdout, "hessian [after corrections]\n" );
  if( MATRIX_PRINT_PLUS( stdout, *hessian_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG_PLUS__ */


  /* deallocations */

#ifdef __DEBUG_PLUS__

  if( RVECTOR_FREE( hessian_corrections_eigenvalues ) ) info=1;

  if( MATRIX_FREE( hessian_corrections_eigenvectors ) ) info=1;

#endif /* __DEBUG_PLUS__ */

  if( MATRIX_FREE( hessian_corrections ) ) info=1;

  if( VECTOR_FREE( F_loc ) ) info=1;


  return info;

}

//------------------------------------------

/* compute_phononic_dos */

int compute_phononic_dos( const constants constants, state_p state_p ){

  /* constants */
  int     N_coor;
  /* state */
  rvector_p   omega_p;
  rvector_p   phononic_dos_p;
  /* dummies */
  int     i, i_tmp;
  double  E_min, E_max, E_tmp;
  int     info=0;


  N_coor          =  constants.N_coor;
  omega_p         = &(state_p->phonons.omega);
  phononic_dos_p  = &(state_p->phonons.phononic_dos);


  for( i=0; i<N_coor; i++ ){

    // search for the E_min
    i_tmp = i;
    E_min = omega_p->rvector[ i_tmp ];
    while( i_tmp > 0 ){
      
      i_tmp--;
      E_tmp = omega_p->rvector[ i_tmp ];
      if( E_tmp + EPS < E_min ){
	E_min = E_tmp;
	break;
      }

    } /* end while */

    // search for the E_max
    i_tmp = i;
    E_max = omega_p->rvector[ i_tmp ];
    while( i_tmp < N_coor -1 ){
      
      i_tmp++;
      E_tmp = omega_p->rvector[ i_tmp ];
      if( E_tmp - EPS > E_max ){
	E_max = E_tmp;
	break;
      }

    } /* end while */       

    E_tmp = ONEO2 *( E_max -E_min );

    // at the ends
    if( 0 == i || N_coor-1 == i ){
      E_tmp *= TWO;
    }

    phononic_dos_p->rvector[ i ] = ONE /E_tmp;

    // fprintf( stdout, "%d   %le   %le\n", i,  omega_p->rvector[ i ], phononic_dos_p->rvector[ i ] );

  } /* end i loop */


  return info;

}  

//------------------------------------------
