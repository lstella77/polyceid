
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
#include "PolyCEID_hamiltonian_single_CHAIN_MOD2.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* H_matrix_single_CHAIN_MOD2 update */

int H_matrix_single_CHAIN_MOD2_update( const constants constants, state_p state_p, config_p config_p, matrix_p matrix_p ){

  /* constants */
  int        N_electrons;
  int        N_atoms;
  double     a_spacing;
  //  double     delta_energy;
  double     k_harmonic;
  double     t_hopping;
  double     alpha_coupling;
  double     dimerisation;
  /* state */
  rvector_p  positions_p;
  /* dummies */
  int        index;
  int        i;
  int        info=0;
  double     distance;
  double     dummy;
  double     dummy2;
  double     extreem1;
  double     extreem2;
  double     sigma;


  N_electrons       =  constants.N_electrons;
  N_atoms           =  constants.N_atoms;
  a_spacing         =  constants.hamiltonian.par_extra[0];
  // delta_energy      =  constants.hamiltonian.par_extra[1];
  k_harmonic        =  constants.hamiltonian.par_extra[2]/(double) N_electrons; // WARNING: rescaling
  t_hopping         =  constants.hamiltonian.par_extra[3];
  alpha_coupling    =  constants.hamiltonian.par_extra[4];
  dimerisation      =  constants.hamiltonian.par_extra[5];
  extreem1          =  constants.hamiltonian.par_extra[6];
  extreem2          =  constants.hamiltonian.par_extra[7];

  positions_p       = &config_p->atoms.positions;


  if( !info ){

    dummy = 0.0e0;

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1; // Important

    //    fprintf( stdout, "extreems: %le, %le\n", extreem1, extreem2 );

    // first par_extra term

    i=-1;

    distance = ( positions_p->rvector[ 0 ] - extreem1 );

    sigma = dimerisation;

    //    fprintf( stdout, "i = %d, sigma = %le\n", i, sigma );

#ifdef __DEBUG__

    if( distance < EPS ){
      
      fprintf( stderr, "ERROR: the distance between particle %d and %d is not positive [%le]\n", i, i+1, distance ); 
      fflush( stderr );

      info=1;
      
    }

#endif /* __DEBUG__ */

    if( alpha_coupling *( distance -a_spacing )< 1.0e0 +sigma ){

      dummy2 = t_hopping *( 1.0e0 -alpha_coupling *( distance -a_spacing ) +sigma ); 
      
      //
      index = ELECTRON_SINGLE_INDEX( i+1, i+2 );
      
      matrix_p->matrix[ index ].z[0] = -dummy2;

      //
      index = ELECTRON_SINGLE_INDEX( i+2, i+1 );

      matrix_p->matrix[ index ].z[0] = -dummy2;

    }
    else{
	
      fprintf( stdout, "WARNING: jump at time %le\n", state_p->time );
	  	
    }
      
    /* harmonic term update */
    dummy += ( distance -a_spacing ) *( distance -a_spacing );

    // second par_extra term

    i=N_atoms-1;

    distance = ( extreem2 -positions_p->rvector[ N_atoms-1 ]);

    sigma = dimerisation *( 2*(i%2) -1 );

    //    fprintf( stdout, "i = %d, sigma = %le\n", i, sigma );

#ifdef __DEBUG__

    if( distance < EPS ){
      
      fprintf( stderr, "ERROR: the distance between particle %d and %d is not positive [%le]\n", i, i+1, distance ); 
      fflush( stderr );

      info=1;

    }

#endif /* __DEBUG__ */

    if( alpha_coupling *( distance -a_spacing )< 1.0e0 +sigma ){
      
      dummy2 = t_hopping *( 1.0e0 -alpha_coupling *( distance -a_spacing ) +sigma ); 
      
	//
      index = ELECTRON_SINGLE_INDEX( i+1, i+2 );

      matrix_p->matrix[ index ].z[0] = -dummy2;

	//
      index = ELECTRON_SINGLE_INDEX( i+2, i+1 );

      matrix_p->matrix[ index ].z[0] = -dummy2;

    }
    else{
      
      fprintf( stdout, "WARNING: jump at time %le\n", state_p->time );
	  	
    }
      
    /* harmonic term update */
    dummy += ( distance -a_spacing ) *( distance -a_spacing );

    // main loop 

    for( i=0;i<N_atoms-1;i++){
      
      distance = ( positions_p->rvector[ i+1 ] - positions_p->rvector[ i ] );
      
      sigma = dimerisation *( 2*(i%2) -1 );

      //      fprintf( stdout, "i = %d, sigma = %le\n", i, sigma );


#ifdef __DEBUG__

      if( distance < EPS ){

	fprintf( stderr, "ERROR: the distance between particle %d and %d is not positive [%le]\n", i, i+1, distance ); 
	fflush( stderr );

	info=1;

	break;

      }

#endif /* __DEBUG__ */

      //      fprintf( stdout, "i = %d, distance = %le, a_spacing = %le\n", i, distance, a_spacing );
      
      if( alpha_coupling *( distance -a_spacing )< 1.0e0 +sigma ){

	dummy2 = t_hopping *( 1.0e0 -alpha_coupling *( distance -a_spacing ) +sigma ); 

	//
	index = ELECTRON_SINGLE_INDEX( i+1, i+2 );

	matrix_p->matrix[ index ].z[0] = -dummy2;

	//
	index = ELECTRON_SINGLE_INDEX( i+2, i+1 );

	matrix_p->matrix[ index ].z[0] = -dummy2;

      }
      else{
	
	fprintf( stdout, "WARNING: jump at time %le\n", state_p->time );
	  	
      }
      
      /* harmonic term update */
      dummy += ( distance -a_spacing ) *( distance -a_spacing );

    } /* i loop */


    if( positions_p->rvector[ 0 ] - extreem1 < EPS || extreem2 - positions_p->rvector[ N_atoms -1 ] < EPS ){

      fprintf( stderr, "ERROR: invalid boundary conditions is not positive \n" ); 
      fflush( stderr );
	
      info=1;

    }


    /* Classical part */
    for( i=0; i<N_atoms+2; i++ ){
      
      matrix_p->matrix[ ELECTRON_SINGLE_INDEX( i, i ) ].z[0] = 0.5e0 *k_harmonic *dummy; //WARNING: diagonal only 
      
    }
    
  }
  

  return info;

}

//------------------------------------------

/* F_matrix_single_CHAIN_MOD2 update */

int F_matrix_single_CHAIN_MOD2_update( const constants constants, state_p state_p, config_p config_p, int r, matrix_p matrix_p ){

  /* constants */
  int        N_electrons;
  int        N_atoms;
  double     a_spacing;
  //  double     delta_energy;
  double     k_harmonic;
  double     t_hopping;
  double     alpha_coupling;
  double     dimerisation;
  /* state */
  rvector_p  positions_p;
  /* dummies */
  int        index;
  int        i;
  double     distance;
  double     dummy;
  int        info=0;
  double     extreem1;
  double     extreem2;
  double     sigma;


  N_electrons       =  constants.N_electrons;
  N_atoms           =  constants.N_atoms;
  a_spacing         =  constants.hamiltonian.par_extra[0];
  // delta_energy      =  constants.hamiltonian.par_extra[1];
  k_harmonic        =  constants.hamiltonian.par_extra[2]/(double) N_electrons; // WARNING: rescaling
  t_hopping         =  constants.hamiltonian.par_extra[3];
  alpha_coupling    =  constants.hamiltonian.par_extra[4];
  dimerisation      =  constants.hamiltonian.par_extra[5];
  extreem1          =  constants.hamiltonian.par_extra[6];
  extreem2          =  constants.hamiltonian.par_extra[7];

  positions_p       = &config_p->atoms.positions;


  if( !info ){

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1; // Important!

    dummy = t_hopping *alpha_coupling;
    
    //
    if( r > 0  ){
      
      distance = ( positions_p->rvector[ r ] - positions_p->rvector[ r-1 ] );
    }      
    else{

      distance = ( positions_p->rvector[ 0 ] - extreem1 );

    }

    sigma = dimerisation *( 1 - 2*(r%2) );

    if( alpha_coupling *( distance -a_spacing ) < 1.0e0 +sigma ){
	
      //
      index = ELECTRON_SINGLE_INDEX( r, r+1 );
      
      matrix_p->matrix[ index ].z[0] = -dummy;
	
      //
      index = ELECTRON_SINGLE_INDEX( r+1, r );
	
      matrix_p->matrix[ index ].z[0] = -dummy;
	
    }

    //  
    if( r < N_atoms-1  ){
      
      distance = ( positions_p->rvector[ r+1 ] - positions_p->rvector[ r ] );

    }      
    else{

      distance = ( extreem2 - positions_p->rvector[ r ] );

    }

    sigma = dimerisation *( 2*(r%2) -1 );
    
    if( alpha_coupling *( distance -a_spacing ) < 1.0e0 +sigma ){
	
      //
      index = ELECTRON_SINGLE_INDEX( r+1, r+2 );
      
      matrix_p->matrix[ index ].z[0] = dummy;
      
      //
      index = ELECTRON_SINGLE_INDEX( r+2, r+1 );
      
      matrix_p->matrix[ index ].z[0] = dummy;
      
    }

 

    if( r==0 ){

      dummy = k_harmonic *( positions_p->rvector[ r+1 ] -2.0e0 *positions_p->rvector[ r ] +extreem1 );

    }
    else if( r == N_atoms-1 ){


      dummy = k_harmonic *( extreem2 -2.0e0 *positions_p->rvector[ r ] +positions_p->rvector[ r-1 ] );

    }
    else{
      
      dummy = k_harmonic *( positions_p->rvector[ r+1 ] - 2.0e0 *positions_p->rvector[ r ] +positions_p->rvector[ r-1 ] );
      
    }

    for( i=0;i<N_atoms+2;i++){

      index = ELECTRON_SINGLE_INDEX( i, i );
	    
      matrix_p->matrix[ index ].z[0] = dummy;
	
    }
      
  }
  

  return info;

}

//------------------------------------------

/* K_matrix_single_CHAIN_MOD2 update */

int K_matrix_single_CHAIN_MOD2_update( const constants constants, state_p state_p, config_p config_p, int r, int s, matrix_p matrix_p ){

  /* constants */
  int        N_electrons;
  int        N_atoms;
  //  double     delta_energy;
  double     k_harmonic;
  //  double     t_hopping;
  //double     alpha_coupling;
  /* state */

  /* dummies */
  int        i;
  int        index;
  int        info=0;


  N_electrons     =  constants.N_electrons;
  N_atoms         =  constants.N_atoms;
  // delta_energy    =  constants.hamiltonian.par_extra[1];
  k_harmonic      =  constants.hamiltonian.par_extra[2]/(double) N_electrons; // WARNING: rescaling
  // t_hopping       =  constants.hamiltonian.par_extra[3];
  // alpha_coupling  =  constants.hamiltonian.par_extra[4];


  if( !info ){

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1; // Important
    
    for( i=0;i<N_atoms+2;i++){

      index = ELECTRON_SINGLE_INDEX( i, i );

      //
      if( s == r+1 || s == r-1 ){
	  
	matrix_p->matrix[ index ].z[0] = -k_harmonic;
	  
      }

      //
      if( s == r ){
	    
	matrix_p->matrix[ index ].z[0] = 2.0e0 *k_harmonic;
	  
      }

    } /* end i loop*/
	
  }


  return info;

}

//------------------------------------------
//------------------------------------------

/* utilities */

//------------------------------------------

int Hamiltonian_single_CHAIN_MOD2_parameters_check( const constants constants ){

  /*constants */
  int     sdim;
  int     N_atoms;
  int     N_levels_single;
  /* dummies */
  int info=0;

  
  sdim             =  constants.spacial_dimension; 
  N_atoms          = constants.N_atoms;
  N_levels_single  = constants.N_levels_single;


#ifdef __DEBUG__

  fprintf( stdout, "#----------\n" );
  fprintf( stdout, "HAMILTONIAN_SINGLE: using CHAIN_MOD2\n" );

#endif /* __DEBUG__ */


  if( sdim != 1 ){

    fprintf( stderr, "ERROR: this Hamiltonian_single works only for sdim = 1.\n" );
    fflush( stderr );

    info=1;

  }


  if( N_levels_single != N_atoms +2 ){

    fprintf( stderr, "ERROR: this Hamiltonian_single works only for N_levels_single = N_atoms +2.\n" );
    fflush( stderr );

    info=1;

  }


  if( N_atoms < 2 ){

    fprintf( stderr, "ERROR: this Hamiltonian_single works only for N_atoms > 1.\n" );
    fflush( stderr );

    info=1;

  }


  if( constants.hamiltonian.N_par_extra != 8 ){

    fprintf( stderr, "ERROR: this Hamiltonian_single works only with 8 parameters: a_spacing, delta_energy, k_harmonic, t_hopping, alpha_coupling, dimerisation, extreme1, extreme2\n" );
    fflush( stderr );

    info=1;

  }


#ifdef __DEBUG__

  if( !info ){

    fprintf( stdout, "# a_spacing      = %le\n", constants.hamiltonian.par_extra[ 0 ] );
    fprintf( stdout, "# delta_energy   = %le\n", constants.hamiltonian.par_extra[ 1 ] );
    fprintf( stdout, "# k_harmonic     = %le\n", constants.hamiltonian.par_extra[ 2 ] );
    fprintf( stdout, "# t_hopping      = %le\n", constants.hamiltonian.par_extra[ 3 ] );
    fprintf( stdout, "# alpha_coupling = %le\n", constants.hamiltonian.par_extra[ 4 ] );
    fprintf( stdout, "# dimerisation   = %le\n", constants.hamiltonian.par_extra[ 5 ] );
    fprintf( stdout, "# extreem1       = %le\n", constants.hamiltonian.par_extra[ 6 ] );
    fprintf( stdout, "# extreem2       = %le\n", constants.hamiltonian.par_extra[ 7 ] );
    fprintf( stdout, "#----------\n" );

  }

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------

/* H_matrix_single_CHAIN_MOD2_check */

int H_matrix_single_CHAIN_MOD2_check( const constants constants, rvector_p positions_p, matrix_p matrix_p ){


  /* constants */
  int        N_electrons;
  int        N_atoms;
  double     a_spacing;
  //  double     delta_energy;
  double     k_harmonic;
  double     t_hopping;
  double     alpha_coupling;
  double     dimerisation;
  /* dummies */
  int        index;
  int        i;
  int        info=0;
  double     distance;
  double     dummy;
  double     dummy2;
  double     extreem1;
  double     extreem2;
  double     sigma;


  N_electrons       =  constants.N_electrons;
  N_atoms           =  constants.N_atoms;
  a_spacing         =  constants.hamiltonian.par_extra[0];
  // delta_energy      =  constants.hamiltonian.par_extra[1];
  k_harmonic        =  constants.hamiltonian.par_extra[2]/(double) N_electrons; // WARNING: rescaling
  t_hopping         =  constants.hamiltonian.par_extra[3];
  alpha_coupling    =  constants.hamiltonian.par_extra[4];
  dimerisation      =  constants.hamiltonian.par_extra[5];
  extreem1          =  constants.hamiltonian.par_extra[6];
  extreem2          =  constants.hamiltonian.par_extra[7];

  if( !info ){

    dummy = 0.0e0;
    
    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1; // Important

    //    fprintf( stdout, "extreems: %le, %le [chk]\n", extreem1, extreem2 );

    // first par_extra term

    i=-1;

    distance = ( positions_p->rvector[ 0 ] - extreem1 );

    sigma = dimerisation;

    //      fprintf( stdout, "i = %d, sigma = %le\n", i, sigma );

#ifdef __DEBUG__

    if( distance < EPS ){
      
      fprintf( stderr, "ERROR: the distance between particle %d and %d is not positive [%le]\n", i, i+1, distance ); 
      fflush( stderr );

      info=1;
      
    }

#endif /* __DEBUG__ */

    if( alpha_coupling *( distance -a_spacing )< 1.0e0 +sigma ){

      dummy2 = t_hopping *( 1.0e0 -alpha_coupling *( distance -a_spacing ) +sigma ); 
      
      //
      index = ELECTRON_SINGLE_INDEX( i+1, i+2 );
      
      matrix_p->matrix[ index ].z[0] = -dummy2;

      //
      index = ELECTRON_SINGLE_INDEX( i+2, i+1 );

      matrix_p->matrix[ index ].z[0] = -dummy2;

    }

      
    /* harmonic term update */
    dummy += ( distance -a_spacing ) *( distance -a_spacing );

    // first par_extra term

    i=N_atoms-1;

    distance = ( extreem2 -positions_p->rvector[ N_atoms-1 ]);

    sigma = dimerisation *( 2*(i%2) -1 );

    //      fprintf( stdout, "i = %d, sigma = %le\n", i, sigma );

#ifdef __DEBUG__

    if( distance < EPS ){
      
      fprintf( stderr, "ERROR: the distance between particle %d and %d is not positive [%le]\n", i, i+1, distance ); 
      fflush( stderr );

      info=1;

    }

#endif /* __DEBUG__ */

    if( alpha_coupling *( distance -a_spacing )< 1.0e0 +sigma ){
      
      dummy2 = t_hopping *( 1.0e0 -alpha_coupling *( distance -a_spacing ) +sigma ); 
      
	//
      index = ELECTRON_SINGLE_INDEX( i+1, i+2 );

      matrix_p->matrix[ index ].z[0] = -dummy2;

	//
      index = ELECTRON_SINGLE_INDEX( i+2, i+1 );

      matrix_p->matrix[ index ].z[0] = -dummy2;

    }
      
    /* harmonic term update */
    dummy += ( distance -a_spacing ) *( distance -a_spacing );

    // main loop 

    for( i=0;i<N_atoms-1;i++){
      
      distance = ( positions_p->rvector[ i+1 ] - positions_p->rvector[ i ] );

      sigma = dimerisation *( 2*(i%2) -1 );

      //      fprintf( stdout, "i = %d, sigma = %le\n", i, sigma );


#ifdef __DEBUG__

      if( distance < EPS ){

	fprintf( stderr, "ERROR: the distance between particle %d and %d is not positive [chk, %le]\n", i, i+1, distance ); 
	fflush( stderr );

	info=1;

	break;

      }

#endif /* __DEBUG__ */

      //      fprintf( stdout, "i = %d, distance = %le, a_spacing = %le\n", i, distance, a_spacing );
      
      if( alpha_coupling *( distance -a_spacing )< 1.0e0 +sigma ){

	dummy2 = t_hopping *( 1.0e0 -alpha_coupling *( distance -a_spacing ) +sigma );

	//
	index = ELECTRON_SINGLE_INDEX( i+1, i+2 );

	matrix_p->matrix[ index ].z[0] = -dummy2;

	//
	index = ELECTRON_SINGLE_INDEX( i+2, i+1 );

	matrix_p->matrix[ index ].z[0] = -dummy2;

      }
      
      /* harmonic term update */
      dummy += ( distance -a_spacing ) *( distance -a_spacing );

    } /* i loop */
    

    if( positions_p->rvector[ 0 ] - extreem1 < EPS || extreem2 - positions_p->rvector[ N_atoms -1 ] < EPS ){

      fprintf( stderr, "ERROR: invalid boundary conditions is not positive \n" ); 
      fflush( stderr );
	
      info=1;

    }


    /* Classical part */
    for( i=0; i<N_atoms+2; i++ ){
      
      matrix_p->matrix[ ELECTRON_SINGLE_INDEX( i, i ) ].z[0] = 0.5e0 *k_harmonic *dummy; //WARNING: diagonal only 
      
    }
    
  }

    
  return info;

}

//------------------------------------------

/* classical_dipole_single_CHAIN_MOD2_update */

int classical_dipole_single_CHAIN_MOD2_update( const constants constants, state_p state_p, config_p config_p, int comp, matrix_p matrix_p ){


  /* constants */
  int        N_atoms;
  /* state */
  rvector_p  positions_p;
  rvector_p  centre_of_mass_p;
  /* dummies */
  int        index;
  int        i;
  int        info=0;
  double     extreem1;
  double     extreem2;


  N_atoms           =  constants.N_atoms;
  extreem1          =  constants.hamiltonian.par_extra[6];
  extreem2          =  constants.hamiltonian.par_extra[7];

  positions_p       = &(config_p->atoms.positions);
  centre_of_mass_p  = &(config_p->atoms.centre_of_mass);


  if( !info ){

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1; // Important

    //    fprintf( stdout, "extreems: %le, %le [chk]\n", extreem1, extreem2 );

    // first par_extra term

    index = ELECTRON_SINGLE_INDEX( 0, 0 );
      
    matrix_p->matrix[ index ].z[0] = extreem1 -( centre_of_mass_p->rvector[ comp ] );

      
    // second par_extra term

    index = ELECTRON_SINGLE_INDEX( N_atoms +1, N_atoms +1 );
      
    matrix_p->matrix[ index ].z[0] = extreem2 -( centre_of_mass_p->rvector[ comp ] );


    // main loop 

    for( i=0;i<N_atoms;i++){
      
      index = ELECTRON_SINGLE_INDEX( i+1, i+1 );
 
      matrix_p->matrix[ index ].z[0] = positions_p->rvector[ i ] -( centre_of_mass_p->rvector[ comp ] );


    } /* i loop */
    

    if( positions_p->rvector[ 0 ] - extreem1 < EPS || extreem2 - positions_p->rvector[ N_atoms -1 ] < EPS ){

      fprintf( stderr, "ERROR: invalid boundary conditions is not positive \n" ); 
      fflush( stderr );
	
      info=1;

    }


  }

    
  return info;

}

//------------------------------------------
