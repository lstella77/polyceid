
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
#include "PolyCEID_energies.h"



/*********************
  FUNCTIONS & MACROS
*********************/

int energies_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             N_chain;
  /* state */
  rvector_p       kinetic_energy_system_p;
  rvector_p       potential_energy_system_p;
  rvector_p       kinetic_energy_thermostat_p;
  rvector_p       potential_energy_thermostat_p;
  double*         kinetic_energy_system_correction_p;
  double*         potential_energy_system_correction_p;
  double*         total_energy_system_p;
  double*         pseudo_energy_system_p;
  /* dummies */
  int i;
  int info=0;


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "DOING: energies_update\n");

#endif /* __DEBUG_PLUS__ */


  N_chain                              =  constants.N_chain;

  kinetic_energy_system_p              = &(state_p->observables.kinetic_energy_system);
  potential_energy_system_p            = &(state_p->observables.potential_energy_system);
  
  kinetic_energy_thermostat_p          = &(state_p->observables.kinetic_energy_thermostat);
  potential_energy_thermostat_p        = &(state_p->observables.potential_energy_thermostat);
  
  kinetic_energy_system_correction_p   = &(state_p->observables.kinetic_energy_system_correction);
  potential_energy_system_correction_p = &(state_p->observables.potential_energy_system_correction);
  
  total_energy_system_p                = &(state_p->observables.total_energy_system);
  pseudo_energy_system_p               = &(state_p->observables.pseudo_energy_system);


  if( !info ){

    /* kinetic energy */
    if( KINETIC_ENERGY_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }

  kinetic_energy_system_p->rvector[ 3 ]  = kinetic_energy_system_p->rvector[ 0 ] +
                                           kinetic_energy_system_p->rvector[ 1 ] +
                                           kinetic_energy_system_p->rvector[ 2 ] +
                                          *kinetic_energy_system_correction_p;


  if( !info && N_chain ){

    /* kinetic energy thermostat */
    if( KINETIC_ENERGY_THERMOSTAT_UPDATE( constants, *state_p, *config_p ) ) info=1;
    
    /* Important! */
    kinetic_energy_thermostat_p->rvector[ N_chain ] = ZERO;

    for( i=0; i<N_chain; i++ ){

      kinetic_energy_thermostat_p->rvector[ N_chain ] +=  kinetic_energy_thermostat_p->rvector[ i ];

    }        
   
  }        


  if( !info ){

    /* potential energy */
    if( POTENTIAL_ENERGY_UPDATE( constants, *state_p, *config_p ) ) info=1;

  }

  potential_energy_system_p->rvector[ 3 ] = potential_energy_system_p->rvector[ 0 ] +
                                            potential_energy_system_p->rvector[ 1 ] +
                                            potential_energy_system_p->rvector[ 2 ] +
                                           *potential_energy_system_correction_p;


  if( !info && N_chain ){

    /* potential energy thermostat */
    if( POTENTIAL_ENERGY_THERMOSTAT_UPDATE( constants, *state_p, *config_p ) ) info=1;
    
    /* Important! */
    potential_energy_thermostat_p->rvector[ N_chain ] = ZERO;

    for( i=0; i<N_chain; i++ ){

      potential_energy_thermostat_p->rvector[ N_chain ] +=  potential_energy_thermostat_p->rvector[ i ];

    }        
   
  }        


  if( !info ){

    /* total energy */
    *total_energy_system_p = ( kinetic_energy_system_p->rvector[ 3 ] ) +
                             ( potential_energy_system_p->rvector[ 3 ] );

  }


  if( !info && N_chain ){

    /* pseudo energy */
    *pseudo_energy_system_p = *total_energy_system_p +
                             ( kinetic_energy_thermostat_p->rvector[ N_chain ] ) +
                             ( potential_energy_thermostat_p->rvector[ N_chain ] );

  }


  return info;

}

//------------------------------------------
//------------------------------------------

int kinetic_energy_update( const constants constants, state_p state_p, config_p config_p ){

  // WARNING: it's probably a good idea to restric the summations over N_coor_red;


  /* constants */
  int             N_atoms;
  int             sdim;
  /* state */
  rvector_p       momenta_p;
  rvector_p       masses_p;
  rvector_p       kinetic_energy_atoms_p;
  rvector_p       kinetic_energy_system_p;
  matrix_array_p  mu01_p;
  matrix_array_p  mu02_p;
  /* dummies */
  int             index;
  int             i_atoms, i_coor;
  double          dummy1, dummy2, dummy3;
  int             info=0;


  N_atoms                        =  constants.N_atoms;
  sdim                           =  constants.spacial_dimension; 
  masses_p                       = &config_p->atoms.masses;

  momenta_p                      = &(config_p->atoms.momenta);
  kinetic_energy_atoms_p         = &(state_p->observables.kinetic_energy_atoms);
  kinetic_energy_system_p        = &(state_p->observables.kinetic_energy_system);
  mu01_p                         = &(config_p->electrons.mu01);
  mu02_p                         = &(config_p->electrons.mu02);


  /* set to zero */
  if( RVECTOR_ZERO( *kinetic_energy_system_p ) ) info=1; 

  for( i_atoms=0; i_atoms<N_atoms; i_atoms++ ){

    /* set to zero */
    dummy3 = dummy2 = dummy1 = ZERO;

    for( i_coor=0; i_coor<sdim; i_coor++ ){

      index = i_coor +sdim *i_atoms;

      if( masses_p->rvector[ index ] > EPS ){ 

        dummy1 += ( momenta_p->rvector[ index ] ) *( momenta_p->rvector[ index ] ) /( masses_p->rvector[ index ] );

#ifndef __NO_CEID__

        dummy2 += ( momenta_p->rvector[ index ] ) *REAL( MATRIX_TRACE( mu01_p->array[ index ] ) ) /( masses_p->rvector[ index ] );

        dummy3 +=  REAL( MATRIX_TRACE( mu02_p->array[ COORDINATE_INDEX( index, index ) ] ) ) /( masses_p->rvector[ index ] );

#endif /* __NO_CEID__ */

      }

    } /* end coor loop */


    kinetic_energy_atoms_p->rvector[ 3 *i_atoms ]    = ONEO2 *dummy1;

    kinetic_energy_atoms_p->rvector[ 3 *i_atoms +1 ] = dummy2;

    kinetic_energy_atoms_p->rvector[ 3 *i_atoms +2 ] = ONEO2 *dummy3;


    kinetic_energy_system_p->rvector[ 0 ] += kinetic_energy_atoms_p->rvector[ 3 *i_atoms ];

    kinetic_energy_system_p->rvector[ 1 ] += kinetic_energy_atoms_p->rvector[ 3 *i_atoms +1 ];

    kinetic_energy_system_p->rvector[ 2 ] += kinetic_energy_atoms_p->rvector[ 3 *i_atoms +2 ];


  } /* end atoms loop */ 


  return info;

}


//------------------------------------------

int potential_energy_update( const constants constants, state_p state_p, config_p config_p ){

  // WARNING: it's probably a good idea to restric the summations over N_coor_red;


  /* constants */
  int             N_coor;
  /* state */
  rvector_p       potential_energy_system_p;
  matrix_p        mu00_p;
  matrix_array_p  mu10_p;
  matrix_array_p  mu20_p;
  matrix_p        H_matrix_p;
  matrix_array_p  F_matrix_p;
  matrix_array_p  K_matrix_p;
  matrix_p        dummy_matrix1_p;
  /* dummies */
  int             info=0;
#ifndef __NO_CEID__  
  int             i_coor1, i_coor2;
#endif /* __NO_CEID__ */  
  double          dummy1, dummy2, dummy3;


  N_coor      =  constants.N_coor;

  potential_energy_system_p        = &(state_p->observables.potential_energy_system);
  mu00_p                           = &(config_p->electrons.mu00);
  mu10_p                           = &(config_p->electrons.mu10);
  mu20_p                           = &(config_p->electrons.mu20);
  H_matrix_p                       = &(config_p->electrons.H_matrix);
  F_matrix_p                       = &(config_p->electrons.F_matrix);
  K_matrix_p                       = &(config_p->electrons.K_matrix);
  dummy_matrix1_p                  = &(state_p->dummy_matrix1);


  if( MATRIX_MATRIX_PRODUCT( *H_matrix_p, *mu00_p, *dummy_matrix1_p  ) ) info=1;

  /*
    fprintf( stdout, "H_matrix\n" );
    if( MATRIX_PRINT_PLUS( stdout, *H_matrix_p ) ) info = 1;
    fprintf( stdout, "\n" );
  */

  dummy1 = REAL( MATRIX_TRACE( *dummy_matrix1_p ) );

  /* set to zero */
  dummy3 = dummy2 = ZERO;

#ifndef __NO_CEID__ 
 
  for( i_coor1=0; i_coor1<N_coor; i_coor1++ ){

    if( MATRIX_MATRIX_PRODUCT( F_matrix_p->array[ i_coor1 ], mu10_p->array[ i_coor1 ], *dummy_matrix1_p  ) ) info=1;

    /*
      fprintf( stdout, "i_coor1 = %d\n", i_coor1 );
      
      fprintf( stdout, "F_matrix\n" );
      if( MATRIX_PRINT_PLUS( stdout, F_matrix_p->array[ i_coor1 ] ) ) info = 1;
      fprintf( stdout, "\n" );
    */

    dummy2 += REAL( MATRIX_TRACE( *dummy_matrix1_p ) );

#ifndef __ENERGY_MOD__ 

    for( i_coor2=0; i_coor2<N_coor; i_coor2++ ){

      
	//  fprintf( stdout, "i_coor1 = %d, i_coor2 = %d\n", i_coor1, i_coor2 );

	/*
	fprintf( stdout, "K_matrix\n" );
	if( MATRIX_PRINT_PLUS( stdout, K_matrix_p->array[ COORDINATE_INDEX( i_coor1, i_coor2 ) ] ) ) info = 1;
	fprintf( stdout, "\n" );

	fprintf( stdout, "mu20\n" );
	if( MATRIX_PRINT_PLUS( stdout, mu20_p->array[ COORDINATE_INDEX( i_coor1, i_coor2 ) ] ) ) info = 1;
	fprintf( stdout, "\n" );
	*/

      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ COORDINATE_INDEX( i_coor1, i_coor2 ) ], 
				 mu20_p->array[ COORDINATE_INDEX( i_coor1, i_coor2 ) ], *dummy_matrix1_p ) ) info=1;

      /*
	fprintf( stdout, "dummy_matrix1\n" );
	if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p ) ) info = 1;
	fprintf( stdout, "\n" );
      */

      /*
        fprintf( stdout, "norm K_matrix  %15.9le \n", MATRIX_NORM(  K_matrix_p->array[ COORDINATE_INDEX( i_coor1, i_coor2 ) ] ) );
        fprintf( stdout, "norm mu02      %15.9le \n", MATRIX_NORM(  mu20_p->array[ COORDINATE_INDEX( i_coor1, i_coor2 ) ] ) );
        fprintf( stdout, "norm dummy1        %15.9le \n", MATRIX_NORM( *dummy_matrix1_p ) );
      */

      //      fprintf( stdout, "dummy3: [before] = %le\n", dummy3 );

      dummy3 += REAL( MATRIX_TRACE( *dummy_matrix1_p ) );

      //      fprintf( stdout, "dummy3: [after] = %le\n", dummy3 );

    } /* end i_coor2 loop */

#else /* __ENERGY_MOD__ */

      /*
      fprintf( stdout, "K_matrix\n" );
      if( MATRIX_PRINT_PLUS( stdout, K_matrix_p->array[ COORDINATE_INDEX( i_coor1, i_coor1 ) ] ) ) info = 1;
      fprintf( stdout, "\n" );

      fprintf( stdout, "mu20\n" );
      if( MATRIX_PRINT_PLUS( stdout, mu20_p->array[ COORDINATE_INDEX( i_coor1, i_coor1 ) ] ) ) info = 1;
      fprintf( stdout, "\n" );
      */

    if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ COORDINATE_INDEX( i_coor1, i_coor1 ) ], 
				 mu20_p->array[ COORDINATE_INDEX( i_coor1, i_coor1 ) ], *dummy_matrix1_p ) ) info=1;

    /*
      fprintf( stdout, "dummy_matrix1\n" );
      if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p ) ) info = 1;
      fprintf( stdout, "\n" );
    */

    /*
      fprintf( stdout, "norm K_matrix  %15.9le \n", MATRIX_NORM(  K_matrix_p->array[ COORDINATE_INDEX( i_coor1, i_coor1 ) ] ) );
      fprintf( stdout, "norm mu02      %15.9le \n", MATRIX_NORM(  mu20_p->array[ COORDINATE_INDEX( i_coor1, i_coor1 ) ] ) );
      fprintf( stdout, "norm dummy1        %15.9le \n", MATRIX_NORM( *dummy_matrix1_p ) );
    */

    //      fprintf( stdout, "dummy3: [before] = %le\n", dummy3 );

    dummy3 += REAL( MATRIX_TRACE( *dummy_matrix1_p ) );

    //      fprintf( stdout, "dummy3: [after] = %le\n", dummy3 );

#endif /* __ENERGY_MOD__ */

  } /* end i_coor1 loop */

#endif /* __NO_CEID__ */  
  
  potential_energy_system_p->rvector[ 0 ] = dummy1;

  potential_energy_system_p->rvector[ 1 ] = -dummy2;


  potential_energy_system_p->rvector[ 2 ] = ONEO2 *dummy3;


  return info;

}

//------------------------------------------

int kinetic_energy_thermostat_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        N_chain;
  rvector    thermostat_masses;
  /* state */
  rvector_p  kinetic_energy_thermostat_p;
  rvector_p  momenta_thermostat_p;
  /* dummies */
  int i;
  int info=0;


  N_chain                      =  constants.N_chain;
  thermostat_masses            =  constants.thermostat_masses;

  kinetic_energy_thermostat_p  = &(state_p->observables.kinetic_energy_thermostat);
  momenta_thermostat_p         = &(config_p->thermostat.momenta);

  /*
    fprintf( stdout, "# thermostat_masses \n" );
    if( RVECTOR_PRINT_PLUS( stdout, thermostat_masses ) ) info=1;

    fprintf( stdout, "# momenta_thermostat [before]\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *momenta_thermostat_p ) ) info=1;
  */
  
  for( i=0; i<N_chain; i++ ){

    kinetic_energy_thermostat_p->rvector[ i ] = ONEO2 *(  momenta_thermostat_p->rvector[ i ] ) *(  momenta_thermostat_p->rvector[ i ] )
                                               /( thermostat_masses.rvector[ i ] );

  }

  /*
    fprintf( stdout, "# momenta_thermostat [after]\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *momenta_thermostat_p ) ) info=1;
  */
  
  
  return info;
  
}	

//------------------------------------------

int potential_energy_thermostat_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        N_coor;
  int        N_chain;
  double     temperature;
  /* state */
  rvector_p  potential_energy_thermostat_p;
  rvector_p  positions_thermostat_p;
  /* dummies */
  int i;
  int info=0;


  N_coor                         =  constants.N_coor;
  N_chain                        =  constants.N_chain;
  temperature                    =  constants.temperature;

  potential_energy_thermostat_p  = &(state_p->observables.potential_energy_thermostat);
  positions_thermostat_p         = &(config_p->thermostat.positions);


  /* WARNING: first term multiplied by the DOF */
  potential_energy_thermostat_p->rvector[ 0 ] = ( (double) N_coor ) *temperature *( positions_thermostat_p->rvector[ 0 ] );
  
  for( i=1; i<N_chain; i++ ){ //WARNING: starting from 1!

    potential_energy_thermostat_p->rvector[ i ] = temperature *( positions_thermostat_p->rvector[ i ] );

  }

  
  return info;
  
}	

//------------------------------------------
