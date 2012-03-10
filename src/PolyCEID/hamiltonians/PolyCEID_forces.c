
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
#include "PolyCEID_forces.h"



/*********************
  FUNCTIONS & MACROS
*********************/

/* forces_update */

int forces_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int              N_coor;
  /* state */
  rvector_p        forces_p;
  rvector_p        forces_cons_p;
  matrix_p         dummy_matrix1_p;
  matrix_array_p   F_matrix_p;
  matrix_array_p   K_matrix_p;
  matrix_p         mu00_p;
  matrix_array_p   mu10_p;
  /* dummies */
  int              i_coor;
#ifdef __FORCES_NON_CONS__  
  double           friction;
#endif /* __FORCES_NON_CONS__ */ 
  int              j_coor;
  int              K_matrix_index;
  double           dummy;
  int              info=0;


  N_coor               =  constants.N_coor;

  forces_p             = &(config_p->atoms.forces );
  forces_cons_p        = &(config_p->atoms.forces_cons );
  dummy_matrix1_p      = &(state_p->dummy_matrix1);
  F_matrix_p           = &(config_p->electrons.F_matrix);
  K_matrix_p           = &(config_p->electrons.K_matrix);
  mu00_p               = &(config_p->electrons.mu00);
  mu10_p               = &(config_p->electrons.mu10);

  
  for( i_coor=0; i_coor<N_coor; i_coor++ ){

    if( MATRIX_MATRIX_PRODUCT( *mu00_p, F_matrix_p->array[ i_coor ], *dummy_matrix1_p ) ) info=1;

    forces_cons_p->rvector[ i_coor ] = REAL( MATRIX_TRACE( *dummy_matrix1_p ) );
    
    if( !constants.flag_Ehrenfest ){

      for( j_coor=0; j_coor<N_coor; j_coor++ ){

        K_matrix_index = COORDINATE_INDEX( i_coor, j_coor );

        if( MATRIX_MATRIX_PRODUCT( mu10_p->array[ j_coor ], K_matrix_p->array[ K_matrix_index ], *dummy_matrix1_p ) ) info=1;

        dummy = REAL( MATRIX_TRACE( *dummy_matrix1_p ) );

        forces_cons_p->rvector[ i_coor ] -= dummy; //WARNING: minus sign

      } /* end j_coor loop */

    }

  } /* end i_coor loop */


  /* Compute non-consevative part */
  if( RVECTOR_COPY( *forces_p, *forces_cons_p) ) info=1;

#ifdef __FORCES_NON_CONS__

  if( constants.N_chain ){

    friction = ( config_p->thermostat.momenta.rvector[ 0 ] ) /( constants.thermostat_masses.rvector[ 0 ] );

    for( i_coor=0; i_coor< N_coor; i_coor++ ){

       forces_p->rvector[ i_coor ] -= friction *( config_p->atoms.momenta.rvector[ i_coor ] );

    }
    
  }        

#endif /* __FORCES_NON_CONS__ */


  return info;

}

//------------------------------------------
//------------------------------------------

/* forces_thermostat_update */

int forces_thermostat_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int        sdim;
  int        N_atoms;
  int        N_coor;
  int        N_chain;
  rvector    masses;
  rvector    thermostat_masses;
  double     temperature;
  /* state */
  rvector_p  momenta_p;
  rvector_p  momenta_thermostat_p;
  rvector_p  forces_thermostat_p;
  /* dummies */
  int        i;
  int        i_atoms;
  int        i_coor;
  double     ekin2;
  int        info=0;


  N_atoms                   =  constants.N_atoms;
  N_coor                    =  constants.N_coor;
  N_chain                   =  constants.N_chain;
  masses                    =  config_p->atoms.masses;
  thermostat_masses         =  constants.thermostat_masses;
  temperature               =  constants.temperature;
  sdim                      =  constants.spacial_dimension; 

  momenta_p                 = &(config_p->atoms.momenta);
  momenta_thermostat_p      = &(config_p->thermostat.momenta);
  forces_thermostat_p       = &(config_p->thermostat.forces);
  

  /* WARNING: important! */
  ekin2 = 0.0e0;
    
  for( i_atoms=0; i_atoms<N_atoms; i_atoms++ ){
	      
    for( i_coor=0; i_coor<sdim; i_coor++ ){
			  
      ekin2 += ( momenta_p->rvector[ i_coor +sdim *i_atoms ] ) *
               ( momenta_p->rvector[ i_coor +sdim *i_atoms ] ) /( masses.rvector[ i_coor +sdim *i_atoms ] );
						   
    } /* end coor loop */
						       
  } /* end i_atoms loop */
							 
	 
  // fprintf( stdout, "ekin2 = %le [1]\n", ekin2  );
							   
  /* compute the first component of the thermkostat forces */
  forces_thermostat_p->rvector[ 0 ] = ( ekin2 -( (double) N_coor ) *temperature );
							       
  /* compute the othe components */
  for( i=1; i<N_chain; i++ ){ //WARNING: i>=1
					   
    forces_thermostat_p->rvector[ i ] = ( ( momenta_thermostat_p->rvector[ i -1 ] ) *
                                              ( momenta_thermostat_p->rvector[ i -1 ] ) /( thermostat_masses.rvector[ i -1 ]) -temperature );
									       
  }


  return info;

}

//------------------------------------------
