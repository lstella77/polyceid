
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
#include "PolyCEID_verlet.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int velocity_verlet_ions( const constants constants, state_p state_p, config_p config_tmp_p, config_p config_def_p, double ratio ){

  /* constants */
  int       N_coor;
  int       N_chain;
  double    dt;
  /* state */
  rvector_p positions_p;
  rvector_p momenta_p;
  rvector_p forces_cons_p;
  rvector_p masses_p;
  /* dummies */
  int i;
  int sdim;
  int info=0;


  N_coor             =  constants.N_coor;
  N_chain            =  constants.N_chain;
  sdim               =  constants.spacial_dimension;
  dt                 =  ratio *constants.dt; // WARNING: notice the use of ratio here
  masses_p           = &config_def_p->atoms.masses;

  positions_p        = &(config_tmp_p->atoms.positions);
  momenta_p          = &(config_tmp_p->atoms.momenta);
  forces_cons_p      = &(config_tmp_p->atoms.forces_cons);


  /* Nose-Hoover chain --- first half */
  if( !info && N_chain ){

    // fprintf( stdout, "DOING: Nose-Hoover chain --- first half\n" ); 

    if( NOSE_HOOVER_CHAIN( constants, *state_p, *config_tmp_p, ratio )  ) info=1;

  } /* end N_chain conditional */


  /* Verlet --- first half */
  if( !info ){

    for( i=0; i<N_coor; i++ ){

      if( config_def_p->atoms.mask.ivector[ i ] ){

        momenta_p->rvector[ i ]   += 0.5e0 *( forces_cons_p->rvector[ i ] ) *dt;

        positions_p->rvector[ i ] += ( momenta_p->rvector[ i ] ) *dt /( masses_p->rvector[ i ] ); // WARNING: it's using the updated momenta! 
      
      }

    }

  }


  /* Verlet --- second half */
  if( !info ){

    if( DISTANCES_UPDATE( constants, *state_p, *config_tmp_p ) )        info=1;

#ifdef __K_MATRIX_UPDATE__ 

    if( HAMILTONIAN_MANY_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1; 

#else /* __K_MATRIX_UPDATE__ */

    if( MATRIX_ARRAY_COPY( config_tmp_p->electrons.K_matrix, config_def_p->electrons.K_matrix ) ) info=1;

    if( HAMILTONIAN_MANY_AUX_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1; 

#endif /* __K_MATRIX_UPDATE__ */

  }

  if( !info ){

    for( i=0; i<N_coor; i++ ){

      if( config_def_p->atoms.mask.ivector[ i ] ){

        momenta_p->rvector[ i ]   += 0.5e0 *( forces_cons_p->rvector[ i ] ) *dt; // WARNING: those are the new forces

      }  

    }

  }


  /* Nose-Hoover chain --- second half */
  if( !info && N_chain ){

    /* Thermostat forces update */
    if( FORCES_THERMOSTAT_UPDATE( constants, *state_p, *config_tmp_p ) ) info=1;

    // fprintf( stdout, "DOING: Nose-Hoover chain --- second half\n" ); 

    if( NOSE_HOOVER_CHAIN( constants, *state_p, *config_tmp_p, ratio )  ) info=1;

  }        

  // Apply PBC, in the case
  if( constants.flag_periodic_boundary_condition ){
                
    for( i=0; i<N_coor; i++ ){

       positions_p->rvector[ i ] = PBC_DIFF( positions_p->rvector[ i ], constants.cell_dim.rvector[ N_coor %sdim ] );

    }   

  }               



  return info;

}

//------------------------------------------

//BUGFIX: this was originally TMP
int nose_hoover_chain( const constants constants, state_p state_p, config_p config_p, double ratio ){

  /* parameters */
  const int    N_yosh =  5;
  const double w[5] = { +0.414490771794376, +0.414490771794376, -0.657963087177503, +0.414490771794376, +0.414490771794376 };
  //double w[5];
  /* constants */
  int        N_atoms;
  int        N_coor;
  int        N_chain;
  int        N_resn;
  double     dt;
  rvector    masses;
  rvector    thermostat_masses;
  double     temperature;
  /* state */
  rvector_p  momenta_p;
  rvector_p  positions_thermostat_p;
  rvector_p  momenta_thermostat_p;
  rvector_p  forces_thermostat_p;
  /* dummies */
  double      wdt[5];
  double     wdt2[5];
  double     wdt4[5];
  double     wdt8[5];
  int        i;
  int        sdim;
  int        i_resn;
  int        i_yosh;
  int        i_atoms;
  int        i_coor;
  double     scale;
  double     ekin2;
  double     aa;
  int        info=0;


  N_atoms                    =  constants.N_atoms;
  N_coor                     =  constants.N_coor;
  N_chain                    =  constants.N_chain;
  N_resn                     =  constants.N_chain_steps;
  dt                         =  ratio *constants.dt; // WARNING: notice the use of ratio here
  masses                     =  config_p->atoms.masses;
  thermostat_masses          =  constants.thermostat_masses;
  temperature                =  constants.temperature;
  sdim                       =  constants.spacial_dimension;

  momenta_p                  = &(config_p->atoms.momenta);
  positions_thermostat_p     = &(config_p->thermostat.positions);
  momenta_thermostat_p       = &(config_p->thermostat.momenta);
  forces_thermostat_p        = &(config_p->thermostat.forces);


  /*
    w[0] = w[1] = w[3] = w[4] =  1.0e0/( 4.0e0 -pow( 4.0e0, ( 1.0e0 /3.0e0 ) ) );
    w[2] = 1.0e0 -4.0e0 *w[0];
  */

  
  for( i_yosh=0; i_yosh < N_yosh; i_yosh++ ){

     wdt[ i_yosh ] =            w[ i_yosh ] *dt /( (double) N_resn );
    wdt2[ i_yosh ] =  0.5e0  *wdt[ i_yosh ];
    wdt4[ i_yosh ] =  0.5e0 *wdt2[ i_yosh ];
    wdt8[ i_yosh ] =  0.5e0 *wdt4[ i_yosh ];
    
  }
  
  
  /* inistialise scale */
  scale = 1.0e0;

  // fprintf( stdout, "scale = %le [1]\n", scale );

  /* compute twice the kinetic energy of the system */

  /* WARNING: important! */
  ekin2 = 0.0e0;

  for( i_atoms=0; i_atoms<N_atoms; i_atoms++ ){

    for( i_coor=0; i_coor<sdim; i_coor++ ){

        ekin2 += ( momenta_p->rvector[ i_coor +sdim *i_atoms ] ) *
	         ( momenta_p->rvector[ i_coor +sdim *i_atoms ] ) /( masses.rvector[ i_coor +sdim *i_atoms ] );

    } /* end coor loop */
						     
  } /* end i_atoms loop */
						     

  // fprintf( stdout, "ekin2 = %le [1]\n", ekin2  );

  /*
    fprintf( stdout, "forces thermostat [1]\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *forces_thermostat_p ) ) info=1;
  */

  /* compute the first component of the thermkostat forces */
  //forces_thermostat_p->rvector[ 0 ] = ( ekin2 -( (double) N_coor ) *temperature );
  
  /* compute the othe components */
  /*
  for( i=1; i<N_chain; i++ ){ //WARNING: i>=1
  
    forces_thermostat_p->rvector[ i ] = ( ( momenta_thermostat_p->rvector[ i -1 ] ) *
                                              ( momenta_thermostat_p->rvector[ i -1 ] ) /( thermostat_masses.rvector[ i -1 ]) -temperature );
   
  }
  */

  /*
    fprintf( stdout, "forces thermostat [2]\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *forces_thermostat_p ) ) info=1;
  */

  /* Multiple time-step procedure */
  for( i_resn=0; i_resn < N_resn; i_resn++ ){

    // fprintf( stdout, "i_resn = %d over %d\n", i_resn, N_resn );

    for( i_yosh=0; i_yosh < N_yosh; i_yosh++ ){

      // fprintf( stdout, "i_yosh = %d over %d\n", i_yosh, N_yosh );

      /* Update thermostat MOMENTA */

      /* Update last component of the MOMENTA */
      momenta_thermostat_p->rvector[ N_chain -1 ] += ( forces_thermostat_p->rvector[ N_chain -1 ] ) *wdt4[ i_yosh ];
      
      /* Other components in reverse order */
      for( i=0; i<(N_chain -1); i++ ){ //WARNING: i<=N_chain-2!

        // fprintf( stdout, "passing here [1]\n" );
	
	aa=exp( -wdt8[ i_yosh ] * ( momenta_thermostat_p->rvector[ N_chain -1 -i ] ) /( thermostat_masses.rvector[ N_chain -1 -i ] ) );

	momenta_thermostat_p->rvector[ N_chain -2 -i ] = ( momenta_thermostat_p->rvector[ N_chain -2 -i ] ) *aa *aa +
	                                                     wdt4[ i_yosh ] *( forces_thermostat_p->rvector[ N_chain -2 -i ] ) *aa;

      } /* end i loop */       

      /* update ionic MOMENTA */
      aa=exp( -wdt2[ i_yosh ] * ( momenta_thermostat_p->rvector[ 0 ] ) /( thermostat_masses.rvector[ 0 ] ) );

      scale *= aa;

      // fprintf( stdout, "scale = %le [2]\n", scale );

      /* compute the first component of the thermkostat forces */
      forces_thermostat_p->rvector[ 0 ] = ( scale *scale *ekin2 - ( (double) N_coor ) *temperature );
  
  
      /* Update thermosta positions */ 
      for( i=0; i<N_chain; i++ ){

        positions_thermostat_p->rvector[ i ] += ( momenta_thermostat_p->rvector[ i ] ) /( thermostat_masses.rvector[ i ] ) *wdt2[ i_yosh ];

      } /* end i loop */       


      /* Update the thermostat MOMENTA and forces */
      for( i=0; i<N_chain-1; i++ ){ //WARNING: i<=N_chain-2

        // fprintf( stdout, "passing here [2]\n" );
	
        aa = exp( -wdt8[ i_yosh ] * ( momenta_thermostat_p->rvector[ i +1 ] ) /( thermostat_masses.rvector[ i +1 ] ) );

	momenta_thermostat_p->rvector[ i ] = ( momenta_thermostat_p->rvector[ i ] ) *aa *aa +
	                                         wdt4[ i_yosh ] *( forces_thermostat_p->rvector[ i ] ) *aa;

        forces_thermostat_p->rvector[ i +1 ] = ( ( momenta_thermostat_p->rvector[ i ] ) *
	                                             ( momenta_thermostat_p->rvector[ i ] ) /( thermostat_masses.rvector[ i ] ) -temperature );
  
	
      } /* end i loop */

      /* Update last component of the MOMENTA*/
      momenta_thermostat_p->rvector[ N_chain -1 ] += ( forces_thermostat_p->rvector[ N_chain -1 ] ) *wdt4[ i_yosh ];
      


    } /* end yosh loop */        

  } /* end i_resn loop */


  /* Update ionic MOMENTA */
  for( i=0; i<N_coor; i++ ){

    momenta_p->rvector[ i ] *= scale;

  }        


  return info;

}  

//------------------------------------------
