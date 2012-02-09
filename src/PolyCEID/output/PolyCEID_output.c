
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
#include "PolyCEID_output.h"


/* private varaiables */

FILE*      geometry_fp;
FILE*      positions_fp;
FILE*      momenta_fp;
FILE*      forces_fp;
FILE*      positions_thermostat_fp;
FILE*      momenta_thermostat_fp;
FILE*      forces_thermostat_fp;
FILE*      populations_fp;
FILE*      energies_fp;
FILE*      mu_traces_fp;
FILE*      mu_norms_fp;
#ifdef __DEBUG_PLUS__
FILE*      rho_traces_fp;
FILE*      rho_norms_fp;
#endif /* __DEBUG_PLUS__ */
FILE*      adiabatic_populations_fp;
FILE*      adiabatic_projections_fp;
FILE*      projections_fp;
FILE*      adiabatic_PES_many_fp;
FILE*      adiabatic_PES_single_fp;
FILE*      single_level_populations_Ehrenfest_fp;
FILE*      single_level_populations_adiabatic_fp;
FILE*      nonadiabatic_couplings_fp;
FILE*      nonadiabatic_rates_fp;
FILE*      one_body_electronic_density_matrix_fp;
FILE*      one_body_electronic_hole_matrix_fp;
FILE*      one_body_electronic_particle_matrix_fp;
FILE*      natural_orbitals_fp;
FILE*      hole_orbitals_fp;
FILE*      particle_orbitals_fp;
FILE*      adiabatic_states_fp;
FILE*      electronic_density_states_fp;
FILE*      ionic_density_states_fp;
FILE*      dipoles_many_fp;
FILE*      dipoles_single_fp;


/*********************
  FUNCTIONS & MACROS
*********************/

int print_observables( const constants constants, const state state, const config config ){

  /* constants */
  int N_levels_many;
  int N_chain;
  /* dummies */
  int  info=0;


  N_levels_many = constants.N_levels_many;
  N_chain       = constants.N_chain;
	     

  if( constants.flag_observable_all || constants.flag_observable_geometry ){

    if( PRINT_GEOMETRY( constants, state, config ) )                     info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_positions ){

    if( PRINT_POSITIONS( constants, state, config ) )                   info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_momenta ){

    if( PRINT_MOMENTA( constants, state, config ) )                     info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_forces ){

    if( PRINT_FORCES( constants, state, config ) )                      info=1;

  }  

  if( N_chain ){

    if( constants.flag_observable_all ){

      if( PRINT_POSITIONS_THERMOSTAT( constants, state, config ) )      info=1;

    }  

    if( constants.flag_observable_all ){

      if( PRINT_MOMENTA_THERMOSTAT( constants, state, config ) )        info=1;

    }  

    if( constants.flag_observable_all ){

      if( PRINT_FORCES_THERMOSTAT( constants, state, config ) )         info=1;

    }   

  }        

  if( constants.flag_observable_all || constants.flag_observable_populations ){

    if( PRINT_POPULATIONS( constants, state, config ) )                 info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_mus ){

    if( PRINT_MU_TRACES( constants, state, config ) )                   info=1;

    if( PRINT_MU_NORMS( constants, state, config ) )                    info=1;

  }  

#ifdef __DEBUG_PLUS__

  if( constants.flag_observable_all ){

    if( PRINT_RHO_TRACES( constants, state, config ) )                  info=1;

  }  

  if( constants.flag_observable_all ){

    if( PRINT_RHO_NORMS( constants, state, config ) )                   info=1;

  }  

#endif /* __DEBUG_PLUS__ */

  if( constants.flag_observable_all || constants.flag_observable_energies ){

    if( PRINT_ENERGIES( constants, state ) )                            info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_adiabatic_populations ){

    if( PRINT_ADIABATIC_POPULATIONS( constants, state ) )               info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_projections ){

    if( PRINT_ADIABATIC_PROJECTIONS( constants, state ) )               info=1;

    if( PRINT_PROJECTIONS( state ) )                                    info=1;

  }

  if( constants.flag_observable_all || constants.flag_observable_adiabatic_pes_many ){

    if( PRINT_ADIABATIC_PES_MANY( constants, state ) )                  info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_adiabatic_pes_single ){

    if( PRINT_ADIABATIC_PES_SINGLE( constants, state ) )                info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_single_level_populations ){

    if( PRINT_SINGLE_LEVEL_POPULATIONS_EHRENFEST( constants, state ) )  info=1;

    if( PRINT_SINGLE_LEVEL_POPULATIONS_ADIABATIC( constants, state ) )  info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_nonadiabatic_couplings ){

    if( PRINT_NONADIABATIC_COUPLINGS( constants, state ) )              info=1; 

  }  

  /*
  if( constants.flag_observable_all ){

    if( PRINT_NONADIABATIC_RATES( constants, state ) )                  info=1; 

  }
  */

  if( constants.flag_observable_all || constants.flag_observable_density_matrix ){

    if( PRINT_ONE_BODY_ELECTRONIC_DENSITY_MATRIX( constants, state ) )  info=1;

    if( PRINT_NATURAL_ORBITALS( constants, state ) )                    info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_transition_matrices ){

    if( PRINT_ONE_BODY_ELECTRONIC_HOLE_MATRIX( constants, state ) )     info=1;

    if( PRINT_ONE_BODY_ELECTRONIC_PARTICLE_MATRIX( constants, state ) ) info=1;

    if( PRINT_HOLE_ORBITALS( constants, state ) )                       info=1;

    if( PRINT_PARTICLE_ORBITALS( constants, state ) )                   info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_adiabatic_states ){

    if( PRINT_ADIABATIC_STATES( constants, state ) )                    info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_electronic_density_states ){

    if( PRINT_ELECTRONIC_DENSITY_STATES( constants, state ) )           info=1;

  }  

  if( constants.flag_observable_all || constants.flag_observable_ionic_density_states ){

    if( PRINT_IONIC_DENSITY_STATES( constants, state ) )                info=1;

  }  

  if( ( N_levels_many > 1 && constants.flag_observable_all ) || constants.flag_observable_dipoles_many ){

    if( PRINT_DIPOLES_MANY( constants, state ) )                        info=1;

  }

  if( constants.flag_observable_all || constants.flag_observable_dipoles_single ){

    if( PRINT_DIPOLES_SINGLE( constants, state ) )                       info=1;

  }  


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_positions( const constants constants, const state state, const config config ){

  /* constants */
  int       N_coor;
  /* state */
  double    time;
  const rvector*  positions_p;
  /* dummies */
  int       i;
  //  double    com;
  int       info=0;


  N_coor               =  constants.N_coor;

  time                 =  state.time;
  positions_p          = &(config.atoms.positions);


  fprintf( positions_fp, "%le  ", time );

  for( i=0; i<(N_coor-1); i++ ){

    fprintf( positions_fp, DOUBLE_FORMAT"  ", positions_p->rvector[ i ] );

  }
  fprintf( positions_fp, DOUBLE_FORMAT"\n", positions_p->rvector[ i ] );

  fflush( positions_fp );


  // check only
  /*
    com=0.0e0;
    for( i=0; i<N_coor; i++ ){
    
    com += positions_p->rvector[ i ];
    
    }

    com /= (double) N_coor;
    
    fprintf( stdout, "# delta positions [observables]\n" );
    fprintf( stdout, "# time %le\n", time );

    for( i=0; i<N_coor; i++ ){
    
    fprintf( stdout, "%21.15le\n", positions_p->rvector[ i ] + positions_p->rvector[ N_coor -1 -i ] -2.0e0 *com );
    
    }
    fprintf( stdout, "\n" );
    fflush( stdout );
  */
  
  
  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_momenta( const constants constants, const state state, const config config ){

  /* constants */
  int       N_coor;
  /* state */
  double    time;
  const rvector*  momenta_p;
  /* dummies */
  int       i;
  int       info=0;


  N_coor               =  constants.N_coor;

  time                 =  state.time;
  momenta_p            = &(config.atoms.momenta);


  fprintf( momenta_fp, "%le  ", time );

  for( i=0; i<(N_coor-1); i++ ){

    fprintf( momenta_fp, DOUBLE_FORMAT"  ", momenta_p->rvector[ i ] );

  }
  fprintf( momenta_fp, DOUBLE_FORMAT"\n", momenta_p->rvector[ i ] );

  fflush( momenta_fp );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_forces( const constants constants, const state state, const config config ){

  /* constants */
  int       N_coor;
  /* state */
  double    time;
  const rvector*  forces_p;
  /* dummies */
  int       i;
  //  double    com;
  int       info=0;


  N_coor               =  constants.N_coor;

  time                 =  state.time;
  forces_p             = &(config.atoms.forces);


  fprintf( forces_fp, "%le   ", time );

  for( i=0; i<(N_coor-1); i++ ){

    fprintf( forces_fp, DOUBLE_FORMAT"  ", forces_p->rvector[ i ] );

  }
  fprintf( forces_fp, DOUBLE_FORMAT"\n", forces_p->rvector[ i ] );

  fflush( forces_fp );


  // check only
  /*
    com=0.0e0;
    for( i=0; i<N_coor; i++ ){
    
    com += forces_p->rvector[ i ];
    
    }

    com /= (double) N_coor;
    
    fprintf( stdout, "# delta forces [observables]\n" );
    fprintf( stdout, "# time %le\n", time );
    
  for( i=0; i<N_coor; i++ ){
  
  fprintf( stdout, "%21.15le\n", forces_p->rvector[ i ] + forces_p->rvector[ N_coor -1 -i ] -2.0e0 *com );
  
  }
  fprintf( stdout, "\n" );
  fflush( stdout );
  */


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_positions_thermostat( const constants constants, const state state, const config config ){

  /* constants */
  int       N_chain;
  /* state */
  double    time;
  const rvector*  positions_thermostat_p;
  /* dummies */
  int       i;
  //  double    com;
  int       info=0;


  N_chain                  =  constants.N_chain;

  time                     =  state.time;
  positions_thermostat_p  = &(config.thermostat.positions);


  fprintf( positions_thermostat_fp, "%le  ", time );

  for( i=0; i<(N_chain-1); i++ ){

    fprintf( positions_thermostat_fp, DOUBLE_FORMAT"  ", positions_thermostat_p->rvector[ i ] );

  }
  fprintf( positions_thermostat_fp, DOUBLE_FORMAT"\n", positions_thermostat_p->rvector[ i ] );

  fflush( positions_thermostat_fp );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_momenta_thermostat( const constants constants, const state state, const config config ){

  /* constants */
  int       N_chain;
  /* state */
  double    time;
  const rvector*  momenta_thermostat_p;
  /* dummies */
  int       i;
  //  double    com;
  int       info=0;


  N_chain                  =  constants.N_chain;

  time                     =  state.time;
  momenta_thermostat_p     = &(config.thermostat.momenta);


  fprintf( momenta_thermostat_fp, "%le  ", time );

  for( i=0; i<(N_chain-1); i++ ){

    fprintf( momenta_thermostat_fp, DOUBLE_FORMAT"  ", momenta_thermostat_p->rvector[ i ] );

  }
  fprintf( momenta_thermostat_fp, DOUBLE_FORMAT"\n", momenta_thermostat_p->rvector[ i ] );

  fflush( momenta_thermostat_fp );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_forces_thermostat( const constants constants, const state state, const config config ){

  /* constants */
  int       N_chain;
  /* state */
  double    time;
  const rvector*  forces_thermostat_p;
  /* dummies */
  int       i;
  //  double    com;
  int       info=0;


  N_chain                  =  constants.N_chain;

  time                     =  state.time;
  forces_thermostat_p      = &(config.thermostat.forces);


  fprintf( forces_thermostat_fp, "%le  ", time );

  for( i=0; i<(N_chain-1); i++ ){

    fprintf( forces_thermostat_fp, DOUBLE_FORMAT"  ", forces_thermostat_p->rvector[ i ] );

  }
  fprintf( forces_thermostat_fp, DOUBLE_FORMAT"\n", forces_thermostat_p->rvector[ i ] );

  fflush( forces_thermostat_fp );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_populations( const constants constants, const state state, const config config ){

  /* constants */
  int       N_levels_many;
  /* state */
  double    time;
  const matrix*  density_matrix_p;
  double    dummy;
  /* dummies */
  int       i;
  int       info=0;


  N_levels_many     =  constants.N_levels_many;

  time              =  state.time;
  density_matrix_p  = &(config.electrons.mu00);


  fprintf( populations_fp, "%le  ", time); 

  for( i=0; i<(N_levels_many-1); i++ ){

    dummy = REAL( density_matrix_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ] );

    fprintf( populations_fp, DOUBLE_FORMAT"  ", dummy );

  }

  dummy = REAL( density_matrix_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ] );

  fprintf( populations_fp, DOUBLE_FORMAT"\n", dummy );

  fflush( populations_fp );


  return info;

}

//------------------------------------------

int print_energies( const constants constants, const state state ){

  /* constants */
  int      N_chain;
  /* state */
  double   time;
  const rvector*  kinetic_energy_system_p;
  const rvector*  potential_energy_system_p;
  const rvector*  kinetic_energy_thermostat_p;
  const rvector*  potential_energy_thermostat_p;
  double   total_energy_system;
  double   pseudo_energy_system;
  /* dummies */
  int     i;
  int     info=0;


  N_chain                             =  constants.N_chain;
  
  time                                =  state.time;
  kinetic_energy_system_p             = &(state.observables.kinetic_energy_system);
  potential_energy_system_p           = &(state.observables.potential_energy_system);
  kinetic_energy_thermostat_p         = &(state.observables.kinetic_energy_thermostat);
  potential_energy_thermostat_p       = &(state.observables.potential_energy_thermostat);
  total_energy_system                 =  state.observables.total_energy_system;
  pseudo_energy_system                =  state.observables.pseudo_energy_system;


  fprintf( energies_fp, "%le  ", time  );

  for( i=0; i<4; i++ ){

    fprintf( energies_fp, DOUBLE_FORMAT"  ", kinetic_energy_system_p->rvector[ i ]  );

  }
  
  for( i=0; i<4; i++ ){

    fprintf( energies_fp, DOUBLE_FORMAT"  ", potential_energy_system_p->rvector[ i ]  );

  }

  if( N_chain ){


    for( i=0; i<(N_chain +1); i++ ){

      fprintf( energies_fp, DOUBLE_FORMAT"  ", kinetic_energy_thermostat_p->rvector[ i ]  );

    }
  
    for( i=0; i<(N_chain +1); i++ ){

      fprintf( energies_fp, DOUBLE_FORMAT"  ", potential_energy_thermostat_p->rvector[ i ]  );

    }

    fprintf( energies_fp, DOUBLE_FORMAT"  ", total_energy_system );

    fprintf( energies_fp, DOUBLE_FORMAT"\n", pseudo_energy_system);

  }
  else{ 
  
    fprintf( energies_fp, DOUBLE_FORMAT"\n", total_energy_system);

  }

  fflush( energies_fp );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_mu_traces( const constants constants, const state state, const config config ){

  /* constants */
  int             N_coor;
  /* state */
  double          time;
  const matrix*        mu00_p;
  const matrix_array*  mu01_p;
  const matrix_array*  mu10_p;
  const matrix_array*  mu02_p;
  const matrix_array*  mu11_p;
  const matrix_array*  mu20_p;
  /* dummies */
  int             i, dim;
  double          dummy;
  int             info=0;


  N_coor  = constants.N_coor;

  time    =  state.time;
  mu00_p  = &(config.electrons.mu00);
  mu01_p  = &(config.electrons.mu01);
  mu10_p  = &(config.electrons.mu10);
  mu02_p  = &(config.electrons.mu02);
  mu11_p  = &(config.electrons.mu11);
  mu20_p  = &(config.electrons.mu20);


  fprintf( mu_traces_fp, "%le  ", time );

  //
  dummy =  REAL( MATRIX_TRACE( *mu00_p ) );

  fprintf( mu_traces_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;
 
  for( i=0; i<N_coor; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu01_p->array[ i ] ) );

  }

  fprintf( mu_traces_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<N_coor; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu10_p->array[ i ] ) );

  }

  fprintf( mu_traces_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  dim = N_coor *N_coor;

  for( i=0; i<dim; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu02_p->array[ i ] ) );

  }

  fprintf( mu_traces_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<dim; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu11_p->array[ i ] ) );

  }

  fprintf( mu_traces_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<dim; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu20_p->array[ i ] ) );

  }

  fprintf( mu_traces_fp, DOUBLE_FORMAT"\n", dummy );

  fflush( mu_traces_fp );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_mu_norms( const constants constants, const state state, const config config ){

  /* constants */
  int             N_coor;
  /* state */
  double          time;
  const matrix*        mu00_p;
  const matrix_array*  mu01_p;
  const matrix_array*  mu10_p;
  const matrix_array*  mu02_p;
  const matrix_array*  mu11_p;
  const matrix_array*  mu20_p;
  /* dummies */
  int             i, dim;
  double          dummy;
  int             info=0;


  N_coor  = constants.N_coor;

  time    =  state.time;
  mu00_p  = &(config.electrons.mu00);
  mu01_p  = &(config.electrons.mu01);
  mu10_p  = &(config.electrons.mu10);
  mu02_p  = &(config.electrons.mu02);
  mu11_p  = &(config.electrons.mu11);
  mu20_p  = &(config.electrons.mu20);


  fprintf( mu_norms_fp, "%le  ", time );
  
  //
  dummy =  MATRIX_NORM( *mu00_p );

  fprintf( mu_norms_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;
 
  for( i=0; i<N_coor; i++ ){

    dummy +=  MATRIX_NORM( mu01_p->array[ i ] );

  }

  fprintf( mu_norms_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<N_coor; i++ ){

    dummy +=  MATRIX_NORM( mu10_p->array[ i ] );

  }

  fprintf( mu_norms_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  dim = N_coor *N_coor;

  for( i=0; i<dim; i++ ){

    dummy +=  MATRIX_NORM( mu02_p->array[ i ] );

  }

  fprintf( mu_norms_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<dim; i++ ){

    dummy +=  MATRIX_NORM( mu11_p->array[ i ] );

  }

  fprintf( mu_norms_fp, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<dim; i++ ){

    dummy +=  MATRIX_NORM( mu20_p->array[ i ] );

  }

  fprintf( mu_norms_fp, DOUBLE_FORMAT"\n", dummy );

  fflush( mu_norms_fp );


  return info;

}

//------------------------------------------

#ifdef __DEBUG_PLUS__ 

//BUGFIX: this was originally DEF
int print_rho_traces( const constants constants, const state state, const config config ){

  /* constants */
  int            max_rho_index;
  /* state */
  double         time;
  const matrix_array* rho_p;
  int            i;
  double         dummy;
  int            info=0;


  max_rho_index    =  constants.max_rho_index;

  time             =  state.time;
  rho_p            = &(config.electrons.rho);


  fprintf( rho_traces_fp, "%le  ", time );

  for( i=0; i<(max_rho_index-1); i++){

    dummy =  REAL( MATRIX_TRACE( rho_p->array[ i ] ) );

    fprintf( rho_traces_fp, DOUBLE_FORMAT"  ", dummy );

    dummy =  IMAG( MATRIX_TRACE( rho_p->array[ i ] ) );

    fprintf( rho_traces_fp, DOUBLE_FORMAT"  ", dummy );

  }

  dummy =  REAL( MATRIX_TRACE( rho_p->array[ i ] ) );

  fprintf( rho_traces_fp, DOUBLE_FORMAT"  ", dummy );

  dummy =  IMAG( MATRIX_TRACE( rho_p->array[ i ] ) );

  fprintf( rho_traces_fp, DOUBLE_FORMAT"\n", dummy );


  fflush( rho_traces_fp );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_rho_norms( const constants constants, const state state, const config config ){


/* constants */
  int            max_rho_index;
  /* state */
  double         time;
  const matrix_array* rho_p;
  int            i;
  double         dummy;
  int            info=0;


  max_rho_index    =  constants.max_rho_index;

  time             =  state.time;
  rho_p            = &(config.electrons.rho);


  fprintf( rho_norms_fp, "%le  ", time );

  for( i=0; i<(max_rho_index-1); i++){

    dummy = MATRIX_NORM( rho_p->array[ i ] );

    fprintf( rho_norms_fp, DOUBLE_FORMAT"  ", dummy );

  }

  dummy = MATRIX_NORM( rho_p->array[ i ] );

  fprintf( rho_norms_fp, DOUBLE_FORMAT"\n", dummy );


  fflush( rho_norms_fp );


  return info;

}

#endif /* __DEBUG_PLUS__ */

//------------------------------------------

int print_adiabatic_populations( const constants constants, const state state ){

  /* constants */
  int       N_levels_many;
  /* state */
  double    time;
  const rvector*    adiabatic_populations_p;
  /* dummies */
  int       i;
  int       info=0;


  N_levels_many                     =  constants.N_levels_many;

  time                              =  state.time;
  adiabatic_populations_p           = &(state.adiabatic_populations);


  fprintf( adiabatic_populations_fp, "%le  ", time );


  for( i=0; i<N_levels_many-1; i++ ){

    fprintf( adiabatic_populations_fp, DOUBLE_FORMAT"  ", adiabatic_populations_p->rvector[ i ] );

  }
  fprintf( adiabatic_populations_fp, DOUBLE_FORMAT"\n", adiabatic_populations_p->rvector[ i ] );

  fflush( adiabatic_populations_fp );


  return info;

}

//------------------------------------------
//
int print_adiabatic_projections( const constants constants, const state state ){

  /* constants */
  int       N_levels_many;
  /* state */
  double    time;
  const rvector*   electronic_density_eigenvalues_p;
  const matrix*    adiabatic_projection_p;
  /* dummies */
  int       i, j;
  int       index;
  double    dummy;
  int       info=0;


  N_levels_many                     =  constants.N_levels_many;

  time                              = state.time;
  electronic_density_eigenvalues_p  = &(state.electronic_density_eigenvalues);
  adiabatic_projection_p            = &(state.adiabatic_projection);


  fprintf( adiabatic_projections_fp, "# time = %le\n", time );

  //
  fprintf( adiabatic_projections_fp, "# eigenvalues:  " );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( adiabatic_projections_fp, DOUBLE_FORMAT"  ", electronic_density_eigenvalues_p->rvector[ N_levels_many -1 -i ] );

  }

  fprintf( adiabatic_projections_fp, "\n" );

  //
  fprintf( adiabatic_projections_fp, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( adiabatic_projections_fp, "#%d  ", i+1 );

    for( j=0; j<N_levels_many; j++ ){

      index = ELECTRON_MANY_INDEX( j, N_levels_many -1 -i );

      dummy = CMPLX_NORM( adiabatic_projection_p->matrix[ index ] );

      fprintf( adiabatic_projections_fp, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( adiabatic_projections_fp, "\n" );

  } /* i loop */

  //
  fprintf( adiabatic_projections_fp, "\n\n" );

  fflush( adiabatic_projections_fp );


  return info;

}

//------------------------------------------

int print_projections( const state state ){

  /* state */
  double    time;
  double    rho_norm;
  double    rho_dot_norm;
  double    rho_projection;
  /* dummies */
  int       info=0;


  time             = state.time;
  rho_norm         = state.rho_norm;
  rho_dot_norm     = state.rho_dot_norm;
  rho_projection   = state.rho_projection;


  fprintf( projections_fp, "%le  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"\n", time, rho_norm, rho_dot_norm, rho_projection );

  fflush( projections_fp );


  return info;

}

//------------------------------------------

int print_adiabatic_PES_many( const constants constants, const state state ){

  /* constants */
  int       N_levels_many;
  /* state */
  double    time;
  const rvector*  adiabatic_PES_many_p;
  /* dummies */
  int       i;
  int       info=0;


  N_levels_many        = constants.N_levels_many;

  time                 = state.time;
  adiabatic_PES_many_p = &(state.adiabatic_PES_many);


  fprintf( adiabatic_PES_many_fp, "%le  ", time );

  for( i=0; i<N_levels_many-1; i++ ){

    fprintf( adiabatic_PES_many_fp, DOUBLE_FORMAT"  ", adiabatic_PES_many_p->rvector[ i ] );

  }

  fprintf( adiabatic_PES_many_fp, DOUBLE_FORMAT"\n", adiabatic_PES_many_p->rvector[ i ] );

  fflush( adiabatic_PES_many_fp );


  return info;

}

//------------------------------------------

int print_adiabatic_PES_single( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const rvector* adiabatic_PES_single_p;
  /* dummies */
  int       i;
  int       info=0;


  N_levels_single        =  constants.N_levels_single;

  time                   =  state.time;
  adiabatic_PES_single_p = &(state.adiabatic_PES_single);


  fprintf( adiabatic_PES_single_fp, "%le  ", time );

  for( i=0; i<N_levels_single-1; i++ ){

    fprintf( adiabatic_PES_single_fp, DOUBLE_FORMAT"  ", adiabatic_PES_single_p->rvector[ i ] );

  }

  fprintf( adiabatic_PES_single_fp, DOUBLE_FORMAT"\n", adiabatic_PES_single_p->rvector[ i ] );

  fflush( adiabatic_PES_single_fp );


  return info;

}

//------------------------------------------

int print_single_level_populations_Ehrenfest( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const rvector* single_level_populations_Ehrenfest_p;
  /* dummies */
  int       i;
  int       info=0;


  N_levels_single                      =  constants.N_levels_single;

  time                                 =  state.time;
  single_level_populations_Ehrenfest_p = &(state.single_level_populations_Ehrenfest);


  fprintf( single_level_populations_Ehrenfest_fp, "%le  ", time );

  for( i=0; i<N_levels_single-1; i++ ){

    fprintf( single_level_populations_Ehrenfest_fp, DOUBLE_FORMAT"  ", single_level_populations_Ehrenfest_p->rvector[ i ] );

  }

  fprintf( single_level_populations_Ehrenfest_fp, DOUBLE_FORMAT"\n", single_level_populations_Ehrenfest_p->rvector[ i ] );

  fflush( single_level_populations_Ehrenfest_fp );


  return info;

}

//------------------------------------------

int print_single_level_populations_adiabatic( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const rvector* single_level_populations_adiabatic_p;
  /* dummies */
  int       i;
  int       info=0;


  N_levels_single                      =  constants.N_levels_single;

  time                                 =  state.time;
  single_level_populations_adiabatic_p = &(state.single_level_populations_adiabatic);


  fprintf( single_level_populations_adiabatic_fp, "%le  ", time );

  for( i=0; i<N_levels_single-1; i++ ){

    fprintf( single_level_populations_adiabatic_fp, DOUBLE_FORMAT"  ", single_level_populations_adiabatic_p->rvector[ i ] );

  }

  fprintf( single_level_populations_adiabatic_fp, DOUBLE_FORMAT"\n", single_level_populations_adiabatic_p->rvector[ i ] );

  fflush( single_level_populations_adiabatic_fp );


  return info;

}

//------------------------------------------

int print_nonadiabatic_couplings( const constants constants, const state state ){

  /* constants */
  int       N_coor;
  /* state */
  double    time;
  const rvector* nonadiabatic_coupling_p;
  /* dummies */
  int       i;
  int       info=0;


  N_coor                  =  constants.N_coor;

  time                    =  state.time;
  nonadiabatic_coupling_p = &(state.nonadiabatic_coupling);


  fprintf( nonadiabatic_couplings_fp, "%le  ", time );

  for( i=0; i<N_coor-1; i++ ){

    fprintf( nonadiabatic_couplings_fp, DOUBLE_FORMAT"  ", nonadiabatic_coupling_p->rvector[ i ] );

  }

  fprintf( nonadiabatic_couplings_fp, DOUBLE_FORMAT"\n", nonadiabatic_coupling_p->rvector[ i ] );

  fflush( nonadiabatic_couplings_fp );


  return info;

}

//------------------------------------------

int print_nonadiabatic_rates( const constants constants, const state state ){

  /* constants */
  int       N_coor;
  /* state */
  double    time;
  const rvector* nonadiabatic_rate_p;
  /* dummies */
  int       i;
  int       info=0;


  N_coor                  =  constants.N_coor;

  time                    =  state.time;
  nonadiabatic_rate_p     = &(state.nonadiabatic_rate);


  fprintf( nonadiabatic_rates_fp, "%le  ", time );

  for( i=0; i<N_coor-1; i++ ){

    fprintf( nonadiabatic_rates_fp, DOUBLE_FORMAT"  ", nonadiabatic_rate_p->rvector[ i ] );

  }

  fprintf( nonadiabatic_rates_fp, DOUBLE_FORMAT"\n", nonadiabatic_rate_p->rvector[ i ] );

  fflush( nonadiabatic_rates_fp );


  return info;

}

//------------------------------------------

int print_one_body_electronic_density_matrix( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const matrix*  one_body_electronic_density_matrix_p;
  /* dummies */
  int       i, j;
  int       index;
  double    dummy;
  int       info=0;


  N_levels_single                      =  constants.N_levels_single;

  time                                 =  state.time;
  one_body_electronic_density_matrix_p = &(state.one_body_electronic_density_matrix);


  fprintf( one_body_electronic_density_matrix_fp, "# time = %le\n", time );

  for( i=0; i<N_levels_single; i++ ){

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( i, j );

      fprintf( one_body_electronic_density_matrix_fp, DOUBLE_FORMAT"  ", CMPLX_NORM( one_body_electronic_density_matrix_p->matrix[ index ]) );

    } /* j loop */

    fprintf( one_body_electronic_density_matrix_fp, "\n" );

  } /* i loop */

  //
  dummy = REAL( MATRIX_TRACE( *one_body_electronic_density_matrix_p ) );

  fprintf( one_body_electronic_density_matrix_fp, "# trace = %le  ", dummy );

  dummy = IMAG( MATRIX_TRACE( *one_body_electronic_density_matrix_p ) );

  fprintf( one_body_electronic_density_matrix_fp, "%le\n", dummy );


  //
  dummy  = MATRIX_NORM( *one_body_electronic_density_matrix_p );

  fprintf( one_body_electronic_density_matrix_fp, "# norm = %le\n", dummy );

  //
  fprintf( one_body_electronic_density_matrix_fp, "\n\n" );

  fflush( one_body_electronic_density_matrix_fp );


  return info;

}
//------------------------------------------

int print_one_body_electronic_hole_matrix( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const matrix*  one_body_electronic_hole_matrix_p;
  /* dummies */
  int       i, j;
  int       index;
  //double    dummy;
  int       info=0;


  N_levels_single                   =  constants.N_levels_single;

  time                              =  state.time;
  one_body_electronic_hole_matrix_p = &(state.one_body_electronic_hole_matrix);


  fprintf( one_body_electronic_hole_matrix_fp, "# time = %le\n", time );

  for( i=0; i<N_levels_single; i++ ){

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( i, j );

      /* This is a tentative spectral function of the transition |GS> -->|instantaneous> */
      fprintf( one_body_electronic_hole_matrix_fp, DOUBLE_FORMAT"  ", CMPLX_NORM( one_body_electronic_hole_matrix_p->matrix[ index ] ) );

    } /* j loop */

    fprintf( one_body_electronic_hole_matrix_fp, "\n" );

  } /* i loop */

  /*
  dummy = REAL( MATRIX_TRACE( *one_body_electronic_hole_matrix_p ) );

  fprintf( one_body_electronic_hole_matrix_fp, "# trace = %le ", dummy );

  dummy = IMAG( MATRIX_TRACE( *one_body_electronic_hole_matrix_p ) );

  fprintf( one_body_electronic_hole_matrix_fp, "%le \n", dummy );

  dummy  = MATRIX_NORM( *one_body_electronic_hole_matrix_p );

  fprintf( one_body_electronic_hole_matrix_fp, "# norm = %le \n", dummy );
  */

  //
  fprintf( one_body_electronic_hole_matrix_fp, "\n\n" );

  fflush( one_body_electronic_hole_matrix_fp );


  return info;

}

//------------------------------------------

int print_one_body_electronic_particle_matrix( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const matrix*  one_body_electronic_particle_matrix_p;
  /* dummies */
  int       i, j;
  int       index;
  //double    dummy;
  int       info=0;


  N_levels_single                       =  constants.N_levels_single;

  time                                  =  state.time;
  one_body_electronic_particle_matrix_p = &(state.one_body_electronic_particle_matrix);


  fprintf( one_body_electronic_particle_matrix_fp, "# time = %le\n", time );

  for( i=0; i<N_levels_single; i++ ){

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( i, j );

      /* This is a tentative spectral function of the transition |GS> -->|instantaneous> */
      fprintf( one_body_electronic_particle_matrix_fp, DOUBLE_FORMAT"  ", CMPLX_NORM( one_body_electronic_particle_matrix_p->matrix[ index ] ) );

    } /* j loop */

    fprintf( one_body_electronic_particle_matrix_fp, "\n" );

  } /* i loop */

  /*
  dummy = REAL( MATRIX_TRACE( *one_body_electronic_particle_matrix_p ) );

  fprintf( one_body_electronic_particle_matrix_fp, "# trace = %le ", dummy );

  dummy = IMAG( MATRIX_TRACE( *one_body_electronic_particle_matrix_p ) );

  fprintf( one_body_electronic_particle_matrix_fp, "%le \n", dummy );

  dummy  = MATRIX_NORM( *one_body_electronic_particle_matrix_p );

  fprintf( one_body_electronic_particle_matrix_fp, "# norm = %le \n", dummy );
  */

  //
  fprintf( one_body_electronic_particle_matrix_fp, "\n\n" );

  fflush( one_body_electronic_particle_matrix_fp );


  return info;

}

//------------------------------------------

int print_natural_orbitals( const constants constants, state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const matrix*  one_body_electronic_density_matrix_p;
  matrix_p  dummy_matrix_single1_p;
  rvector_p dummy_rvector_single_p;
  /* dummies */
  int       i, j;
  int       index;
  double    dummy;
  int       info=0;


  N_levels_single                      =  constants.N_levels_single;

  time                                 =  state.time;
  one_body_electronic_density_matrix_p = &(state.one_body_electronic_density_matrix);
  dummy_matrix_single1_p               = &(state.dummy_matrix_single1);
  dummy_rvector_single_p               = &(state.dummy_rvector_single);


  // diagonalisation
  if( DIAGONALISATION( *one_body_electronic_density_matrix_p, *dummy_matrix_single1_p, *dummy_rvector_single_p ) ) info=1;


  fprintf( natural_orbitals_fp, "# time = %le\n", time );

  //
  fprintf( natural_orbitals_fp, "# eigenvalues:  " );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( natural_orbitals_fp, DOUBLE_FORMAT"  ", dummy_rvector_single_p->rvector[ N_levels_single -1 -i ] );

  }

  fprintf( natural_orbitals_fp, "\n" );

  //
  fprintf( natural_orbitals_fp, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( natural_orbitals_fp, "#%d   ", i+1 );

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( j, N_levels_single -1 -i );  // WARNING: note the reverse ordering

      dummy = CMPLX_NORM( dummy_matrix_single1_p->matrix[ index ] );

      fprintf( natural_orbitals_fp, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( natural_orbitals_fp, "\n" );

  } /* i loop */

  //
  fprintf( natural_orbitals_fp, "\n\n" );

  fflush( natural_orbitals_fp );


  return info;

}

//------------------------------------------

int print_hole_orbitals( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const matrix*  hole_orbitals_p;
  const rvector* hole_populations_p;
  /* dummies */
  int       i, j;
  int       index;
  double    dummy;
  int       info=0;


  N_levels_single                         =  constants.N_levels_single;

  time                                    =  state.time;
  hole_orbitals_p                         = &(state.hole_orbitals);
  hole_populations_p                      = &(state.hole_populations);


  fprintf( hole_orbitals_fp, "# time = %le\n", time );

  //
  fprintf( hole_orbitals_fp, "# eigenvalues:  " );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( hole_orbitals_fp, DOUBLE_FORMAT"  ", hole_populations_p->rvector[ N_levels_single -1 -i ] );

  }

  fprintf( hole_orbitals_fp, "\n" );

  //
  fprintf( hole_orbitals_fp, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( hole_orbitals_fp, "#%d  ", i+1 );

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( j, N_levels_single -1 -i );  // WARNING: note the reverse ordering

      dummy = CMPLX_NORM( hole_orbitals_p->matrix[ index ] );

      fprintf( hole_orbitals_fp, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( hole_orbitals_fp, "\n" );

  } /* i loop */

  //
  fprintf( hole_orbitals_fp, "\n\n" );

  fflush( hole_orbitals_fp );


  return info;

}

//------------------------------------------

int print_particle_orbitals( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const matrix*  particle_orbitals_p;
  const rvector* particle_populations_p;
  /* dummies */
  int       i, j;
  int       index;
  double    dummy;
  int       info=0;


  N_levels_single                         =  constants.N_levels_single;

  time                                    =  state.time;
  particle_orbitals_p                     = &(state.particle_orbitals);
  particle_populations_p                  = &(state.particle_populations);


  fprintf( particle_orbitals_fp, "# time = %le\n", time );

  //
  fprintf( particle_orbitals_fp, "# eigenvalues:  " );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( particle_orbitals_fp, DOUBLE_FORMAT"  ", particle_populations_p->rvector[ N_levels_single -1 -i ] );

  }

  fprintf( particle_orbitals_fp, "\n" );

  //
  fprintf( particle_orbitals_fp, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( particle_orbitals_fp, "#%d  ", i+1 );

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( j, N_levels_single -1 -i );  // WARNING: note the reverse ordering

      dummy = CMPLX_NORM( particle_orbitals_p->matrix[ index ] );

      fprintf( particle_orbitals_fp, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( particle_orbitals_fp, "\n" );

  } /* i loop */

  //
  fprintf( particle_orbitals_fp, "\n\n" );

  fflush( particle_orbitals_fp );


  return info;

}

//------------------------------------------

int print_adiabatic_states( const constants constants, const state state ){

  /* constants */
  int       N_levels_many;
  /* state */
  double    time;
  const rvector*   adiabatic_PES_many_p;
  const matrix*    adiabatic_states_many_p;
  /* dummies */
  int       i, j;
  int       index;
  double    dummy;
  int       info=0;


  N_levels_many            =  constants.N_levels_many;

  time                     =  state.time;
  adiabatic_PES_many_p     = &(state.adiabatic_PES_many);
  adiabatic_states_many_p  = &(state.adiabatic_states_many);


  fprintf( adiabatic_states_fp, "# time = %le\n", time );

  //
  fprintf( adiabatic_states_fp, "# eigenvalues:  " );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( adiabatic_states_fp, DOUBLE_FORMAT"  ", adiabatic_PES_many_p->rvector[ i ] );

  }

  fprintf( adiabatic_states_fp, "\n" );

  //
  fprintf( adiabatic_states_fp, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( adiabatic_states_fp, "#%d  ", i+1 );

    for( j=0; j<N_levels_many; j++ ){

      index = ELECTRON_MANY_INDEX( j, i );

      dummy = CMPLX_NORM( adiabatic_states_many_p->matrix[ index ] );

      fprintf( adiabatic_states_fp, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( adiabatic_states_fp, "\n" );

  } /* i loop */

  //
  fprintf( adiabatic_states_fp, "\n\n" );

  fflush( adiabatic_states_fp );


  return info;

}

//------------------------------------------

int print_electronic_density_states( const constants constants, const state state ){

  /* constants */
  int       N_levels_many;
  /* state */
  double    time;
  const rvector*   electronic_density_eigenvalues_p;
  const matrix*    electronic_density_eigenvectors_p;
  /* dummies */
  int       i, j;
  int       index;
  double    dummy;
  int       info=0;


  N_levels_many                     =  constants.N_levels_many;

  time                              =  state.time;
  electronic_density_eigenvalues_p  = &(state.electronic_density_eigenvalues);
  electronic_density_eigenvectors_p = &(state.electronic_density_eigenvectors);


  fprintf( electronic_density_states_fp, "# time = %le\n", time );

  //
  fprintf( electronic_density_states_fp, "# eigenvalues:  " );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( electronic_density_states_fp, DOUBLE_FORMAT"  ", electronic_density_eigenvalues_p->rvector[ N_levels_many -1 -i ] );

  }

  fprintf( electronic_density_states_fp, "\n" );

  //
  fprintf( electronic_density_states_fp, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( electronic_density_states_fp, "#%d  ", i+1 );

    for( j=0; j<N_levels_many; j++ ){

      index = ELECTRON_MANY_INDEX( j, N_levels_many -1 -i );

      dummy = CMPLX_NORM( electronic_density_eigenvectors_p->matrix[ index ] );

      fprintf( electronic_density_states_fp, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( electronic_density_states_fp, "\n" );

  } /* i loop */

  //
  fprintf( electronic_density_states_fp, "\n\n" );

  fflush( electronic_density_states_fp );


  return info;

}

//------------------------------------------

int print_ionic_density_states( const constants constants, const state state ){

  /* constants */
  int       sqrt_max_rho_index;
  /* state */
  double    time;
  const rvector*   ionic_density_eigenvalues_p;
  const matrix*    ionic_density_eigenvectors_p;
  /* dummies */
  int       i, j;
  int       index;
  double    dummy;
  int       info=0;


  sqrt_max_rho_index           =  constants.sqrt_max_rho_index;

  time                         =  state.time;
  ionic_density_eigenvalues_p  = &(state.ionic_density_eigenvalues);
  ionic_density_eigenvectors_p = &(state.ionic_density_eigenvectors);


  fprintf( ionic_density_states_fp, "# time = %le\n", time );

  //
  fprintf( ionic_density_states_fp, "# eigenvalues:  " );

  for( i=0; i<sqrt_max_rho_index; i++ ){

    fprintf( ionic_density_states_fp, DOUBLE_FORMAT"  ", ionic_density_eigenvalues_p->rvector[ sqrt_max_rho_index -1 -i ] );

  }

  fprintf( ionic_density_states_fp, "\n" );

  //
  fprintf( ionic_density_states_fp, "# eigenvectors squared =\n" );

  for( i=0; i<sqrt_max_rho_index; i++ ){

    fprintf( ionic_density_states_fp, "#%d  ", i+1 );

    for( j=0; j<sqrt_max_rho_index; j++ ){

      index = RHO_INDEX_AUX( j, sqrt_max_rho_index -1 -i );

      dummy = CMPLX_NORM( ionic_density_eigenvectors_p->matrix[ index ] );

      fprintf( ionic_density_states_fp, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( ionic_density_states_fp, "\n" );

  } /* i loop */

  //
  fprintf( ionic_density_states_fp, "\n\n" );

  fflush( ionic_density_states_fp );


  return info;

}

//------------------------------------------

int print_dipoles_many( const constants constants, const state state ){

  /* constants */
  int       N_levels_many;
  /* states */
  double    time;
  const rvector*   oscillator_strength_many_p;
  const rvector*   oscillator_frequency_many_p;
  const rvector*   adiabatic_populations_p;
  /* dummies */
  int       i, j;
  int       index;
  int       index_ratio;
  double    ratio;
  int       info=0;
  

  N_levels_many                =  constants.N_levels_many;

  time                         =  state.time;
  oscillator_strength_many_p   = &(state.observables.oscillator_strength_many);
  oscillator_frequency_many_p  = &(state.observables.oscillator_frequency_many);
  adiabatic_populations_p      = &(state.adiabatic_populations);


  fprintf( dipoles_many_fp, "%le", time );


  index=0;

  for( i=0; i<N_levels_many; i++ ){

    for( j=(i+1); j<N_levels_many; j++ ){
    
       // index_ratio = i; // Absorption
       index_ratio = j; // Emission

       ratio = adiabatic_populations_p->rvector[ index_ratio ];

       fprintf( dipoles_many_fp, "  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"", oscillator_frequency_many_p->rvector[ index ], ratio *oscillator_strength_many_p->rvector[ index ], oscillator_strength_many_p->rvector[ index ] );

       index++;

     } /* end loop j */  

  } /* end loop i */

  fprintf( dipoles_many_fp, "\n" );


  return info;

}

//------------------------------------------

int print_dipoles_single( const constants constants, const state state ){

  /* constants */
  int       N_levels_single;
  /* state */
  double    time;
  const rvector*   hole_populations_p;
  const rvector*   oscillator_strength_single_p;
  const rvector*   oscillator_frequency_single_p;
  /* dummies */
  double    ratio;
  int       i;
  int       info=0;


  N_levels_single               =  constants.N_levels_single;

  time                          =  state.time;
  hole_populations_p            = &(state.hole_populations);
  oscillator_strength_single_p  = &(state.observables.oscillator_strength_single);
  oscillator_frequency_single_p = &(state.observables.oscillator_frequency_single);


  fprintf( dipoles_single_fp, "%le  ", time );

  for( i=0; i<(N_levels_single-1); i++ ){
                  
    ratio = ( hole_populations_p->rvector[ N_levels_single -1 -i ] );

    fprintf( dipoles_single_fp, DOUBLE_FORMAT"  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"  ", oscillator_frequency_single_p->rvector[ N_levels_single -1 -i ], ratio *oscillator_strength_single_p->rvector[ N_levels_single -1 -i ], oscillator_strength_single_p->rvector[ N_levels_single -1 -i ] );

  }

  ratio = ( hole_populations_p->rvector[ N_levels_single -1 -i ] );

  fprintf( dipoles_single_fp, DOUBLE_FORMAT"  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"\n", oscillator_frequency_single_p->rvector[ N_levels_single -1 -i ], ratio *oscillator_strength_single_p->rvector[ N_levels_single -1 -i ], oscillator_strength_single_p->rvector[ N_levels_single -1 -i ] );


  fflush( dipoles_single_fp );


  return info;

}

//------------------------------------------

//BUGFIX: thsi was originally DEF
int print_geometry( const constants constants, const state state, const config config ){

  /* constants */
  int        N_atoms;
  int        sdim;
  /* state */
  double     time;
  const rvector*  positions_p;
  /* dummies */
  int        i_atom;
  int        info=0;


  sdim = constants.spacial_dimension;

#ifdef __DEBUG__

  if( sdim > 3 ){

    fprintf( stderr, "ERROR: sdim >3 is not admitted! [What are you trying to simulate?]\n" );
    fflush( stderr );

    info=1;

  }

#endif /* __DEBUG__ */


  N_atoms          =  constants.N_atoms;
  sdim             =  constants.spacial_dimension; 

  time             =  state.time;
  positions_p      = &(config.atoms.positions);

  fprintf( geometry_fp, "%d\n", N_atoms );
  fprintf( geometry_fp, "# time = %le\n", time );

  if( !info ){

    for( i_atom=0; i_atom<N_atoms; i_atom++ ){

      fprintf( geometry_fp, "%s  "DOUBLE_FORMAT"  ", "C", positions_p->rvector[ i_atom *sdim ] );
      
      if( sdim > 1 ){

	fprintf( geometry_fp, DOUBLE_FORMAT"  ", positions_p->rvector[ 1 +i_atom *sdim ] );

      }
      else{
      
	fprintf( geometry_fp, DOUBLE_FORMAT"  ", 0.0e0  );
      
      }

      if( sdim > 2 ){

	fprintf( geometry_fp, DOUBLE_FORMAT"\n", positions_p->rvector[ 2 +i_atom *sdim ] );
      
      }
      else{
      
	fprintf( geometry_fp, DOUBLE_FORMAT"\n", 0.0e0  );
      
      }

    } /* end for loop */

  } /* end info if*/

  fflush(  geometry_fp );


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int open_file( FILE** fpp, char* filename ){

  /* dummies */
  int info=0;


  *fpp = fopen( filename, "wt" );

  if( FILE_CHECK( *fpp, open_file ) ){

    fprintf( stderr, "ERROR occurred when opening file %s.\n", filename );

    fflush( stderr );

    info=1;

  }


  return info;

}	

//------------------------------------------

int close_file( FILE** fpp, char* filename ){

  /* dummies */
  int info=0;


  if( FILE_CHECK( *fpp, close_file ) ){

    fprintf( stderr, "ERROR occurred when closing file %s.\n", filename );

    fflush( stderr );

    info=1;

  }

  fclose( *fpp );


  return info;

}	

//------------------------------------------

int output_files_opening( const constants constants ){

  /* constants */
  int       N_levels_many;
  int       N_chain;
  /* dummies */
  char      buffer[MAX_STRING_LENGTH];
  int       info=0;


  N_levels_many = constants.N_levels_many;
  N_chain       = constants.N_chain;


  /* geometry file */
  if( constants.flag_observable_all ){

    sprintf( buffer, "%s.xyz", constants.output_label );

    OPEN_FILE( geometry_fp, buffer );	

  }


  /* positions file */
  if( constants.flag_observable_all || constants.flag_observable_positions ){

    sprintf( buffer, "positions_%s.dat", constants.output_label );

    OPEN_FILE( positions_fp, buffer );	

  }


  /* momenta file */
  if( constants.flag_observable_all || constants.flag_observable_momenta ){

    sprintf( buffer, "momenta_%s.dat", constants.output_label );

    OPEN_FILE( momenta_fp, buffer );	

  }


  /* forces file */
  if( constants.flag_observable_all || constants.flag_observable_forces ){

    sprintf( buffer, "forces_%s.dat", constants.output_label );

    OPEN_FILE( forces_fp, buffer );	

  }


  if( N_chain ){

    /* positions thermostat file */
    if( constants.flag_observable_all ){

      sprintf( buffer, "positions_thermostat_%s.dat", constants.output_label );

      OPEN_FILE( positions_thermostat_fp, buffer );	

    }


    /* momenta thermostat file */
    if( constants.flag_observable_all ){

      sprintf( buffer, "momenta_thermostat_%s.dat", constants.output_label );

      OPEN_FILE( momenta_thermostat_fp, buffer );	

    }


    /* forces thermostat file */
    if( constants.flag_observable_all ){

      sprintf( buffer, "forces_thermostat_%s.dat", constants.output_label );

      OPEN_FILE( forces_thermostat_fp, buffer );	

    }

  }


  /* populations file */
  if( constants.flag_observable_all || constants.flag_observable_populations ){

    sprintf( buffer, "populations_%s.dat", constants.output_label );

    OPEN_FILE( populations_fp, buffer );	

  }


  if( constants.flag_observable_all || constants.flag_observable_mus ){

    /* mu traces file */
    sprintf( buffer, "mu_traces_%s.dat", constants.output_label );

    OPEN_FILE( mu_traces_fp, buffer );	

    /* mu norms file */
    sprintf( buffer, "mu_norms_%s.dat", constants.output_label );

    OPEN_FILE( mu_norms_fp, buffer );	

  }


#ifdef __DEBUG_PLUS__

  /* rho traces file */
  if( constants.flag_observable_all ){

    sprintf( buffer, "rho_traces_%s.dat", constants.output_label );

    OPEN_FILE( rho_traces_fp, buffer );	

  }


  /* rho norms file */
  if( constants.flag_observable_all ){

    sprintf( buffer, "rho_norms_%s.dat", constants.output_label );

    OPEN_FILE( rho_norms_fp, buffer );	

  }

#endif /* __DEBUG_PLUS__ */


  /* energies file */
  if( constants.flag_observable_all || constants.flag_observable_energies ){

    sprintf( buffer, "energies_%s.dat", constants.output_label );

    OPEN_FILE( energies_fp, buffer );	

  }
    

  /* adiabatic_populations file */
  if( constants.flag_observable_all || constants.flag_observable_adiabatic_populations ){

    sprintf( buffer, "adiabatic_populations_%s.dat", constants.output_label );

    OPEN_FILE( adiabatic_populations_fp, buffer );	

  }


  if( constants.flag_observable_all || constants.flag_observable_projections ){

    /* adiabatic_projections file */
    sprintf( buffer, "adiabatic_projections_%s.dat", constants.output_label );

    OPEN_FILE( adiabatic_projections_fp, buffer );	

    /* projections file */
    sprintf( buffer, "projections_%s.dat", constants.output_label );

    OPEN_FILE( projections_fp, buffer );	

  }


  /* adiabatic_PES_many file */
  if( constants.flag_observable_all || constants.flag_observable_adiabatic_pes_many ){

    sprintf( buffer, "adiabatic_PES_many_%s.dat", constants.output_label );

    OPEN_FILE( adiabatic_PES_many_fp, buffer );	

  }


  /* adiabatic_PES_single file */
  if( constants.flag_observable_all || constants.flag_observable_adiabatic_pes_single ){

    sprintf( buffer, "adiabatic_PES_single_%s.dat", constants.output_label );

    OPEN_FILE( adiabatic_PES_single_fp, buffer );	

  }


  if( constants.flag_observable_all || constants.flag_observable_single_level_populations ){

    /* single_level_populations_Ehrenfest file */
    sprintf( buffer, "single_level_populations_Ehrenfest_%s.dat", constants.output_label );

    OPEN_FILE( single_level_populations_Ehrenfest_fp, buffer );	
  
    /* single_level_populations_adiabatic file */
    sprintf( buffer, "single_level_populations_adiabatic_%s.dat", constants.output_label );

    OPEN_FILE( single_level_populations_adiabatic_fp, buffer );	

  }


  /* nonadiabatic_couplings file */
  if( constants.flag_observable_all || constants.flag_observable_nonadiabatic_couplings ){

    sprintf( buffer, "nonadiabatic_couplings_%s.dat", constants.output_label );

    OPEN_FILE( nonadiabatic_couplings_fp, buffer );	

  }


  /* nonadiabatic_rates file */
  /*
  if( constants.flag_observable_all ){

      sprintf( buffer, "nonadiabatic_rates_%s.dat", constants.output_label );

      OPEN_FILE( nonadiabatic_rates_fp, buffer );	

  }
  */


  if( constants.flag_observable_all || constants.flag_observable_density_matrix ){

    /* one_body_electronic_density_matrix file */
    sprintf( buffer, "one_body_electronic_density_matrix_%s.dat", constants.output_label );

    OPEN_FILE( one_body_electronic_density_matrix_fp, buffer );	

    /* natural_orbitals file */
    sprintf( buffer, "natural_orbitals_%s.dat", constants.output_label );

    OPEN_FILE( natural_orbitals_fp, buffer );	

  }


  if( constants.flag_observable_all || constants.flag_observable_transition_matrices ){

    /* one_body_electronic_hole_matrix file */
    sprintf( buffer, "one_body_electronic_hole_matrix_%s.dat", constants.output_label );

    OPEN_FILE( one_body_electronic_hole_matrix_fp, buffer );	

    /* one_body_electronic_particle_matrix file */
    sprintf( buffer, "one_body_electronic_particle_matrix_%s.dat", constants.output_label );

    OPEN_FILE( one_body_electronic_particle_matrix_fp, buffer );	

    /* hole_orbitals file */
    sprintf( buffer, "hole_orbitals_%s.dat", constants.output_label );

    OPEN_FILE( hole_orbitals_fp, buffer );	

    /* particle_orbitals file */
    sprintf( buffer, "particle_orbitals_%s.dat", constants.output_label );

    OPEN_FILE( particle_orbitals_fp, buffer );	

  }


  /* adiabatic_states file */
  if( constants.flag_observable_all || constants.flag_observable_adiabatic_states ){

    sprintf( buffer, "adiabatic_states_%s.dat", constants.output_label );

    OPEN_FILE( adiabatic_states_fp, buffer );	

  }


  /* electronic_density_states file */
  if( constants.flag_observable_all || constants.flag_observable_electronic_density_states ){

    sprintf( buffer, "electronic_density_states_%s.dat", constants.output_label );

    OPEN_FILE( electronic_density_states_fp, buffer );	

  }


  /* ionic_density_states file */
  if( constants.flag_observable_all ){

    sprintf( buffer, "ionic_density_states_%s.dat", constants.output_label );

    OPEN_FILE( ionic_density_states_fp, buffer );	

  }


  /* dipoles many file */
  if( ( N_levels_many > 1 && constants.flag_observable_all ) || constants.flag_observable_dipoles_many ){

    sprintf( buffer, "dipoles_many_%s.dat", constants.output_label );

    OPEN_FILE( dipoles_many_fp, buffer );	

  }


  /* dipoles single file */
  if( constants.flag_observable_all || constants.flag_observable_dipoles_single ){

    sprintf( buffer, "dipoles_single_%s.dat", constants.output_label );

    OPEN_FILE( dipoles_single_fp, buffer );	

  }


  return info;

}

//------------------------------------------

int output_files_closing( const constants constants ){

  /* constants */
  int       N_levels_many;
  int       N_chain;
  /* dummies */
  char      buffer[MAX_STRING_LENGTH];
  int       info=0;

                          
  N_levels_many = constants.N_levels_many;
  N_chain       = constants.N_chain;

  /* geometry file */
  if( constants.flag_observable_all ){

    sprintf( buffer, "%s.geometry", constants.output_label );

    CLOSE_FILE( geometry_fp, buffer );	

  }


  /* positions file */
  if( constants.flag_observable_all || constants.flag_observable_positions ){

    sprintf( buffer, "positions_%s.dat", constants.output_label );

    CLOSE_FILE( positions_fp, buffer );	

  }


  /* momenta file */
  if( constants.flag_observable_all || constants.flag_observable_momenta ){

    sprintf( buffer, "momenta_%s.dat", constants.output_label );

    CLOSE_FILE( momenta_fp, buffer );	

  }


  /* forces file */
  if( constants.flag_observable_all || constants.flag_observable_forces ){

    sprintf( buffer, "forces_%s.dat", constants.output_label );

    CLOSE_FILE( forces_fp, buffer );	

  }


  if( N_chain ){

    /* positions thermostat file */
    if( constants.flag_observable_all ){

      sprintf( buffer, "positions_thermostat_%s.dat", constants.output_label );

      CLOSE_FILE( positions_thermostat_fp, buffer );	

    }


    /* momenta thermostat file */
    if( constants.flag_observable_all ){

      sprintf( buffer, "momenta_thermostat_%s.dat", constants.output_label );

      CLOSE_FILE( momenta_thermostat_fp, buffer );	

    }


    /* forces thermostat file */
    if( constants.flag_observable_all ){

      sprintf( buffer, "forces_thermostat_%s.dat", constants.output_label );

      CLOSE_FILE( forces_thermostat_fp, buffer );	

    }

  }


  /* populations file */
  if( constants.flag_observable_all || constants.flag_observable_populations ){

    sprintf( buffer, "populations_%s.dat", constants.output_label );

    CLOSE_FILE( populations_fp, buffer );	

  }


  if( constants.flag_observable_all || constants.flag_observable_mus ){

    /* mu traces file */
    sprintf( buffer, "mu_traces_%s.dat", constants.output_label );

    CLOSE_FILE( mu_traces_fp, buffer );	

    /* mu norms file */
    sprintf( buffer, "mu_norms_%s.dat", constants.output_label );

    CLOSE_FILE( mu_norms_fp, buffer );	

  }


#ifdef __DEBUG_PLUS__

  /* rho traces file */
  if( constants.flag_observable_all ){

    sprintf( buffer, "rho_traces_%s.dat", constants.output_label );

    CLOSE_FILE( rho_traces_fp, buffer );	

  }


  /* rho norms file */
  if( constants.flag_observable_all ){

    sprintf( buffer, "rho_norms_%s.dat", constants.output_label );

    CLOSE_FILE( rho_norms_fp, buffer );	

  }

#endif /* __DEBUG_PLUS__ */


  /* energies file */
  if( constants.flag_observable_all || constants.flag_observable_energies ){

    sprintf( buffer, "energies_%s.dat", constants.output_label );

    CLOSE_FILE( energies_fp, buffer );	

  }
    

  /* adiabatic_populations file */
  if( constants.flag_observable_all || constants.flag_observable_adiabatic_populations ){

    sprintf( buffer, "adiabatic_populations_%s.dat", constants.output_label );

    CLOSE_FILE( adiabatic_populations_fp, buffer );	

  }


  if( constants.flag_observable_all || constants.flag_observable_projections ){

    /* adiabatic_projections file */
    sprintf( buffer, "adiabatic_projections_%s.dat", constants.output_label );

    CLOSE_FILE( adiabatic_projections_fp, buffer );	

    /* projections file */
    sprintf( buffer, "projections_%s.dat", constants.output_label );

    CLOSE_FILE( projections_fp, buffer );	

  }


  /* adiabatic_PES_many file */
  if( constants.flag_observable_all || constants.flag_observable_adiabatic_pes_many ){

    sprintf( buffer, "adiabatic_PES_many_%s.dat", constants.output_label );

    CLOSE_FILE( adiabatic_PES_many_fp, buffer );	

  }


  /* adiabatic_PES_single file */
  if( constants.flag_observable_all || constants.flag_observable_adiabatic_pes_single ){

    sprintf( buffer, "adiabatic_PES_single_%s.dat", constants.output_label );

    CLOSE_FILE( adiabatic_PES_single_fp, buffer );	

  }


  if( constants.flag_observable_all || constants.flag_observable_single_level_populations ){

    /* single_level_populations_Ehrenfest file */
    sprintf( buffer, "single_level_populations_Ehrenfest_%s.dat", constants.output_label );

    CLOSE_FILE( single_level_populations_Ehrenfest_fp, buffer );	
  
    /* single_level_populations_adiabatic file */
    sprintf( buffer, "single_level_populations_adiabatic_%s.dat", constants.output_label );

    CLOSE_FILE( single_level_populations_adiabatic_fp, buffer );	

  }


  /* nonadiabatic_couplings file */
  if( constants.flag_observable_all || constants.flag_observable_nonadiabatic_couplings ){

    sprintf( buffer, "nonadiabatic_couplings_%s.dat", constants.output_label );

    CLOSE_FILE( nonadiabatic_couplings_fp, buffer );	

  }


  /* nonadiabatic_rates file */
  /*
  if( constants.flag_observable_all ){

      sprintf( buffer, "nonadiabatic_rates_%s.dat", constants.output_label );

      CLOSE_FILE( nonadiabatic_rates_fp, buffer );	

  }
  */


  if( constants.flag_observable_all || constants.flag_observable_density_matrix ){

    /* one_body_electronic_density_matrix file */
    sprintf( buffer, "one_body_electronic_density_matrix_%s.dat", constants.output_label );

    CLOSE_FILE( one_body_electronic_density_matrix_fp, buffer );	

    /* natural_orbitals file */
    sprintf( buffer, "natural_orbitals_%s.dat", constants.output_label );

    CLOSE_FILE( natural_orbitals_fp, buffer );	

  }


  if( constants.flag_observable_all || constants.flag_observable_transition_matrices ){

    /* one_body_electronic_hole_matrix file */
    sprintf( buffer, "one_body_electronic_hole_matrix_%s.dat", constants.output_label );

    CLOSE_FILE( one_body_electronic_hole_matrix_fp, buffer );	

    /* one_body_electronic_particle_matrix file */
    sprintf( buffer, "one_body_electronic_particle_matrix_%s.dat", constants.output_label );

    CLOSE_FILE( one_body_electronic_particle_matrix_fp, buffer );	

    /* hole_orbitals file */
    sprintf( buffer, "hole_orbitals_%s.dat", constants.output_label );

    CLOSE_FILE( hole_orbitals_fp, buffer );	

    /* particle_orbitals file */
    sprintf( buffer, "particle_orbitals_%s.dat", constants.output_label );

    CLOSE_FILE( particle_orbitals_fp, buffer );	

  }


  /* adiabatic_states file */
  if( constants.flag_observable_all || constants.flag_observable_adiabatic_states ){

    sprintf( buffer, "adiabatic_states_%s.dat", constants.output_label );

    CLOSE_FILE( adiabatic_states_fp, buffer );	

  }


  /* electronic_density_states file */
  if( constants.flag_observable_all || constants.flag_observable_electronic_density_states ){

    sprintf( buffer, "electronic_density_states_%s.dat", constants.output_label );

    CLOSE_FILE( electronic_density_states_fp, buffer );	

  }


  /* ionic_density_states file */
  if( constants.flag_observable_all ){

    sprintf( buffer, "ionic_density_states_%s.dat", constants.output_label );

    CLOSE_FILE( ionic_density_states_fp, buffer );	

  }


  /* dipoles many file */
  if( ( N_levels_many > 1 && constants.flag_observable_all ) || constants.flag_observable_dipoles_many ){

    sprintf( buffer, "dipoles_many_%s.dat", constants.output_label );

    CLOSE_FILE( dipoles_many_fp, buffer );	

  }


  /* dipoles single file */
  if( constants.flag_observable_all || constants.flag_observable_dipoles_single ){

    sprintf( buffer, "dipoles_single_%s.dat", constants.output_label );

    CLOSE_FILE( dipoles_single_fp, buffer );	

  }


  return info;

}

//------------------------------------------
