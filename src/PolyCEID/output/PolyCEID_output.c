
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
#include "PolyCEID_output.h"


/* private varaiables */

FILE*      xyz_file_p;
FILE*      positions_file_p;
FILE*      momenta_file_p;
FILE*      forces_file_p;
FILE*      positions_thermostat_file_p;
FILE*      momenta_thermostat_file_p;
FILE*      forces_thermostat_file_p;
FILE*      populations_file_p;
FILE*      energies_file_p;
FILE*      mu_trace_file_p;
FILE*      mu_norm_file_p;
#ifdef __DEBUG_PLUS__
FILE*      rho_trace_file_p;
FILE*      rho_norm_file_p;
#endif /* __DEBUG_PLUS__ */
FILE*      adiabatic_populations_file_p;
FILE*      adiabatic_projection_file_p;
FILE*      projection_file_p;
FILE*      adiabatic_PES_many_file_p;
FILE*      adiabatic_PES_single_file_p;
FILE*      single_level_populations_Ehrenfest_file_p;
FILE*      single_level_populations_adiabatic_file_p;
FILE*      nonadiabatic_coupling_file_p;
FILE*      nonadiabatic_rate_file_p;
FILE*      one_body_electronic_density_matrix_file_p;
FILE*      one_body_electronic_hole_matrix_file_p;
FILE*      one_body_electronic_particle_matrix_file_p;
FILE*      natural_orbitals_file_p;
FILE*      hole_orbitals_file_p;
FILE*      particle_orbitals_file_p;
FILE*      adiabatic_states_file_p;
FILE*      electronic_density_states_file_p;
FILE*      ionic_density_states_file_p;
FILE*      dipole_many_file_p;
FILE*      dipole_single_file_p;


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
	     

  if( PRINT_POSITIONS( constants, state, config ) )                   info=1;

  if( PRINT_MOMENTA( constants, state, config ) )                     info=1;

  if( PRINT_FORCES( constants, state, config ) )                      info=1;

  if( N_chain ){

    if( PRINT_POSITIONS_THERMOSTAT( constants, state, config ) )      info=1;

    if( PRINT_MOMENTA_THERMOSTAT( constants, state, config ) )        info=1;

    if( PRINT_FORCES_THERMOSTAT( constants, state, config ) )         info=1;

  }        

  if( PRINT_POPULATIONS( constants, state, config ) )                 info=1;

  if( PRINT_MU_TRACE( constants, state, config ) )                    info=1;

  if( PRINT_MU_NORM( constants, state, config ) )                     info=1;

#ifdef __DEBUG_PLUS__

  if( PRINT_RHO_TRACE( constants, state, config ) )                   info=1;

  if( PRINT_RHO_NORM( constants, state, config ) )                    info=1;

#endif /* __DEBUG_PLUS__ */

  if( PRINT_ENERGIES( constants, state ) )                            info=1;

  if( PRINT_ADIABATIC_POPULATIONS( constants, state ) )               info=1;

  if( PRINT_ADIABATIC_PROJECTION( constants, state ) )                info=1;

  if( PRINT_PROJECTION( state ) )                                     info=1;

  if( PRINT_ADIABATIC_PES_MANY( constants, state ) )                  info=1;

  if( PRINT_ADIABATIC_PES_SINGLE( constants, state ) )                info=1;

  if( PRINT_SINGLE_LEVEL_POPULATIONS_EHRENFEST( constants, state ) )  info=1;

  if( PRINT_SINGLE_LEVEL_POPULATIONS_ADIABATIC( constants, state ) )  info=1;

  if( PRINT_NONADIABATIC_COUPLING( constants, state ) )               info=1; 

  // if( PRINT_NONADIABATIC_RATE( constants, state ) )                   info=1; 

  if( PRINT_ONE_BODY_ELECTRONIC_DENSITY_MATRIX( constants, state ) )  info=1;

  if( PRINT_ONE_BODY_ELECTRONIC_HOLE_MATRIX( constants, state ) )     info=1;

  if( PRINT_ONE_BODY_ELECTRONIC_PARTICLE_MATRIX( constants, state ) ) info=1;

  if( PRINT_NATURAL_ORBITALS( constants, state ) )                    info=1;

  if( PRINT_HOLE_ORBITALS( constants, state ) )                       info=1;

  if( PRINT_PARTICLE_ORBITALS( constants, state ) )                   info=1;

  if( PRINT_ADIABATIC_STATES( constants, state ) )                    info=1;

  if( PRINT_ELECTRONIC_DENSITY_STATES( constants, state ) )           info=1;

  if( PRINT_IONIC_DENSITY_STATES( constants, state ) )                info=1;

  if( N_levels_many > 1 ){

    if( PRINT_DIPOLE_MANY( constants, state ) )                       info=1;

  }

  if( PRINT_DIPOLE_SINGLE( constants, state ) )                       info=1;

  if( PRINT_FRAME( constants, state, config ) )                       info=1;


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


  fprintf( positions_file_p, "%le  ", time );

  for( i=0; i<(N_coor-1); i++ ){

    fprintf( positions_file_p, DOUBLE_FORMAT"  ", positions_p->rvector[ i ] );

  }
  fprintf( positions_file_p, DOUBLE_FORMAT"\n", positions_p->rvector[ i ] );

  fflush( positions_file_p );


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


  fprintf( momenta_file_p, "%le  ", time );

  for( i=0; i<(N_coor-1); i++ ){

    fprintf( momenta_file_p, DOUBLE_FORMAT"  ", momenta_p->rvector[ i ] );

  }
  fprintf( momenta_file_p, DOUBLE_FORMAT"\n", momenta_p->rvector[ i ] );

  fflush( momenta_file_p );


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


  fprintf( forces_file_p, "%le   ", time );

  for( i=0; i<(N_coor-1); i++ ){

    fprintf( forces_file_p, DOUBLE_FORMAT"  ", forces_p->rvector[ i ] );

  }
  fprintf( forces_file_p, DOUBLE_FORMAT"\n", forces_p->rvector[ i ] );

  fflush( forces_file_p );


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


  fprintf( positions_thermostat_file_p, "%le  ", time );

  for( i=0; i<(N_chain-1); i++ ){

    fprintf( positions_thermostat_file_p, DOUBLE_FORMAT"  ", positions_thermostat_p->rvector[ i ] );

  }
  fprintf( positions_thermostat_file_p, DOUBLE_FORMAT"\n", positions_thermostat_p->rvector[ i ] );

  fflush( positions_thermostat_file_p );


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


  fprintf( momenta_thermostat_file_p, "%le  ", time );

  for( i=0; i<(N_chain-1); i++ ){

    fprintf( momenta_thermostat_file_p, DOUBLE_FORMAT"  ", momenta_thermostat_p->rvector[ i ] );

  }
  fprintf( momenta_thermostat_file_p, DOUBLE_FORMAT"\n", momenta_thermostat_p->rvector[ i ] );

  fflush( momenta_thermostat_file_p );


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


  fprintf( forces_thermostat_file_p, "%le  ", time );

  for( i=0; i<(N_chain-1); i++ ){

    fprintf( forces_thermostat_file_p, DOUBLE_FORMAT"  ", forces_thermostat_p->rvector[ i ] );

  }
  fprintf( forces_thermostat_file_p, DOUBLE_FORMAT"\n", forces_thermostat_p->rvector[ i ] );

  fflush( forces_thermostat_file_p );


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


  fprintf( populations_file_p, "%le  ", time); 

  for( i=0; i<(N_levels_many-1); i++ ){

    dummy = REAL( density_matrix_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ] );

    fprintf( populations_file_p, DOUBLE_FORMAT"  ", dummy );

  }

  dummy = REAL( density_matrix_p->matrix[ ELECTRON_MANY_INDEX( i, i ) ] );

  fprintf( populations_file_p, DOUBLE_FORMAT"\n", dummy );

  fflush( populations_file_p );


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


  fprintf( energies_file_p, "%le  ", time  );

  for( i=0; i<4; i++ ){

    fprintf( energies_file_p, DOUBLE_FORMAT"  ", kinetic_energy_system_p->rvector[ i ]  );

  }
  
  for( i=0; i<4; i++ ){

    fprintf( energies_file_p, DOUBLE_FORMAT"  ", potential_energy_system_p->rvector[ i ]  );

  }

  if( N_chain ){


    for( i=0; i<(N_chain +1); i++ ){

      fprintf( energies_file_p, DOUBLE_FORMAT"  ", kinetic_energy_thermostat_p->rvector[ i ]  );

    }
  
    for( i=0; i<(N_chain +1); i++ ){

      fprintf( energies_file_p, DOUBLE_FORMAT"  ", potential_energy_thermostat_p->rvector[ i ]  );

    }

    fprintf( energies_file_p, DOUBLE_FORMAT"  ", total_energy_system );

    fprintf( energies_file_p, DOUBLE_FORMAT"\n", pseudo_energy_system);

  }
  else{ 
  
    fprintf( energies_file_p, DOUBLE_FORMAT"\n", total_energy_system);

  }

  fflush( energies_file_p );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_mu_trace( const constants constants, const state state, const config config ){

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


  fprintf( mu_trace_file_p, "%le  ", time );

  //
  dummy =  REAL( MATRIX_TRACE( *mu00_p ) );

  fprintf( mu_trace_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;
 
  for( i=0; i<N_coor; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu01_p->array[ i ] ) );

  }

  fprintf( mu_trace_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<N_coor; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu10_p->array[ i ] ) );

  }

  fprintf( mu_trace_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  dim = N_coor *N_coor;

  for( i=0; i<dim; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu02_p->array[ i ] ) );

  }

  fprintf( mu_trace_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<dim; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu11_p->array[ i ] ) );

  }

  fprintf( mu_trace_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<dim; i++ ){

    dummy +=  REAL( MATRIX_TRACE( mu20_p->array[ i ] ) );

  }

  fprintf( mu_trace_file_p, DOUBLE_FORMAT"\n", dummy );

  fflush( mu_trace_file_p );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_mu_norm( const constants constants, const state state, const config config ){

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


  fprintf( mu_norm_file_p, "%le  ", time );
  
  //
  dummy =  MATRIX_NORM( *mu00_p );

  fprintf( mu_norm_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;
 
  for( i=0; i<N_coor; i++ ){

    dummy +=  MATRIX_NORM( mu01_p->array[ i ] );

  }

  fprintf( mu_norm_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<N_coor; i++ ){

    dummy +=  MATRIX_NORM( mu10_p->array[ i ] );

  }

  fprintf( mu_norm_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  dim = N_coor *N_coor;

  for( i=0; i<dim; i++ ){

    dummy +=  MATRIX_NORM( mu02_p->array[ i ] );

  }

  fprintf( mu_norm_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<dim; i++ ){

    dummy +=  MATRIX_NORM( mu11_p->array[ i ] );

  }

  fprintf( mu_norm_file_p, DOUBLE_FORMAT"  ", dummy );

  //
  dummy=0.0e0;

  for( i=0; i<dim; i++ ){

    dummy +=  MATRIX_NORM( mu20_p->array[ i ] );

  }

  fprintf( mu_norm_file_p, DOUBLE_FORMAT"\n", dummy );

  fflush( mu_norm_file_p );


  return info;

}

//------------------------------------------

#ifdef __DEBUG_PLUS__ 

//BUGFIX: this was originally DEF
int print_rho_trace( const constants constants, const state state, const config config ){

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


  fprintf( rho_trace_file_p, "%le  ", time );

  for( i=0; i<(max_rho_index-1); i++){

    dummy =  REAL( MATRIX_TRACE( rho_p->array[ i ] ) );

    fprintf( rho_trace_file_p, DOUBLE_FORMAT"  ", dummy );

    dummy =  IMAG( MATRIX_TRACE( rho_p->array[ i ] ) );

    fprintf( rho_trace_file_p, DOUBLE_FORMAT"  ", dummy );

  }

  dummy =  REAL( MATRIX_TRACE( rho_p->array[ i ] ) );

  fprintf( rho_trace_file_p, DOUBLE_FORMAT"  ", dummy );

  dummy =  IMAG( MATRIX_TRACE( rho_p->array[ i ] ) );

  fprintf( rho_trace_file_p, DOUBLE_FORMAT"\n", dummy );


  fflush( rho_trace_file_p );


  return info;

}

//------------------------------------------

//BUGFIX: this was originally DEF
int print_rho_norm( const constants constants, const state state, const config config ){


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


  fprintf( rho_norm_file_p, "%le  ", time );

  for( i=0; i<(max_rho_index-1); i++){

    dummy = MATRIX_NORM( rho_p->array[ i ] );

    fprintf( rho_norm_file_p, DOUBLE_FORMAT"  ", dummy );

  }

  dummy = MATRIX_NORM( rho_p->array[ i ] );

  fprintf( rho_norm_file_p, DOUBLE_FORMAT"\n", dummy );


  fflush( rho_norm_file_p );


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


  fprintf( adiabatic_populations_file_p, "%le  ", time );


  for( i=0; i<N_levels_many-1; i++ ){

    fprintf( adiabatic_populations_file_p, DOUBLE_FORMAT"  ", adiabatic_populations_p->rvector[ i ] );

  }
  fprintf( adiabatic_populations_file_p, DOUBLE_FORMAT"\n", adiabatic_populations_p->rvector[ i ] );

  fflush( adiabatic_populations_file_p );


  return info;

}

//------------------------------------------
//
int print_adiabatic_projection( const constants constants, const state state ){

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


  fprintf( adiabatic_projection_file_p, "# time = %le\n", time );

  //
  fprintf( adiabatic_projection_file_p, "# eigenvalues:  " );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( adiabatic_projection_file_p, DOUBLE_FORMAT"  ", electronic_density_eigenvalues_p->rvector[ N_levels_many -1 -i ] );

  }

  fprintf( adiabatic_projection_file_p, "\n" );

  //
  fprintf( adiabatic_projection_file_p, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( adiabatic_projection_file_p, "#%d  ", i+1 );

    for( j=0; j<N_levels_many; j++ ){

      index = ELECTRON_MANY_INDEX( j, N_levels_many -1 -i );

      dummy = CMPLX_NORM( adiabatic_projection_p->matrix[ index ] );

      fprintf( adiabatic_projection_file_p, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( adiabatic_projection_file_p, "\n" );

  } /* i loop */

  //
  fprintf( adiabatic_projection_file_p, "\n\n" );

  fflush( adiabatic_projection_file_p );


  return info;

}

//------------------------------------------

int print_projection( const state state ){

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


  fprintf( projection_file_p, "%le  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"\n", time, rho_norm, rho_dot_norm, rho_projection );

  fflush( projection_file_p );


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


  fprintf( adiabatic_PES_many_file_p, "%le  ", time );

  for( i=0; i<N_levels_many-1; i++ ){

    fprintf( adiabatic_PES_many_file_p, DOUBLE_FORMAT"  ", adiabatic_PES_many_p->rvector[ i ] );

  }

  fprintf( adiabatic_PES_many_file_p, DOUBLE_FORMAT"\n", adiabatic_PES_many_p->rvector[ i ] );

  fflush( adiabatic_PES_many_file_p );


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


  fprintf( adiabatic_PES_single_file_p, "%le  ", time );

  for( i=0; i<N_levels_single-1; i++ ){

    fprintf( adiabatic_PES_single_file_p, DOUBLE_FORMAT"  ", adiabatic_PES_single_p->rvector[ i ] );

  }

  fprintf( adiabatic_PES_single_file_p, DOUBLE_FORMAT"\n", adiabatic_PES_single_p->rvector[ i ] );

  fflush( adiabatic_PES_single_file_p );


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


  fprintf( single_level_populations_Ehrenfest_file_p, "%le  ", time );

  for( i=0; i<N_levels_single-1; i++ ){

    fprintf( single_level_populations_Ehrenfest_file_p, DOUBLE_FORMAT"  ", single_level_populations_Ehrenfest_p->rvector[ i ] );

  }

  fprintf( single_level_populations_Ehrenfest_file_p, DOUBLE_FORMAT"\n", single_level_populations_Ehrenfest_p->rvector[ i ] );

  fflush( single_level_populations_Ehrenfest_file_p );


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


  fprintf( single_level_populations_adiabatic_file_p, "%le  ", time );

  for( i=0; i<N_levels_single-1; i++ ){

    fprintf( single_level_populations_adiabatic_file_p, DOUBLE_FORMAT"  ", single_level_populations_adiabatic_p->rvector[ i ] );

  }

  fprintf( single_level_populations_adiabatic_file_p, DOUBLE_FORMAT"\n", single_level_populations_adiabatic_p->rvector[ i ] );

  fflush( single_level_populations_adiabatic_file_p );


  return info;

}

//------------------------------------------

int print_nonadiabatic_coupling( const constants constants, const state state ){

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


  fprintf( nonadiabatic_coupling_file_p, "%le  ", time );

  for( i=0; i<N_coor-1; i++ ){

    fprintf( nonadiabatic_coupling_file_p, DOUBLE_FORMAT"  ", nonadiabatic_coupling_p->rvector[ i ] );

  }

  fprintf( nonadiabatic_coupling_file_p, DOUBLE_FORMAT"\n", nonadiabatic_coupling_p->rvector[ i ] );

  fflush( nonadiabatic_coupling_file_p );


  return info;

}

//------------------------------------------

int print_nonadiabatic_rate( const constants constants, const state state ){

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


  fprintf( nonadiabatic_rate_file_p, "%le  ", time );

  for( i=0; i<N_coor-1; i++ ){

    fprintf( nonadiabatic_rate_file_p, DOUBLE_FORMAT"  ", nonadiabatic_rate_p->rvector[ i ] );

  }

  fprintf( nonadiabatic_rate_file_p, DOUBLE_FORMAT"\n", nonadiabatic_rate_p->rvector[ i ] );

  fflush( nonadiabatic_rate_file_p );


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


  fprintf( one_body_electronic_density_matrix_file_p, "# time = %le\n", time );

  for( i=0; i<N_levels_single; i++ ){

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( i, j );

      fprintf( one_body_electronic_density_matrix_file_p, DOUBLE_FORMAT"  ", CMPLX_NORM( one_body_electronic_density_matrix_p->matrix[ index ]) );

    } /* j loop */

    fprintf( one_body_electronic_density_matrix_file_p, "\n" );

  } /* i loop */

  //
  dummy = REAL( MATRIX_TRACE( *one_body_electronic_density_matrix_p ) );

  fprintf( one_body_electronic_density_matrix_file_p, "# trace = %le  ", dummy );

  dummy = IMAG( MATRIX_TRACE( *one_body_electronic_density_matrix_p ) );

  fprintf( one_body_electronic_density_matrix_file_p, "%le\n", dummy );


  //
  dummy  = MATRIX_NORM( *one_body_electronic_density_matrix_p );

  fprintf( one_body_electronic_density_matrix_file_p, "# norm = %le\n", dummy );

  //
  fprintf( one_body_electronic_density_matrix_file_p, "\n\n" );

  fflush( one_body_electronic_density_matrix_file_p );


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


  fprintf( one_body_electronic_hole_matrix_file_p, "# time = %le\n", time );

  for( i=0; i<N_levels_single; i++ ){

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( i, j );

      /* This is a tentative spectral function of the transition |GS> -->|instantaneous> */
      fprintf( one_body_electronic_hole_matrix_file_p, DOUBLE_FORMAT"  ", CMPLX_NORM( one_body_electronic_hole_matrix_p->matrix[ index ] ) );

    } /* j loop */

    fprintf( one_body_electronic_hole_matrix_file_p, "\n" );

  } /* i loop */

  /*
  dummy = REAL( MATRIX_TRACE( *one_body_electronic_hole_matrix_p ) );

  fprintf( one_body_electronic_hole_matrix_file_p, "# trace = %le ", dummy );

  dummy = IMAG( MATRIX_TRACE( *one_body_electronic_hole_matrix_p ) );

  fprintf( one_body_electronic_hole_matrix_file_p, "%le \n", dummy );

  dummy  = MATRIX_NORM( *one_body_electronic_hole_matrix_p );

  fprintf( one_body_electronic_hole_matrix_file_p, "# norm = %le \n", dummy );
  */

  //
  fprintf( one_body_electronic_hole_matrix_file_p, "\n\n" );

  fflush( one_body_electronic_hole_matrix_file_p );


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


  fprintf( one_body_electronic_particle_matrix_file_p, "# time = %le\n", time );

  for( i=0; i<N_levels_single; i++ ){

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( i, j );

      /* This is a tentative spectral function of the transition |GS> -->|instantaneous> */
      fprintf( one_body_electronic_particle_matrix_file_p, DOUBLE_FORMAT"  ", CMPLX_NORM( one_body_electronic_particle_matrix_p->matrix[ index ] ) );

    } /* j loop */

    fprintf( one_body_electronic_particle_matrix_file_p, "\n" );

  } /* i loop */

  /*
  dummy = REAL( MATRIX_TRACE( *one_body_electronic_particle_matrix_p ) );

  fprintf( one_body_electronic_particle_matrix_file_p, "# trace = %le ", dummy );

  dummy = IMAG( MATRIX_TRACE( *one_body_electronic_particle_matrix_p ) );

  fprintf( one_body_electronic_particle_matrix_file_p, "%le \n", dummy );

  dummy  = MATRIX_NORM( *one_body_electronic_particle_matrix_p );

  fprintf( one_body_electronic_particle_matrix_file_p, "# norm = %le \n", dummy );
  */

  //
  fprintf( one_body_electronic_particle_matrix_file_p, "\n\n" );

  fflush( one_body_electronic_particle_matrix_file_p );


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


  fprintf( natural_orbitals_file_p, "# time = %le\n", time );

  //
  fprintf( natural_orbitals_file_p, "# eigenvalues:  " );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( natural_orbitals_file_p, DOUBLE_FORMAT"  ", dummy_rvector_single_p->rvector[ N_levels_single -1 -i ] );

  }

  fprintf( natural_orbitals_file_p, "\n" );

  //
  fprintf( natural_orbitals_file_p, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( natural_orbitals_file_p, "#%d   ", i+1 );

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( j, N_levels_single -1 -i );  // WARNING: note the reverse ordering

      dummy = CMPLX_NORM( dummy_matrix_single1_p->matrix[ index ] );

      fprintf( natural_orbitals_file_p, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( natural_orbitals_file_p, "\n" );

  } /* i loop */

  //
  fprintf( natural_orbitals_file_p, "\n\n" );

  fflush( natural_orbitals_file_p );


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


  fprintf( hole_orbitals_file_p, "# time = %le\n", time );

  //
  fprintf( hole_orbitals_file_p, "# eigenvalues:  " );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( hole_orbitals_file_p, DOUBLE_FORMAT"  ", hole_populations_p->rvector[ N_levels_single -1 -i ] );

  }

  fprintf( hole_orbitals_file_p, "\n" );

  //
  fprintf( hole_orbitals_file_p, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( hole_orbitals_file_p, "#%d  ", i+1 );

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( j, N_levels_single -1 -i );  // WARNING: note the reverse ordering

      dummy = CMPLX_NORM( hole_orbitals_p->matrix[ index ] );

      fprintf( hole_orbitals_file_p, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( hole_orbitals_file_p, "\n" );

  } /* i loop */

  //
  fprintf( hole_orbitals_file_p, "\n\n" );

  fflush( hole_orbitals_file_p );


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


  fprintf( particle_orbitals_file_p, "# time = %le\n", time );

  //
  fprintf( particle_orbitals_file_p, "# eigenvalues:  " );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( particle_orbitals_file_p, DOUBLE_FORMAT"  ", particle_populations_p->rvector[ N_levels_single -1 -i ] );

  }

  fprintf( particle_orbitals_file_p, "\n" );

  //
  fprintf( particle_orbitals_file_p, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_single; i++ ){

    fprintf( particle_orbitals_file_p, "#%d  ", i+1 );

    for( j=0; j<N_levels_single; j++ ){

      index = ELECTRON_SINGLE_INDEX( j, N_levels_single -1 -i );  // WARNING: note the reverse ordering

      dummy = CMPLX_NORM( particle_orbitals_p->matrix[ index ] );

      fprintf( particle_orbitals_file_p, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( particle_orbitals_file_p, "\n" );

  } /* i loop */

  //
  fprintf( particle_orbitals_file_p, "\n\n" );

  fflush( particle_orbitals_file_p );


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


  fprintf( adiabatic_states_file_p, "# time = %le\n", time );

  //
  fprintf( adiabatic_states_file_p, "# eigenvalues:  " );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( adiabatic_states_file_p, DOUBLE_FORMAT"  ", adiabatic_PES_many_p->rvector[ i ] );

  }

  fprintf( adiabatic_states_file_p, "\n" );

  //
  fprintf( adiabatic_states_file_p, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( adiabatic_states_file_p, "#%d  ", i+1 );

    for( j=0; j<N_levels_many; j++ ){

      index = ELECTRON_MANY_INDEX( j, i );

      dummy = CMPLX_NORM( adiabatic_states_many_p->matrix[ index ] );

      fprintf( adiabatic_states_file_p, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( adiabatic_states_file_p, "\n" );

  } /* i loop */

  //
  fprintf( adiabatic_states_file_p, "\n\n" );

  fflush( adiabatic_states_file_p );


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


  fprintf( electronic_density_states_file_p, "# time = %le\n", time );

  //
  fprintf( electronic_density_states_file_p, "# eigenvalues:  " );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( electronic_density_states_file_p, DOUBLE_FORMAT"  ", electronic_density_eigenvalues_p->rvector[ N_levels_many -1 -i ] );

  }

  fprintf( electronic_density_states_file_p, "\n" );

  //
  fprintf( electronic_density_states_file_p, "# eigenvectors squared =\n" );

  for( i=0; i<N_levels_many; i++ ){

    fprintf( electronic_density_states_file_p, "#%d  ", i+1 );

    for( j=0; j<N_levels_many; j++ ){

      index = ELECTRON_MANY_INDEX( j, N_levels_many -1 -i );

      dummy = CMPLX_NORM( electronic_density_eigenvectors_p->matrix[ index ] );

      fprintf( electronic_density_states_file_p, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( electronic_density_states_file_p, "\n" );

  } /* i loop */

  //
  fprintf( electronic_density_states_file_p, "\n\n" );

  fflush( electronic_density_states_file_p );


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


  fprintf( ionic_density_states_file_p, "# time = %le\n", time );

  //
  fprintf( ionic_density_states_file_p, "# eigenvalues:  " );

  for( i=0; i<sqrt_max_rho_index; i++ ){

    fprintf( ionic_density_states_file_p, DOUBLE_FORMAT"  ", ionic_density_eigenvalues_p->rvector[ sqrt_max_rho_index -1 -i ] );

  }

  fprintf( ionic_density_states_file_p, "\n" );

  //
  fprintf( ionic_density_states_file_p, "# eigenvectors squared =\n" );

  for( i=0; i<sqrt_max_rho_index; i++ ){

    fprintf( ionic_density_states_file_p, "#%d  ", i+1 );

    for( j=0; j<sqrt_max_rho_index; j++ ){

      index = RHO_INDEX_AUX( j, sqrt_max_rho_index -1 -i );

      dummy = CMPLX_NORM( ionic_density_eigenvectors_p->matrix[ index ] );

      fprintf( ionic_density_states_file_p, DOUBLE_FORMAT"  ", dummy *dummy );

    } /* j loop */

    fprintf( ionic_density_states_file_p, "\n" );

  } /* i loop */

  //
  fprintf( ionic_density_states_file_p, "\n\n" );

  fflush( ionic_density_states_file_p );


  return info;

}

//------------------------------------------

int print_dipole_many( const constants constants, const state state ){

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


  fprintf( dipole_many_file_p, "%le", time );


  index=0;

  for( i=0; i<N_levels_many; i++ ){

    for( j=(i+1); j<N_levels_many; j++ ){
    
       // index_ratio = i; // Absorption
       index_ratio = j; // Emission

       ratio = adiabatic_populations_p->rvector[ index_ratio ];

       fprintf( dipole_many_file_p, "  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"", oscillator_frequency_many_p->rvector[ index ], ratio *oscillator_strength_many_p->rvector[ index ], oscillator_strength_many_p->rvector[ index ] );

       index++;

     } /* end loop j */  

  } /* end loop i */

  fprintf( dipole_many_file_p, "\n" );


  return info;

}

//------------------------------------------

int print_dipole_single( const constants constants, const state state ){

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


  fprintf( dipole_single_file_p, "%le  ", time );

  for( i=0; i<(N_levels_single-1); i++ ){
                  
    ratio = ( hole_populations_p->rvector[ N_levels_single -1 -i ] );

    fprintf( dipole_single_file_p, DOUBLE_FORMAT"  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"  ", oscillator_frequency_single_p->rvector[ N_levels_single -1 -i ], ratio *oscillator_strength_single_p->rvector[ N_levels_single -1 -i ], oscillator_strength_single_p->rvector[ N_levels_single -1 -i ] );

  }

  ratio = ( hole_populations_p->rvector[ N_levels_single -1 -i ] );

  fprintf( dipole_single_file_p, DOUBLE_FORMAT"  "DOUBLE_FORMAT"  "DOUBLE_FORMAT"\n", oscillator_frequency_single_p->rvector[ N_levels_single -1 -i ], ratio *oscillator_strength_single_p->rvector[ N_levels_single -1 -i ], oscillator_strength_single_p->rvector[ N_levels_single -1 -i ] );


  fflush( dipole_single_file_p );


  return info;

}

//------------------------------------------

//BUGFIX: thsi was originally DEF
int print_position_ion_frame( const constants constants, const state state, const config config ){

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

  fprintf( xyz_file_p, "%d\n", N_atoms );
  fprintf( xyz_file_p, "# time = %le\n", time );

  if( !info ){

    for( i_atom=0; i_atom<N_atoms; i_atom++ ){

      fprintf( xyz_file_p, "%s  "DOUBLE_FORMAT"  ", "C", positions_p->rvector[ i_atom *sdim ] );
      
      if( sdim > 1 ){

	fprintf( xyz_file_p, DOUBLE_FORMAT"  ", positions_p->rvector[ 1 +i_atom *sdim ] );

      }
      else{
      
	fprintf( xyz_file_p, DOUBLE_FORMAT"  ", 0.0e0  );
      
      }

      if( sdim > 2 ){

	fprintf( xyz_file_p, DOUBLE_FORMAT"\n", positions_p->rvector[ 2 +i_atom *sdim ] );
      
      }
      else{
      
	fprintf( xyz_file_p, DOUBLE_FORMAT"\n", 0.0e0  );
      
      }

    } /* end for loop */

  } /* end info if*/

  fflush(  xyz_file_p );


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int output_files_opening( const constants constants ){

  /* constants */
  int       N_levels_many;
  int       N_chain;
  /* dummies */
  char      buffer[MAX_STRING_LENGTH];
  char*     output_label;
  int       info=0;


  N_levels_many = constants.N_levels_many;
  N_chain       = constants.N_chain;
  output_label  = (char*)constants.output_label;


  /* xyz file */
  strcpy(  buffer, output_label );
  strcat(  buffer, ".xyz" );

  xyz_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( xyz_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }


  /* positions file */
  strcpy(  buffer, "positions_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  positions_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( positions_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( positions_file_p, constants ) ) info=1;


  /* momenta file */
  strcpy(  buffer, "momenta_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  momenta_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( momenta_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( momenta_file_p, constants ) ) info=1;


  /* forces file */
  strcpy(  buffer, "forces_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  forces_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( forces_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( forces_file_p, constants ) ) info=1;


  if( N_chain ){

    /* positions thermostat file */
    strcpy(  buffer, "positions_thermostat_" );
    strcat(  buffer, output_label );
    strcat(  buffer, ".dat" );

    positions_thermostat_file_p = fopen( buffer, "wt");
    if( FILE_CHECK( positions_thermostat_file_p, output_files_opening ) ){
      fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
      fflush( stderr);
      info=1;
    }

    // if( CONSTANTS_VERBOSE_PRINT( positions_thermostat_file_p, constants ) ) info=1;


    /* momenta thermostat file */
    strcpy(  buffer, "momenta_thermostat_" );
    strcat(  buffer, output_label );
    strcat(  buffer, ".dat" );

    momenta_thermostat_file_p = fopen( buffer, "wt");
    if( FILE_CHECK( momenta_thermostat_file_p, output_files_opening ) ){
      fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
      fflush( stderr);
      info=1;
    }

    // if( CONSTANTS_VERBOSE_PRINT( momenta_thermostat_file_p, constants ) ) info=1;


    /* forces thermostat file */
    strcpy(  buffer, "forces_thermostat_" );
    strcat(  buffer, output_label );
    strcat(  buffer, ".dat" );

    forces_thermostat_file_p = fopen( buffer, "wt");
    if( FILE_CHECK( forces_thermostat_file_p, output_files_opening ) ){
      fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
      fflush( stderr);
      info=1;
    }

   // if( CONSTANTS_VERBOSE_PRINT( forces_thermostat_file_p, constants ) ) info=1;

  }


  /* populations file */
  strcpy(  buffer, "populations_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  populations_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( populations_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( populations_file_p, constants ) ) info=1;


  /* energies file */
  strcpy(  buffer, "energies_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  energies_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( energies_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( energies_file_p, constants ) ) info=1;


  /* mu trace file */
  strcpy(  buffer, "mu_trace_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  mu_trace_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( mu_trace_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( mu_trace_file_p, constants ) ) info=1;


  /* mu norm file */
  strcpy(  buffer, "mu_norm_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  mu_norm_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( mu_norm_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( mu_norm_file_p, constants ) ) info=1;

#ifdef __DEBUG_PLUS__


  /* rho trace file */
  strcpy(  buffer, "rho_trace_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  rho_trace_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( rho_trace_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( rho_trace_file_p, constants ) ) info=1;


  /* rho norm file */
  strcpy(  buffer, "rho_norm_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  rho_norm_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( rho_norm_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( rho_norm_file_p, constants ) ) info=1;

#endif /* __DEBUG_PLUS__ */


  /* adiabatic_populations file */
  strcpy(  buffer, "adiabatic_populations_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  adiabatic_populations_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( adiabatic_populations_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( adiabatic_populations_file_p, constants ) ) info=1;


  /* adiabatic_projection file */
  strcpy(  buffer, "adiabatic_projection_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  adiabatic_projection_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( adiabatic_projection_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( adiabatic_projection_file_p, constants ) ) info=1;


  /* projection file */
  strcpy(  buffer, "projection_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  projection_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( projection_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( projection_file_p, constants ) ) info=1;


  /* adiabatic_PES_many file */
  strcpy(  buffer, "adiabatic_PES_many_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  adiabatic_PES_many_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( adiabatic_PES_many_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( adiabatic_PES_many_file_p, constants ) ) info=1;


  /* adiabatic_PES_single file */
  strcpy(  buffer, "adiabatic_PES_single_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  adiabatic_PES_single_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( adiabatic_PES_single_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( adiabatic_PES_single_file_p, constants ) ) info=1;


  /* single_level_populations_Ehrenfest file */
  strcpy(  buffer, "single_level_populations_Ehrenfest_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  single_level_populations_Ehrenfest_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( single_level_populations_Ehrenfest_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( single_level_populations_Ehrenfest_file_p, constants ) ) info=1;


  /* single_level_populations_adiabatic file */
  strcpy(  buffer, "single_level_populations_adiabatic_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  single_level_populations_adiabatic_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( single_level_populations_adiabatic_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( single_level_populations_adiabatic_file_p, constants ) ) info=1;


  /* nonadiabatic_coupling file */
  strcpy(  buffer, "nonadiabatic_coupling_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  nonadiabatic_coupling_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( nonadiabatic_coupling_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( nonadiabatic_coupling_file_p, constants ) ) info=1;


  /* nonadiabatic_rate file */
  /*
  strcpy(  buffer, "nonadiabatic_rate_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  nonadiabatic_rate_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( nonadiabatic_rate_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  if( CONSTANTS_VERBOSE_PRINT( nonadiabatic_rate_file_p, constants ) ) info=1;
  */

  /* one_body_electronic_density_matrix file */
  strcpy(  buffer, "one_body_electronic_density_matrix_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  one_body_electronic_density_matrix_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( one_body_electronic_density_matrix_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( one_body_electronic_density_matrix_file_p, constants ) ) info=1;
  

  /* one_body_electronic_hole_matrix file */
  strcpy(  buffer, "one_body_electronic_hole_matrix_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  one_body_electronic_hole_matrix_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( one_body_electronic_hole_matrix_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( one_body_electronic_hole_matrix_file_p, constants ) ) info=1;
  

  /* one_body_electronic_particle_matrix file */
  strcpy(  buffer, "one_body_electronic_particle_matrix_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  one_body_electronic_particle_matrix_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( one_body_electronic_particle_matrix_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( one_body_electronic_particle_matrix_file_p, constants ) ) info=1;
  

  /* natural_orbitals file */
  strcpy(  buffer, "natural_orbitals_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  natural_orbitals_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( natural_orbitals_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( natural_orbitals_file_p, constants ) ) info=1;


  /* hole_orbitals file */
  strcpy(  buffer, "hole_orbitals_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  hole_orbitals_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( hole_orbitals_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( hole_orbitals_file_p, constants ) ) info=1;


  /* particle_orbitals file */
  strcpy(  buffer, "particle_orbitals_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  particle_orbitals_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( particle_orbitals_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( particle_orbitals_file_p, constants ) ) info=1;


  /* adiabatic_states file */
  strcpy(  buffer, "adiabatic_states_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  adiabatic_states_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( adiabatic_states_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( adiabatic_states_file_p, constants ) ) info=1;


  /* electronic_density_states file */
  strcpy(  buffer, "electronic_density_states_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  electronic_density_states_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( electronic_density_states_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( electronic_density_states_file_p, constants ) ) info=1;


  /* ionic_density_states file */
  strcpy(  buffer, "ionic_density_states_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  ionic_density_states_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( ionic_density_states_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( ionic_density_states_file_p, constants ) ) info=1;


  if( N_levels_many > 1 ){

    /* dipole many file */
    strcpy(  buffer, "dipole_many_" );
    strcat(  buffer, output_label );
    strcat(  buffer, ".dat" );

    dipole_many_file_p = fopen( buffer, "wt");
    if( FILE_CHECK( dipole_many_file_p, output_files_opening ) ){
      fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
      fflush( stderr);
      info=1;
    }

    // if( CONSTANTS_VERBOSE_PRINT( dipole_many_file_p, constants ) ) info=1;

  }


  /* dipole single file */
  strcpy(  buffer, "dipole_single_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  dipole_single_file_p = fopen( buffer, "wt");
  if( FILE_CHECK( dipole_single_file_p, output_files_opening ) ){
    fprintf( stderr, "ERROR occurred when opening file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }

  // if( CONSTANTS_VERBOSE_PRINT( dipole_single_file_p, constants ) ) info=1;


  return info;

}

//------------------------------------------

int output_files_closing( const constants constants ){

  /* constants */
  int       N_levels_many;
  int       N_chain;
  /* dummies */
  char      buffer[MAX_STRING_LENGTH];
  char*     output_label;
  int       info=0;

                          
  N_levels_many = constants.N_levels_many;
  N_chain       = constants.N_chain;
  output_label  = (char*)constants.output_label;


  /* xyz file */
  strcpy(  buffer, output_label );
  strcat(  buffer, ".xyz" );

  if( FILE_CHECK( xyz_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( xyz_file_p );


  /* positions file */
  strcpy(  buffer, "positions_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( positions_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( positions_file_p );


  /* momenta file */
  strcpy(  buffer, "momenta_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( momenta_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( momenta_file_p );


  /* forces file */
  strcpy(  buffer, "forces_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( forces_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( forces_file_p );


  if( N_chain ){

    /* positions thermostat file */
    strcpy(  buffer, "positions_thermostat_" );
    strcat(  buffer, output_label );
    strcat(  buffer, ".dat" );

    if( FILE_CHECK( positions_thermostat_file_p, output_files_closing ) ){
      fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
      fflush( stderr);
      info=1;
    }
    fclose( positions_thermostat_file_p );


    /* momenta thermostat file */
    strcpy(  buffer, "momenta_thermostat_" );
    strcat(  buffer, output_label );
    strcat(  buffer, ".dat" );

    if( FILE_CHECK( momenta_thermostat_file_p, output_files_closing ) ){
      fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
      fflush( stderr);
      info=1;
    }
    fclose( momenta_thermostat_file_p );


    /* forces file */
    strcpy(  buffer, "forces_thermostat_" );
    strcat(  buffer, output_label );
    strcat(  buffer, ".dat" );

    if( FILE_CHECK( forces_thermostat_file_p, output_files_closing ) ){
      fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
      fflush( stderr);
      info=1;
    }
    fclose( forces_thermostat_file_p );

  }


  /* populations file */
  strcpy(  buffer, "populations_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( populations_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( populations_file_p );


  /* energies file */
  strcpy(  buffer, "energies_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( energies_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( energies_file_p );


  /* mu trace file */
  strcpy(  buffer, "mu_trace" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( mu_trace_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( mu_trace_file_p );


  /* mu norm file */
  strcpy(  buffer, "mu_norm" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( mu_norm_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( mu_norm_file_p );

#ifdef __DEBUG_PLUS__

  /* rho trace file */
  strcpy(  buffer, "rho_trace" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( rho_trace_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( rho_trace_file_p );


  /* rho norm file */
  strcpy(  buffer, "rho_norm" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( rho_norm_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( rho_norm_file_p );

#endif /* __DEBUG_PLUS__ */


  /* adiabatic_populations file */
  strcpy(  buffer, "adiabatic_populations" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( adiabatic_populations_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( adiabatic_populations_file_p );


  /* adiabatic_projection file */
  strcpy(  buffer, "adiabatic_projection" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( adiabatic_projection_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( adiabatic_projection_file_p );


  /* projection file */
  strcpy(  buffer, "projection" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( projection_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( projection_file_p );


  /* adiabatic_PES_many file */
  strcpy(  buffer, "adiabatic_PES_many_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( adiabatic_PES_many_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( adiabatic_PES_many_file_p );


  /* adiabatic_PES_single file */
  strcpy(  buffer, "adiabatic_PES_single_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( adiabatic_PES_single_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( adiabatic_PES_single_file_p );


  /* single_level_populations_Ehrenfest file */
  strcpy(  buffer, "single_level_populations_Ehrenfest_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( single_level_populations_Ehrenfest_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( single_level_populations_Ehrenfest_file_p );


  /* single_level_populations_adiabatic file */
  strcpy(  buffer, "single_level_populations_adiabatic_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( single_level_populations_adiabatic_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( single_level_populations_adiabatic_file_p );


  /* nonadiabatic_coupling file */
  strcpy(  buffer, "nonadiabatic_coupling_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( nonadiabatic_coupling_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( nonadiabatic_coupling_file_p );


  /* nonadiabatic_rate file */
  /*
  strcpy(  buffer, "nonadiabatic_rate_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( nonadiabatic_rate_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( nonadiabatic_rate_file_p );
  */

  /* one_body_electronic_density_matrix file */
  strcpy(  buffer, "one_body_electronic_density_matrix_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( one_body_electronic_density_matrix_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( one_body_electronic_density_matrix_file_p );


  /* one_body_electronic_hole_matrix file */
  strcpy(  buffer, "one_body_electronic_hole_matrix_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( one_body_electronic_hole_matrix_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( one_body_electronic_hole_matrix_file_p );


  /* one_body_electronic_particle_matrix file */
  strcpy(  buffer, "one_body_electronic_particle_matrix_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( one_body_electronic_particle_matrix_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( one_body_electronic_particle_matrix_file_p );


  /* natural_orbitals file */
  strcpy(  buffer, "natural_orbitals_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( natural_orbitals_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( natural_orbitals_file_p );


  /* hole_orbitals file */
  strcpy(  buffer, "hole_orbitals_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( hole_orbitals_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( hole_orbitals_file_p );


  /* particle_orbitals file */
  strcpy(  buffer, "particle_orbitals_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( particle_orbitals_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( particle_orbitals_file_p );


  /* adiabatic_states file */
  strcpy(  buffer, "adiabatic_states_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( adiabatic_states_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( adiabatic_states_file_p );


  /* electronic_density_states file */
  strcpy(  buffer, "electronic_density_states_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( electronic_density_states_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( electronic_density_states_file_p );


  /* ionic_density_states file */
  strcpy(  buffer, "ionic_density_states_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( ionic_density_states_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( ionic_density_states_file_p );


  if( N_levels_many > 1 ){

    /* dipole many file */
    strcpy(  buffer, "dipole_many_" );
    strcat(  buffer, output_label );
    strcat(  buffer, ".dat" );

    if( FILE_CHECK( dipole_many_file_p, output_files_closing ) ){
      fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
      fflush( stderr);
      info=1;
    }
    fclose( dipole_many_file_p );

  }

  
  /* dipole single file */
  strcpy(  buffer, "dipole_single_" );
  strcat(  buffer, output_label );
  strcat(  buffer, ".dat" );

  if( FILE_CHECK( dipole_single_file_p, output_files_closing ) ){
    fprintf( stderr, "ERROR occurred when closing file %s.\n", buffer );
    fflush( stderr);
    info=1;
  }
  fclose( dipole_single_file_p );


  return info;

}

//------------------------------------------
