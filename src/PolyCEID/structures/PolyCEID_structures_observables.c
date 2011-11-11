
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
#include "PolyCEID_structures_observables.h"


/* global */
int is_thermostat_observables=0;
int is_dipole_many=0;


/* ALLOCATE */

int PolyCEID_observables_allocate( int np, int* pp, observables_p observables_p ){

  /* dummies */
  int sdim;
  int N_atoms;
  int N_levels_single;
  int N_levels_many;
  int N_chain;
  int info=0;


  if( np != NP_OBSERVABLES ) info=1;


  N_atoms         = pp[0];
  N_levels_single = pp[1];
  N_levels_many   = pp[2];
  N_chain         = pp[3];
  sdim            = pp[4];

  if( N_chain ){

    is_thermostat_observables = 1;

  }
  else{

    is_thermostat_observables = 0;

  }

  if( N_levels_many > 1 ){

    is_dipole_many = 1;

  }
  else{

    is_dipole_many = 0;

  }

  
  if( RVECTOR_ALLOCATE( 3 *N_atoms, observables_p->kinetic_energy_atoms ) )                                            info=1;

  if( RVECTOR_ALLOCATE( 4, observables_p->kinetic_energy_system ) )                                                    info=1;

  if( RVECTOR_ALLOCATE( 4, observables_p->potential_energy_system ) )                                                  info=1;

  if( is_thermostat_observables ){
	  
    if( RVECTOR_ALLOCATE( N_chain +1, observables_p->kinetic_energy_thermostat ) )                                     info=1;

    if( RVECTOR_ALLOCATE( N_chain +1, observables_p->potential_energy_thermostat ) )                                   info=1;

  }  

  if( is_dipole_many ){

    if( RVECTOR_ALLOCATE( sdim *N_levels_many *(N_levels_many -1) /2, observables_p->oscillator_strength_many ) )               info=1;

    if( RVECTOR_ALLOCATE( sdim *N_levels_many *(N_levels_many -1) /2, observables_p->oscillator_frequency_many ) )              info=1;

  }

  if( RVECTOR_ALLOCATE( N_levels_single, observables_p->oscillator_strength_single ) )                                 info=1;

  if( RVECTOR_ALLOCATE( N_levels_single, observables_p->oscillator_frequency_single ) )                                info=1;


  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_observables_free( observables_p observables_p ){

  /* dummies */
  int info=0;


  if( RVECTOR_FREE( observables_p->oscillator_strength_single ) )                        info=1;

  if( RVECTOR_FREE( observables_p->oscillator_frequency_single ) )                       info=1;

  if( is_dipole_many ){

    if( RVECTOR_FREE( observables_p->oscillator_strength_many ) )                        info=1;

    if( RVECTOR_FREE( observables_p->oscillator_frequency_many ) )                       info=1;

  }        

  if( RVECTOR_FREE( observables_p->potential_energy_system ) )                           info=1;

  if( is_thermostat_observables ){
	  
    if( RVECTOR_FREE( observables_p->kinetic_energy_thermostat ) )                       info=1;

    if( RVECTOR_FREE( observables_p->potential_energy_thermostat ) )                     info=1;

  }  

  if( RVECTOR_FREE( observables_p->kinetic_energy_system ) )                             info=1;

  if( RVECTOR_FREE( observables_p->kinetic_energy_atoms ) )                              info=1;


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_observables_read( FILE* fp, observables_p observables_p ){

  /* dummies */
  double dummy;
  int info=0;


  if( RVECTOR_READ( fp, observables_p->kinetic_energy_atoms ) )    info=1;

  if( RVECTOR_READ( fp, observables_p->kinetic_energy_system ) )   info=1;

  if( RVECTOR_READ( fp, observables_p->potential_energy_system ) ) info=1;

  if( is_thermostat_observables ){
	  
    if( RVECTOR_READ( fp, observables_p->kinetic_energy_thermostat ) )   info=1;

    if( RVECTOR_READ( fp, observables_p->potential_energy_thermostat ) ) info=1;

  }  

  if( fscanf( fp, "%le", &dummy ) < 1 ) info=1;
  observables_p->kinetic_energy_system_correction = dummy;

  if( fscanf( fp, "%le", &dummy ) < 1 ) info=1;
  observables_p->potential_energy_system_correction = dummy;

  if( fscanf( fp, "%le", &dummy ) < 1 ) info=1;
  observables_p->total_energy_system = dummy;

  if( is_thermostat_observables ){
	  
    if( fscanf( fp, "%le", &dummy ) < 1 ) info=1;
    observables_p->pseudo_energy_system = dummy;

  }  

  //
  if( fscanf( fp, "%le", &dummy ) < 1 ) info=1;
  observables_p->superposition_instantaneous_and_initial_state = dummy;

  if( fscanf( fp, "%le", &dummy ) < 1 ) info=1;
  observables_p->superposition_instantaneous_and_excited_state = dummy;

  //
  if( is_dipole_many ){

    if( RVECTOR_READ( fp, observables_p->oscillator_strength_many ) ) info=1;

    if( RVECTOR_READ( fp, observables_p->oscillator_frequency_many ) ) info=1;

  }

  //
  if( RVECTOR_READ( fp, observables_p->oscillator_strength_single ) ) info=1;

  if( RVECTOR_READ( fp, observables_p->oscillator_frequency_single ) ) info=1;


  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_observables_print( FILE* fp, const observables observables ){

  /* dummies */
  int info=0;


  /* kinetic energy of the atoms */
  if( RVECTOR_PRINT( fp, observables.kinetic_energy_atoms ) ) info=1;

  /* kinetic energy of the system */
  if( RVECTOR_PRINT( fp, observables.kinetic_energy_system ) ) info=1;

  /* potential energy of the system */
  if( RVECTOR_PRINT( fp, observables.potential_energy_system ) ) info=1;

  if( is_thermostat_observables ){
  
    /* kinetic energy of the system */
    if( RVECTOR_PRINT( fp, observables.kinetic_energy_thermostat ) ) info=1;

    /* potential energy of the system */
    if( RVECTOR_PRINT( fp, observables.potential_energy_thermostat ) ) info=1;

  }  

  /* correction to the kinetic energy of the system */
  if( fprintf( fp, DOUBLE_FORMAT"\n", observables.kinetic_energy_system_correction ) < 1) info=1;

  /* correction to the potential energy of the system */
  if( fprintf( fp, DOUBLE_FORMAT"\n", observables.potential_energy_system_correction ) < 1 ) info=1;

  /* total energy of the system */
  if( fprintf( fp, DOUBLE_FORMAT"\n", observables.total_energy_system ) < 1) info=1;

  if( is_thermostat_observables ){
  
    /* total energy of the system */
    if( fprintf( fp, DOUBLE_FORMAT"\n", observables.pseudo_energy_system ) < 1) info=1;
  
  }
  
  //
  /* superposition_instantaneous_and_initial_state */
  if( fprintf( fp, DOUBLE_FORMAT"\n", observables.superposition_instantaneous_and_initial_state ) ) info=1;

  /* superposition_instantaneous_and_excited_state */
  if( fprintf( fp, DOUBLE_FORMAT"\n", observables.superposition_instantaneous_and_excited_state ) ) info=1;

  //
  if( is_dipole_many ){

    /* oscillator_strength_many */
    if( RVECTOR_PRINT( fp, observables.oscillator_strength_many ) ) info=1;

    /* oscillator_frequency_many */
    if( RVECTOR_PRINT( fp, observables.oscillator_frequency_many ) ) info=1;

  }

  //
  /* oscillator_strength_many */
  if( RVECTOR_PRINT( fp, observables.oscillator_strength_single ) ) info=1;

  /* oscillator_frequency_many */
  if( RVECTOR_PRINT( fp, observables.oscillator_frequency_single ) ) info=1;


  fflush( fp );


  return info;

}

//------------------------------------------

/* PRINT VERBOSE */

int PolyCEID_observables_verbose_print( FILE* fp, const observables observables ){

  /* dummies */
  int info=0;


  /* kinetic energy of atoms */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# kinetic_energy_atoms:\n" );

  if( RVECTOR_PRINT_PLUS( fp, observables.kinetic_energy_atoms ) ) info=1;

  /* kinetic energy of system */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# kinetic_energy_system:\n" );

  if( RVECTOR_PRINT_PLUS( fp, observables.kinetic_energy_system ) ) info=1;

  /* potential energy of system */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# potential_energy_system:\n" );

  if( RVECTOR_PRINT_PLUS( fp, observables.potential_energy_system ) ) info=1;

  if( is_thermostat_observables ){
  
    /* kinetic energy of the thermostat */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# kinetic_energy_thermostat:\n" );

    if( RVECTOR_PRINT_PLUS( fp, observables.kinetic_energy_thermostat ) ) info=1;

    /* potential energy of the thermostat */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# potential_energy_thermostat:\n" );

    if( RVECTOR_PRINT_PLUS( fp, observables.potential_energy_thermostat ) ) info=1;

  }  

  /* correction to the kinetic energy of system */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# kinetic_energy_system_correction:\n" );

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", observables.kinetic_energy_system_correction ) < 1 ) info=1;

  /* correction to the potential energy of system */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# potential_energy_system_correction:\n" );

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", observables.potential_energy_system_correction ) < 1 ) info=1;

  /* total energy of system */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# total_energy_system:\n" );

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", observables.total_energy_system ) < 1) info=1;

  if( is_thermostat_observables ){
  
    /* pseudo energy of system */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# pseudo_energy_system:\n" );

    if( fprintf( fp, "# "DOUBLE_FORMAT"\n", observables.pseudo_energy_system ) < 1) info=1;

  }

  //
  /* superposition_instantaneous_and_initial_state  */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# superposition_instantaneous_and_initial_state:\n" );

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", observables.superposition_instantaneous_and_initial_state ) < 1) info=1;

  /* superposition_instantaneous_and_excited_state  */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# superposition_instantaneous_and_excited_state:\n" );

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", observables.superposition_instantaneous_and_excited_state ) < 1) info=1;

  //
  if( is_dipole_many ){

    /* oscillator_strength_many */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# oscillator_strength_many:\n" );

    if( RVECTOR_PRINT_PLUS( fp, observables.oscillator_strength_many ) ) info=1;

    //
    /* oscillator_frequency_many */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# oscillator_frequency_many:\n" );

    if( RVECTOR_PRINT_PLUS( fp, observables.oscillator_frequency_many ) ) info=1;

  }  

  //
  /* oscillator_strength_single */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# oscillator_strength_single:\n" );

  if( RVECTOR_PRINT_PLUS( fp, observables.oscillator_strength_single ) ) info=1;

  //
  /* oscillator_frequency_single */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# oscillator_frequency_single:\n" );

  if( RVECTOR_PRINT_PLUS( fp, observables.oscillator_frequency_single ) ) info=1;

  fflush( fp );


  return info;

}

//------------------------------------------

/* COPY [ o1 = o2 ] */

int PolyCEID_observables_copy( observables_p observables_p, const observables observables ){

  /* dummies */
  int info=0;


  if( RVECTOR_COPY( observables_p->kinetic_energy_atoms, observables.kinetic_energy_atoms ) );

  if( RVECTOR_COPY( observables_p->kinetic_energy_system, observables.kinetic_energy_system ) );

  if( RVECTOR_COPY( observables_p->potential_energy_system, observables.potential_energy_system ) );

  if( is_thermostat_observables ){
  
    if( RVECTOR_COPY( observables_p->kinetic_energy_thermostat, observables.kinetic_energy_thermostat ) );

    if( RVECTOR_COPY( observables_p->potential_energy_thermostat, observables.potential_energy_thermostat ) );

  }  

  observables_p->kinetic_energy_system_correction              = observables.kinetic_energy_system_correction;

  observables_p->potential_energy_system_correction            = observables.potential_energy_system_correction;

  observables_p->total_energy_system                           = observables.total_energy_system;

  if( is_thermostat_observables ){
  
    observables_p->pseudo_energy_system                         = observables.pseudo_energy_system;

  }

  //
  observables_p->superposition_instantaneous_and_initial_state = observables.superposition_instantaneous_and_initial_state;

  observables_p->superposition_instantaneous_and_excited_state = observables.superposition_instantaneous_and_excited_state;

  //
  if( is_dipole_many ){

    if( RVECTOR_COPY( observables_p->oscillator_strength_many, observables.oscillator_strength_many ) );

    //
    if( RVECTOR_COPY( observables_p->oscillator_frequency_many, observables.oscillator_frequency_many ) );

  }  

  //
  if( RVECTOR_COPY( observables_p->oscillator_strength_single, observables.oscillator_strength_single ) );

  //
  if( RVECTOR_COPY( observables_p->oscillator_frequency_single, observables.oscillator_frequency_single ) );


  return info;

}

//------------------------------------------
