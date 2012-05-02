
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
#include "PolyCEID_structures_config.h"

/* global */
int is_thermostat_config=0;


/* ALLOCATE */

int PolyCEID_config_allocate( int np, int* pp, config_p config_p ){

  /* config */
  int N_coor;
  int N_levels_single;
  int N_levels_many;
  int max_rho_index;
  int N_atoms;
  int N_chain;
  int sdim;
  /* dummies */
  int pp_atoms[ NP_ATOMS ];
  int pp_thermostat[ NP_THERMOSTAT ];
  int pp_electrons[ NP_ELECTRONS ];
  int info=0;


  if( np != NP_CONFIG ) info=1;

  /* read parameters */
  N_coor             = pp[0];
  N_levels_single    = pp[1];
  N_levels_many      = pp[2];
  max_rho_index      = pp[3];
  N_atoms            = pp[4];
  N_chain            = pp[5];
  sdim               = pp[6];


  /* atoms */
  pp_atoms[0] = N_atoms;
  pp_atoms[1] = sdim; 

  if( ATOMS_ALLOCATE( NP_ATOMS, pp_atoms, config_p->atoms ) ) info=1;
				   

  /* thermostat */
  if( N_chain ){
	  
    is_thermostat_config=1;
    
  }
  else{
	  
    is_thermostat_config=0;

  }

  if( is_thermostat_config ){
	  
    pp_thermostat[0] = N_chain;

    if( THERMOSTAT_ALLOCATE( NP_THERMOSTAT, pp_thermostat, config_p->thermostat ) ) info=1;

  }  


  /* electrons */
  pp_electrons[0] = N_levels_single;
  pp_electrons[1] = N_levels_many;
  pp_electrons[2] = N_coor;
  pp_electrons[3] = max_rho_index;

  if( ELECTRONS_ALLOCATE( NP_ELECTRONS, pp_electrons, config_p->electrons ) ) info=1;
 

  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_config_free( config_p config_p ){

  /* dummies */
  int info=0;


  /* electrons */
  if( ELECTRONS_FREE( config_p->electrons ) ) info=1;


  /* thermostat */
  if( is_thermostat_config ){
	  
    if( THERMOSTAT_FREE( config_p->thermostat ) ) info=1;

  }  


  /* atoms */
  if( ATOMS_FREE( config_p->atoms ) ) info=1;


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_config_read( FILE* fp, config_p config_p ){

  /* dummies */
  int info=0;


  /* time */
  if( fread( ( void* ) &config_p->time, sizeof( config_p->time ), 1, fp ) < 1 ) info=1;

  /* atoms */
  if( ATOMS_READ( fp, config_p->atoms ) )             info=1;


  /* thermostat */
  if( is_thermostat_config ){
	  
    if( THERMOSTAT_READ( fp, config_p->thermostat ) ) info=1;

  }


  /* electrons */
  if( ELECTRONS_READ( fp, config_p->electrons ) )     info=1;
  

  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_config_print( FILE* fp, const config config ){

  /* dummies */
  int info=0;


  /* time */
  if( fwrite( ( const void* ) &config.time, sizeof( config.time ), 1, fp ) < 1 ) info=1;


  /* atoms */
  if( ATOMS_PRINT( fp, config.atoms ) ) info=1;


  if( is_thermostat_config ){
	  
    /* thermostat */
    if( THERMOSTAT_PRINT( fp, config.thermostat ) ) info=1;

  }  


  /* electrons */
  if( ELECTRONS_PRINT( fp, config.electrons ) ) info=1;


  return info;

}

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_config_verbose_print( FILE* fp, const config config ){

  /* dummies */
  int    info=0;


  /* time */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# time:\n");

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", config.time ) < 1 ) info=1;


  /* atoms */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# atoms:\n");

  if( ATOMS_VERBOSE_PRINT( fp, config.atoms ) ) info=1;


  if( is_thermostat_config ){
	  
    /* thermostat */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# thermostat:\n");

    if( THERMOSTAT_VERBOSE_PRINT( fp, config.thermostat ) ) info=1;

  }


  /* electrons */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# electron:\n");

  if( ELECTRONS_VERBOSE_PRINT( fp, config.electrons ) ) info=1;


  fflush( fp );


  return info;

}

//------------------------------------------

/* COPY [ cf1 = cf2 ] */

int PolyCEID_config_copy( config_p config_p, const config config ){

  /* dummies */
  int info=0;


  /* time */
  config_p->time = config.time;


  /* atoms */
  if( ATOMS_COPY( config_p->atoms, config.atoms ) ) info=1;


  /* thermostat  */
  if( is_thermostat_config ){
  
    if( THERMOSTAT_COPY( config_p->thermostat, config.thermostat ) ) info=1;

  }


  /* electrons */
  if( ELECTRONS_COPY( config_p->electrons, config.electrons ) ) info=1;


  return info;

}

//------------------------------------------

/* COMPARE [ cf1 = cf2 ] */

int PolyCEID_config_compare( const config config1, const config config2 ){

  /* dummies */
  int info=0;


  if( ATOMS_COMPARE( config1.atoms, config2.atoms ) ) info=1;

  if( is_thermostat_config ){

    if( THERMOSTAT_COMPARE( config1.thermostat, config2.thermostat ) ) info=1;

  }

  if( ELECTRONS_COMPARE( config1.electrons, config2.electrons ) ) info=1;


  return info;

}  

//------------------------------------------
