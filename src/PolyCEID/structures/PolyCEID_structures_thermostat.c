
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
#include "PolyCEID_structures_thermostat.h"


/* ALLOCATE */

int PolyCEID_thermostat_allocate( int np, int* pp, thermostat_p thermostat_p ){

  /* thermostat */
  int  N_chain;
  /* dummies */
  int  info=0;


  if( np != NP_THERMOSTAT ) info=1;


  N_chain         = pp[0];


  if( RVECTOR_ALLOCATE( N_chain, thermostat_p->positions ) )           info=1;

  if( RVECTOR_ALLOCATE( N_chain, thermostat_p->momenta ) )             info=1;

  if( RVECTOR_ALLOCATE( N_chain, thermostat_p->forces ) )              info=1;


  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_thermostat_free( thermostat_p thermostat_p ){

  /* dummies */
  int info=0;


  if( RVECTOR_FREE( thermostat_p->forces ) )          info=1;

  if( RVECTOR_FREE( thermostat_p->momenta ) )         info=1;

  if( RVECTOR_FREE( thermostat_p->positions ) )       info=1;


  return info;

}

//------------------------------------------

/* COPY [ a1 = a2 ] */

int PolyCEID_thermostat_copy( thermostat_p thermostat_p, const thermostat thermostat ){

  /* dummies */
  int info=0;


  if( RVECTOR_COPY( thermostat_p->positions, thermostat.positions ) )             info=1;

  if( RVECTOR_COPY( thermostat_p->momenta, thermostat.momenta ) )                 info=1;

  if( RVECTOR_COPY( thermostat_p->forces, thermostat.forces ) )                   info=1;


  return info;

}

//------------------------------------------

/* COMPARE [ a1 = a2 ] */

int PolyCEID_thermostat_compare( const thermostat thermostat1, const thermostat thermostat2 ){

  /* dummies */
  int info=0;


  if( RVECTOR_COMPARE( thermostat1.positions, thermostat2.positions ) ){ 

    fprintf( stderr, "ERROR: positions are not equals\n" );
    fflush( stderr );

    info=1;

  }

  if( RVECTOR_COMPARE( thermostat1.momenta, thermostat2.momenta ) ){ 

    fprintf( stderr, "ERROR: momenta are not equals\n" );
    fflush( stderr );

    info=1;

  }

  if( RVECTOR_COMPARE( thermostat1.forces, thermostat2.forces ) ){ 

    fprintf( stderr, "ERROR: forces are not equals\n" );
    fflush( stderr );

    info=1;

  }


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_thermostat_read( FILE* fp, thermostat_p thermostat_p ){

  /* dummies */
  int info=0;


  if( RVECTOR_READ( fp, thermostat_p->positions ) )       info=1;

  if( RVECTOR_READ( fp, thermostat_p->momenta ) )         info=1;

  if( RVECTOR_READ( fp, thermostat_p->forces ) )          info=1;


  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_thermostat_print( FILE* fp, const thermostat thermostat ){

  /* dummies */
  int info=0;


  /* positions */
  if( RVECTOR_PRINT( fp, thermostat.positions ) ) info=1;

  /* momenta */
  if( RVECTOR_PRINT( fp, thermostat.momenta ) ) info=1;

  /* forces */
  if( RVECTOR_PRINT( fp, thermostat.forces ) ) info=1;

  fflush( fp );


  return info;

}

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_thermostat_verbose_print( FILE* fp, const thermostat thermostat ){

  /* dummies */
  int info=0;


  /* positions */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# positions:\n");

  if( RVECTOR_PRINT_PLUS( fp, thermostat.positions ) ) info=1;

  /* momenta */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# momenta:\n");

  if( RVECTOR_PRINT_PLUS( fp, thermostat.momenta ) ) info=1;

  /* forces */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# forces:\n");

  if( RVECTOR_PRINT_PLUS( fp, thermostat.forces ) ) info=1;

  fflush( fp);


  return info;

}

//------------------------------------------
