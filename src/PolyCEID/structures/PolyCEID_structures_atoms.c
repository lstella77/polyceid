
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
#include "PolyCEID_structures_atoms.h"


/* ALLOCATE */

int PolyCEID_atoms_allocate( int np, int* pp, atoms_p atoms_p ){

  /* atoms */
  int  N_atoms;
  int  sdim;
  int  N_coor;
  /* dummies */
  int  i;
  int  info=0;

  if( np != NP_ATOMS ) info=1;

  N_atoms = pp[ 0 ];
  atoms_p->N_atoms = N_atoms;

  sdim = pp [ 1 ];
  atoms_p->spacial_dimension = sdim;

  N_coor  = N_atoms *sdim;

  // allocating atom names
  atoms_p->names = ( char** ) calloc( (size_t) N_atoms, sizeof( char* ) );

  for( i=0; i<N_atoms; i++ ){

    atoms_p->names[ i ] = ( char* ) calloc( (size_t) ATOM_NAME_LENGTH, sizeof( char ) ); 

    strcpy( atoms_p->names[ i ], "No_name" );

  }        

  if( RVECTOR_ALLOCATE( N_coor, atoms_p->masses ) )              info=1;

  if( RVECTOR_ALLOCATE( N_coor, atoms_p->masses_aux ) )          info=1;

  if( RVECTOR_ALLOCATE( N_coor, atoms_p->positions ) )           info=1;

  if( RVECTOR_ALLOCATE( N_coor, atoms_p->momenta ) )             info=1;

  if( RVECTOR_ALLOCATE( N_coor, atoms_p->forces ) )              info=1;

  if( RVECTOR_ALLOCATE( N_coor, atoms_p->forces_cons ) )         info=1;

  if( RVECTOR_ALLOCATE( N_atoms *N_atoms, atoms_p->distances ) ) info=1;

  if( RVECTOR_ALLOCATE( sdim, atoms_p->centre_of_mass ) )        info=1;

  if( IVECTOR_ALLOCATE( N_coor, atoms_p->mask ) )                info=1;


  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_atoms_free( atoms_p atoms_p ){

  /* dummies */
  int N_atoms;
  int i;
  int info=0;


  N_atoms = atoms_p->N_atoms;

  if( IVECTOR_FREE( atoms_p->mask ) )            info=1;

  if( RVECTOR_FREE( atoms_p->centre_of_mass ) )  info=1;

  if( RVECTOR_FREE( atoms_p->distances ) )       info=1;

  if( RVECTOR_FREE( atoms_p->forces_cons ) )     info=1;

  if( RVECTOR_FREE( atoms_p->forces ) )          info=1;

  if( RVECTOR_FREE( atoms_p->momenta ) )         info=1;

  if( RVECTOR_FREE( atoms_p->positions ) )       info=1;

  if( RVECTOR_FREE( atoms_p->masses_aux ) )      info=1;

  if( RVECTOR_FREE( atoms_p->masses ) )          info=1;

  for( i=0; i<N_atoms; i++ ){

    free( atoms_p->names[ i ] );

  }
                                    
  free( atoms_p->names );


  return info;

}

//------------------------------------------

/* COPY [ a1 = a2 ] */

int PolyCEID_atoms_copy( atoms_p atoms_p, const atoms atoms ){

  /* dummies */
  int N_atoms;
  int i;
  int info=0;


  N_atoms = atoms_p->N_atoms;

  for( i=0; i<N_atoms; i++ ){

    strcpy( atoms_p->names[ i ], atoms.names[ i ] );

  }        

  if( RVECTOR_COPY( atoms_p->positions, atoms.positions ) )           info=1;

  if( RVECTOR_COPY( atoms_p->momenta, atoms.momenta ) )               info=1;

  if( RVECTOR_COPY( atoms_p->forces, atoms.forces ) )                 info=1;

  if( RVECTOR_COPY( atoms_p->forces_cons, atoms.forces_cons ) )       info=1;

  if( RVECTOR_COPY( atoms_p->masses, atoms.masses ) )                 info=1;

  if( RVECTOR_COPY( atoms_p->masses_aux, atoms.masses_aux ) )         info=1;

  atoms_p->mass_tot = atoms.mass_tot;

  atoms_p->mass_ave = atoms.mass_ave;

  if( RVECTOR_COPY( atoms_p->distances, atoms.distances ) )           info=1;

  if( RVECTOR_COPY( atoms_p->centre_of_mass, atoms.centre_of_mass ) ) info=1;

  if( IVECTOR_COPY( atoms_p->mask, atoms.mask ) )                     info=1;


  return info;

}

//------------------------------------------

/* COMPARE [ a1 = a2 ] */

int PolyCEID_atoms_compare( const atoms atoms1, const atoms atoms2 ){

  /* dummies */
  int N_atoms1;
  int N_atoms2;
  int i;
  int info=0;


  N_atoms1 = atoms1.N_atoms;

  N_atoms2 = atoms2.N_atoms;


  if( N_atoms1 != N_atoms2 ){

    fprintf( stderr, "ERROR: N_atoms are not equals\n" );
    fflush( stderr );

    info=1;

  }


  for( i=0; i<N_atoms1; i++ ){

    if( strcmp( atoms1.names[ i ], atoms2.names[ i ] ) ){

      fprintf( stderr, "ERROR: names are not equals\n" );
      fflush( stderr );

      break;

      info=1;

    }  

  }


  if( RVECTOR_COMPARE( atoms1.positions, atoms2.positions ) ){ 

    fprintf( stderr, "ERROR: positions are not equals\n" );
    fflush( stderr );

    info=1;

  }


  if( RVECTOR_COMPARE( atoms1.momenta, atoms2.momenta ) ){ 

    fprintf( stderr, "ERROR: momenta are not equals\n" );
    fflush( stderr );

    info=1;

  }

  if( RVECTOR_COMPARE( atoms1.forces, atoms2.forces ) ){ 

    fprintf( stderr, "ERROR: forces are not equals\n" );
    fflush( stderr );

    info=1;

  }

  if( RVECTOR_COMPARE( atoms1.masses, atoms2.masses ) ){ 

    fprintf( stderr, "ERROR: masses are not equals\n" );
    fflush( stderr );

    info=1;

  }

  if( RVECTOR_COMPARE( atoms1.masses_aux, atoms2.masses_aux ) ){ 

    fprintf( stderr, "ERROR: masses_aux are not equals\n" );
    fflush( stderr );

    info=1;

  }

  if( fabs( atoms1.mass_tot -atoms2.mass_tot ) > EPS ){ 

    fprintf( stderr, "ERROR: mass_tot's are not equals\n" );
    fflush( stderr );

    info=1;

  }

  if( fabs( atoms1.mass_ave -atoms2.mass_ave ) > EPS ){ 

    fprintf( stderr, "ERROR: mass_ave's are not equals\n" );
    fflush( stderr );

    info=1;

  }

  if( RVECTOR_COMPARE( atoms1.distances, atoms2.distances ) ){ 

    fprintf( stderr, "ERROR: distances are not equals\n" );
    fflush( stderr );

    info=1;

  }


  if( RVECTOR_COMPARE( atoms1.centre_of_mass, atoms2.centre_of_mass ) ){ 

    fprintf( stderr, "ERROR: centres of mass are not equals\n" );
    fflush( stderr );

    info=1;

  }


  if( IVECTOR_COMPARE( atoms1.mask, atoms2.mask ) ){ 

    fprintf( stderr, "ERROR: masks are not equals\n" );
    fflush( stderr );

    info=1;

  }


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_atoms_read( FILE* fp, atoms_p atoms_p ){

  /* dummies */
  int N_atoms;
  int i;
  int info=0;


  N_atoms = atoms_p->N_atoms;


  for( i=0; i<N_atoms; i++ ){

    if( fscanf( fp, "%s", atoms_p->names[ i ] ) ) info=1;

  }

  if( RVECTOR_READ( fp, atoms_p->masses ) )          info=1;

  if( RVECTOR_READ( fp, atoms_p->masses_aux ) )      info=1;

  if( fscanf( fp, "%le", &atoms_p->mass_tot ) < 1 )  info=1;

  if( fscanf( fp, "%le", &atoms_p->mass_ave ) < 1 )  info=1;

  if( RVECTOR_READ( fp, atoms_p->positions ) )       info=1;

  if( RVECTOR_READ( fp, atoms_p->momenta ) )         info=1;

  if( RVECTOR_READ( fp, atoms_p->forces ) )          info=1;

  if( RVECTOR_READ( fp, atoms_p->distances ) )       info=1;

  if( RVECTOR_READ( fp, atoms_p->centre_of_mass ) )  info=1;

  if( IVECTOR_READ( fp, atoms_p->mask ) )            info=1;


  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_atoms_print( FILE* fp, const atoms atoms ){

  /* dummies */
  int N_atoms;
  int i;
  int info=0;


  N_atoms = atoms.N_atoms;

  for( i=0; i<N_atoms; i++ ){

    fprintf( fp, "%s  ", atoms.names[ i ] );

  }
  fprintf( fp, "\n");

  /* masses */
  if( RVECTOR_PRINT( fp, atoms.masses ) )         info=1;

  /* masses_aux */
  if( RVECTOR_PRINT( fp, atoms.masses_aux ) )    info=1;

  /* masses_tot */
  if( fprintf( fp, "# %12.5le\n", atoms.mass_tot ) ) info=1;

  /* masses_ave */
  if( fprintf( fp, "# %12.5le\n", atoms.mass_ave ) ) info=1;

  /* positions */
  if( RVECTOR_PRINT( fp, atoms.positions ) )      info=1;

  /* momenta */
  if( RVECTOR_PRINT( fp, atoms.momenta ) )        info=1;

  /* forces */
  if( RVECTOR_PRINT( fp, atoms.forces ) )         info=1;

  /* forces_cons */
  if( RVECTOR_PRINT( fp, atoms.forces_cons ) )    info=1;

  /* distances */
  if( RVECTOR_PRINT( fp, atoms.distances ) )      info=1;

  /* distances */
  if( RVECTOR_PRINT( fp, atoms.centre_of_mass ) ) info=1;

  /* mask */
  if( IVECTOR_PRINT( fp, atoms.mask ) )           info=1;

  fflush( fp );


  return info;

}

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_atoms_verbose_print( FILE* fp, const atoms atoms ){

  /* dummies */
  int N_atoms;
  int i;
  int info=0;


  N_atoms = atoms.N_atoms;

  /* names */
  for( i=0; i<N_atoms; i++ ){

    fprintf( fp, "# %s\n", atoms.names[ i ] );

  }

  /* masses */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# masses:\n");

  if( RVECTOR_PRINT_PLUS( fp, atoms.masses ) ) info=1;

  /* masses_aux */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# masses_aux:\n");

  if( RVECTOR_PRINT_PLUS( fp, atoms.masses_aux ) ) info=1;

  /* masses_tot */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# masses_tot:\n");

  if( fprintf( fp, "# %12.5le\n", atoms.mass_tot ) < 1 ) info=1;

  /* masses_ave */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# masses_ave:\n");

  if( fprintf( fp, "# %12.5le\n", atoms.mass_ave ) < 1 ) info=1;

  /* positions */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# positions:\n");

  if( RVECTOR_PRINT_PLUS( fp, atoms.positions ) ) info=1;

  /* momenta */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# momenta:\n");

  if( RVECTOR_PRINT_PLUS( fp, atoms.momenta ) ) info=1;

  /* forces */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# forces:\n");

  if( RVECTOR_PRINT_PLUS( fp, atoms.forces ) ) info=1;

  /* forces_cons */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# forces_cons:\n");

  if( RVECTOR_PRINT_PLUS( fp, atoms.forces_cons ) ) info=1;

  /* distances */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# distances:\n");

  if( RVECTOR_PRINT_PLUS( fp, atoms.distances ) ) info=1;

  /* centre_of_mass */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# centre of mass:\n");

  if( RVECTOR_PRINT_PLUS( fp, atoms.centre_of_mass ) ) info=1;

  /* mask */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# mask:\n");

  if( IVECTOR_PRINT_PLUS( fp, atoms.mask ) ) info=1;


  fflush( fp);


  return info;

}

//------------------------------------------
