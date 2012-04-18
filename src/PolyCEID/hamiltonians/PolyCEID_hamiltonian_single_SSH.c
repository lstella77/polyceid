
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
#include "PolyCEID_hamiltonian_single_SSH.h"


/*********************
  FUNCTIONS & MACROS
*********************/

/* H_matrix_single_SSH update */

int H_matrix_single_SSH_update( const constants constants, state_p state_p, config_p config_p, matrix_p matrix_p ){

  /* constants */
  int        N_electrons;
  int        N_nodes;
  int*       neighbour_index;
  int*       neighbour_list;
  int*       neighbour_link;
  double**   par_node;
  double**   par_link;
  double*    par_extra;
  /* state */
  rvector_p  distances_p;
  /* dummies */
  double     spacing;
  int        index;
  int        i, j;
  double     distance;
  double     dummy;
  double     E_harmonic;
  int        info=0;


  N_electrons     =  constants.N_electrons;
  N_nodes         =  constants.hamiltonian.N_nodes;
  neighbour_index =  constants.hamiltonian.neighbour_index;
  neighbour_list  =  constants.hamiltonian.neighbour_list;
  neighbour_link  =  constants.hamiltonian.neighbour_link;
  par_node        =  constants.hamiltonian.par_node;
  par_link        =  constants.hamiltonian.par_link;
  par_extra       =  constants.hamiltonian.par_extra;
  distances_p     = &(config_p->atoms.distances);

  if( constants.flag_periodic_boundary_condition ){

    spacing = constants.cell_dim.rvector[ 0 ] /constants.N_atoms;

  }
  else{

    spacing = par_extra[ 0 ];         

  }


  if( !info ){

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1;

    E_harmonic = 0.0e0;

    for( i=0; i<N_nodes; i++ ){
      
      // diagonal part
      index = ELECTRON_SINGLE_INDEX( i, i );

      matrix_p->matrix[ index ].z[0] = par_node[ i ][ 0 ];

      // off-diagonal part  
      for( j=neighbour_index[ i ]; j<neighbour_index[ i+1 ]; j++ ){

	/* electronic term update */
	distance = distances_p->rvector[ ATOM_INDEX( i, neighbour_list[ j ] ) ];

	// fprintf( stdout, "distance = %le, a_spacing = %le\n", distance, a_spacing );

        dummy = par_link[ neighbour_link[ j ] ][ 1 ] *( distance -spacing ) +par_link[ neighbour_link[ j ] ][ 2 ];  

	if( dummy < 1.0e0 ){

	  index = ELECTRON_SINGLE_INDEX( i, neighbour_list[ j ] );

	  matrix_p->matrix[ index ].z[0] = -par_link[ neighbour_link[ j ] ][ 0 ] *( 1.0e0 -dummy );

	}
	else{

	  fprintf( stdout, "[ H_matrix ] WARNING: jump at time %le\n", config_p->time );

	}

	/* harmonic term update */
	E_harmonic += par_link[ neighbour_link[ j ] ][ 3 ] *( distance -spacing ) *( distance -spacing );

      }
      
    }

    /* Classical part */
    for( i=0; i<N_nodes; i++ ){

      matrix_p->matrix[ ELECTRON_SINGLE_INDEX( i, i ) ].z[0] = 0.25e0 *E_harmonic / (double) N_electrons; //WARNING: diagonal only 
      
    }

  }

  /*
    fprintf( stdout, "# H_matrix\n" );
    MATRIX_PRINT_PLUS( stdout, *matrix_p );
    fprintf( stdout, "\n" );
  */


  return info;

}

//------------------------------------------

/* F_matrix_single_SSH update */

int F_matrix_single_SSH_update( const constants constants, state_p state_p, config_p config_p, int r, matrix_p matrix_p ){

  /* constants */
  int        N_electrons;
  int        N_nodes;
  int*       neighbour_index;
  int*       neighbour_list;
  int*       neighbour_link;
  double**   par_link;
  double*    par_extra;
  /* state */
  rvector_p  positions_p;
  rvector_p  distances_p;
  /* dummies */
  double     spacing;
  int        index;
  int        j;
  double     distance;
  double     distance_d;
  double     dummy;
  double     dummy_d;
  double     E_harmonic;
  int        info=0;


  N_electrons     =  constants.N_electrons;
  N_nodes         =  constants.hamiltonian.N_nodes;
  neighbour_index =  constants.hamiltonian.neighbour_index;
  neighbour_list  =  constants.hamiltonian.neighbour_list;
  neighbour_link  =  constants.hamiltonian.neighbour_link;
  par_link        =  constants.hamiltonian.par_link;
  par_extra       =  constants.hamiltonian.par_extra;
  positions_p     = &(config_p->atoms.positions);
  distances_p     = &(config_p->atoms.distances);

  if( constants.flag_periodic_boundary_condition ){

    spacing = constants.cell_dim.rvector[ 0 ] /constants.N_atoms;

  }
  else{

    spacing = par_extra[ 0 ];         

  }


  if( !info ){

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1;

    E_harmonic = 0.0e0;

    // there is NO electronic diagonal part

    // off-diagonal part  
    for( j=neighbour_index[ r ]; j<neighbour_index[ r+1 ]; j++ ){

      /* electronic term update */
      distance = distances_p->rvector[ ATOM_INDEX( r, neighbour_list[ j ] ) ];
        
      if( distance > EPS ){

        // BUGFIX: it currently works just in 1D!
        if( constants.flag_periodic_boundary_condition ){

          distance_d = PBC_DIFF( positions_p->rvector[ r ] -positions_p->rvector[ neighbour_list[ j ] ], constants.cell_dim.rvector[ 0 ] ) /distance;

        }
        else{

          distance_d = ( positions_p->rvector[ r ] -positions_p->rvector[ neighbour_list[ j ] ] ) /distance;

        }

      }
      else{

        distance_d = 1.0;

      }	
      // fprintf( stdout, "distance   = %le, a_spacing = %le\n", distance, a_spacing );
      // fprintf( stdout, "distance_d = %le, a_spacing = %le\n", distance_d, a_spacing );

      dummy = par_link[ neighbour_link[ j ] ][ 1 ] *( distance -spacing ) +par_link[ neighbour_link[ j ] ][ 2 ];  

      dummy_d = par_link[ neighbour_link[ j ] ][ 1 ] *distance_d;  

      if( dummy < 1.0e0 ){

        index = ELECTRON_SINGLE_INDEX( r, neighbour_list[ j ] );

        matrix_p->matrix[ index ].z[0] = -par_link[ neighbour_link[ j ] ][ 0 ] *dummy_d;

        index = ELECTRON_SINGLE_INDEX( neighbour_list[ j ], r );

        matrix_p->matrix[ index ].z[0] = -par_link[ neighbour_link[ j ] ][ 0 ] *dummy_d;

      }
      else{

        fprintf( stdout, "[ F_matrix ] WARNING: jump at time %le\n", config_p->time );

      }

      /* harmonic term update */
      E_harmonic -= par_link[ neighbour_link[ j ] ][ 3 ] *( distance -spacing ) *distance_d;

    }

    /* Classical part */
    for( j=0; j<N_nodes; j++ ){

      matrix_p->matrix[ ELECTRON_SINGLE_INDEX( j, j ) ].z[0] = E_harmonic/ (double) N_electrons; //WARNING: diagonal only 
      
    }

  }
  
  /*
    fprintf( stdout, "# F_matrix[%d]\n", r );
    MATRIX_PRINT_PLUS( stdout, *matrix_p );
    fprintf( stdout, "\n" );
  */


  return info;

}

//------------------------------------------

/* K_matrix_single_SSH update */

int K_matrix_single_SSH_update( const constants constants, state_p state_p, config_p config_p, int r, int s, matrix_p matrix_p ){

  /* constants */
  int        N_electrons;
  int        N_nodes;
  int*       neighbour_index;
  int*       neighbour_list;
  int*       neighbour_link;
  double**   par_link;
  double*    par_extra;
  /* state */
  rvector_p  positions_p;
  rvector_p  distances_p;
  /* dummies */
  double     spacing;
  int        index;
  int        j;
  double     distance;
  double     distance_d;
  double     sign;
  double     dummy;
  double     dummy_d;
  double     E_harmonic;
  int        info=0;


  N_electrons     =  constants.N_electrons;
  N_nodes         =  constants.hamiltonian.N_nodes;
  neighbour_index =  constants.hamiltonian.neighbour_index;
  neighbour_list  =  constants.hamiltonian.neighbour_list;
  neighbour_link  =  constants.hamiltonian.neighbour_link;
  par_link        =  constants.hamiltonian.par_link;
  par_extra       =  constants.hamiltonian.par_extra;
  positions_p     = &(config_p->atoms.positions);
  distances_p     = &(config_p->atoms.distances);

  if( constants.flag_periodic_boundary_condition ){

    spacing = constants.cell_dim.rvector[ 0 ] /constants.N_atoms;

  }
  else{

    spacing = par_extra[ 0 ];         

  }


  if( !info ){

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1;

    E_harmonic = 0.0e0;

    // there is NO electronic diagonal part

    // off-diagonal part  
    for( j=neighbour_index[ r ]; j<neighbour_index[ r+1 ]; j++ ){

      if( s == neighbour_list[ j ] || s == r ){

        /* electronic term update */
        distance = distances_p->rvector[ ATOM_INDEX( r, neighbour_list[ j ] ) ];
        
        if( distance > EPS ){

          // BUGFIX: it currently works just in 1D!
          if( constants.flag_periodic_boundary_condition ){

            distance_d = PBC_DIFF( positions_p->rvector[ r ] -positions_p->rvector[ neighbour_list[ j ] ], constants.cell_dim.rvector[ 0 ] ) /distance;

          }
          else{

            distance_d = ( positions_p->rvector[ r ] -positions_p->rvector[ neighbour_list[ j ] ] ) /distance;

          }

          distance_d = ( distance_d *distance_d - 1.0 ) /distance;

        }
        else{

          distance_d = 0.0;

        }	
        // fprintf( stdout, "distance   = %le, a_spacing = %le\n", distance, a_spacing );
        // fprintf( stdout, "distance_d = %le, a_spacing = %le\n", distance_d, a_spacing );

        dummy = par_link[ neighbour_link[ j ] ][ 1 ] *( distance -spacing ) +par_link[ neighbour_link[ j ] ][ 2 ];  

        dummy_d = par_link[ neighbour_link[ j ] ][ 1 ] *distance_d;  

        sign = -1.0;
        if ( s == r ){

          sign = -sign;

	}

        if( dummy < 1.0e0 ){

          index = ELECTRON_SINGLE_INDEX( r, neighbour_list[ j ] );

          matrix_p->matrix[ index ].z[0] = -sign *par_link[ neighbour_link[ j ] ][ 0 ] *dummy_d;

          index = ELECTRON_SINGLE_INDEX( neighbour_list[ j ], r );

          matrix_p->matrix[ index ].z[0] = -sign *par_link[ neighbour_link[ j ] ][ 0 ] *dummy_d;

        }
        else{

          fprintf( stdout, "[ K_matrix ] WARNING: jump at time %le\n", config_p->time );

        }

        /* harmonic term update */
        E_harmonic += sign *par_link[ neighbour_link[ j ] ][ 3 ] *( 1.0 -spacing *dummy_d );

      }

    }

    /* Classical part */
    for( j=0; j<N_nodes; j++ ){

      matrix_p->matrix[ ELECTRON_SINGLE_INDEX( j, j ) ].z[0] = E_harmonic/ (double) N_electrons; //WARNING: diagonal only 
      
    }

  }

  /*
    fprintf( stdout, "# K_matrix[%d,%d]\n", r, s );
    MATRIX_PRINT_PLUS( stdout, *matrix_p );
    fprintf( stdout, "\n" );
  */


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

int Hamiltonian_single_SSH_parameters_check( const constants constants ){

  /*constants */
  int     N_atoms;
  int     N_levels_single;
  /* dummies */
  int info=0;


  N_atoms          = constants.N_atoms;
  N_levels_single  = constants.N_levels_single;


#ifdef __DEBUG__

  fprintf( stdout, "----------\n" );
  fprintf( stdout, "HAMILTONIAN_SINGLE: using SSH\n" );

#endif /* __DEBUG__ */


  if( N_levels_single != N_atoms ){

    fprintf( stderr, "ERROR: this Hamiltonian_single works only for N_levels_single = N_atoms.\n" );
    fflush( stderr );

    info=1;

  }

  if( !constants.flag_periodic_boundary_condition && constants.hamiltonian.N_par_extra != 1 ){

   fprintf( stderr, "ERROR: this Hamiltonian_single wants an extra parameter (spacing) if PBC are not defined.");
   fflush( stderr );

   info=1;

  }        
  

  return info;

}

//------------------------------------------

/* H_matrix_single_SSH_check */

int H_matrix_single_SSH_check( const constants constants, rvector_p positions_p, matrix_p matrix_p ){

  /* constants */
  int        N_electrons;
  int        N_nodes;
  int*       neighbour_index;
  int*       neighbour_list;
  int*       neighbour_link;
  double**   par_node;
  double**   par_link;
  double*    par_extra;
  /* dummies */
  double     spacing;
  int        index;
  int        i, j;
  double     distance;
  double     dummy;
  double     E_harmonic;
  int        info=0;


  N_electrons     =  constants.N_electrons;
  N_nodes         =  constants.hamiltonian.N_nodes;
  neighbour_index =  constants.hamiltonian.neighbour_index;
  neighbour_list  =  constants.hamiltonian.neighbour_list;
  neighbour_link  =  constants.hamiltonian.neighbour_link;
  par_node        =  constants.hamiltonian.par_node;
  par_link        =  constants.hamiltonian.par_link;
  par_extra       =  constants.hamiltonian.par_extra;

  if( constants.flag_periodic_boundary_condition ){

    spacing = constants.cell_dim.rvector[ 0 ] /constants.N_atoms;

  }
  else{

    spacing = par_extra[ 0 ];         

  }


  if( !info ){

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1;

    E_harmonic = 0.0e0;

    for( i=0; i<N_nodes; i++ ){
      
      // diagonal part
      index = ELECTRON_SINGLE_INDEX( i, i );

      matrix_p->matrix[ index ].z[0] = par_node[ i ][ 0 ];

      // off-diagonal part  
      for( j=neighbour_index[ i ]; j<neighbour_index[ i+1 ]; j++ ){

	/* electronic term update */
        distance = DISTANCE_LOC( constants, i, neighbour_list[ j ], *positions_p  );

	// fprintf( stdout, "distance = %le, a_spacing = %le\n", distance, a_spacing );

        dummy = par_link[ neighbour_link[ j ] ][ 1 ] *( distance -spacing ) +par_link[ neighbour_link[ j ] ][ 2 ];  

	if( dummy < 1.0e0 ){

	  index = ELECTRON_SINGLE_INDEX( i, neighbour_list[ j ] );

	  matrix_p->matrix[ index ].z[0] = -par_link[ neighbour_link[ j ] ][ 0 ] *( 1.0e0 -dummy );

	}

	/* harmonic term update */
	E_harmonic += par_link[ neighbour_link[ j ] ][ 3 ] *( distance -spacing ) *( distance -spacing );

      }
      
    }

    /* Classical part */
    for( i=0; i<N_nodes; i++ ){

      matrix_p->matrix[ ELECTRON_SINGLE_INDEX( i, i ) ].z[0] = 0.25e0 *E_harmonic /(double) N_electrons; //WARNING: diagonal only 
      
    }

  }

  /* 
    fprintf( stdout, "# positions [check]\n" );
    RVECTOR_PRINT_PLUS( stdout, *positions_p );
    fprintf( stdout, "# H_matrix\n" );
    MATRIX_PRINT_PLUS( stdout, *matrix_p );
    fprintf( stdout, "\n" );
  */


  return info;

}

//------------------------------------------

double distance_loc( const constants constants, int atom_i, int atom_j, rvector_p positions_p ){

  /* constants*/
  int    sdim;
  /* dummies */
  int    comp;
  double distance;
  double dummy;
  int    info=0;


  sdim   =  constants.spacial_dimension; 

  if( positions_p->rvector_dim != constants.N_coor ) info=1;


  if( !info ){

    distance = 0.0e0;

    for( comp=0; comp<sdim; comp++ ){

      dummy = ( positions_p->rvector[ comp +atom_i *sdim ] -positions_p->rvector[ comp +atom_j *sdim ] );

      if( constants.flag_periodic_boundary_condition ){

        dummy = PBC_DIFF( dummy, constants.cell_dim.rvector[ comp ] );

      }

      dummy = dummy *dummy;

      distance += dummy;

    }

    distance = sqrt( distance );

  }

  if( info ){

    fprintf( stderr, "ERROR: invalid distance generated\n" );
    fflush( stderr );

    distance = HUGE_VAL;

  }


  return distance;

}

//------------------------------------------

/* classical_dipole_single_SSH_update */

int classical_dipole_single_SSH_update( const constants constants, state_p state_p, config_p config_p, int comp, matrix_p matrix_p ){

  /* constants */
  int        N_atoms;
  int        sdim;
  /* state */
  rvector_p  positions_p;
  rvector_p  centre_of_mass_p;
  /* dummies */
  int        index;
  int        i;
  int        info=0;


  N_atoms           =  constants.N_atoms;
  sdim              =  constants.spacial_dimension;
  positions_p       = &(config_p->atoms.positions);
  centre_of_mass_p  = &(config_p->atoms.centre_of_mass);


  if( !info ){

    /* set matrix to zero */
    if( MATRIX_ZERO( *matrix_p ) ) info=1; // Important

    for( i=0; i<N_atoms; i++){
      
      index = ELECTRON_SINGLE_INDEX( i, i );
 
      matrix_p->matrix[ index ].z[0] = positions_p->rvector[ comp +i*sdim ] -( centre_of_mass_p->rvector[ comp ] );

    } /* i loop */
    
  }

  /* 
    fprintf( stdout, "# positions [dipole]\n" );
    RVECTOR_PRINT_PLUS( stdout, *positions_p );
    fprintf( stdout, "# centre of mass [dipole]\n" );
    RVECTOR_PRINT_PLUS( stdout, *centre_of_mass_p );
    fprintf( stdout, "# dipole_single\n" );
    MATRIX_PRINT_PLUS( stdout, *matrix_p );
    fprintf( stdout, "\n" );
  */


  return info;

}

//------------------------------------------
