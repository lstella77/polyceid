
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
#include "PolyCEID_structures_electrons.h"


/* ALLOCATE */

int PolyCEID_electrons_allocate( int np, int* pp, electrons_p electrons_p ){

  /* electrons */
  int  N_coor; // the original value!
  int  N_levels_single;
  int  N_levels_many;
  int  max_rho_index;
  /* dummies */
  int  dim;
  int  info=0;


  if( np != NP_ELECTRONS ) info=1;


  N_levels_single = pp[0];
  N_levels_many   = pp[1];
  N_coor          = pp[2];
  max_rho_index   = pp[3];


  /* H_matrix */
  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, electrons_p->H_matrix ) ) info=1;

  /* F_matrix */
  dim = N_coor;

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->F_matrix ) ) info=1;

  /* delta_F_matrix */

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->delta_F_matrix ) ) info=1;

  /* K_matrix */
  dim = N_coor *N_coor; //WARNING: can be exceeding large!

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->K_matrix ) ) info=1;

  /* dipole_many */
  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, electrons_p->dipole_many ) ) info=1;

  /* Ehrenfest_frame */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, electrons_p->Ehrenfest_frame ) ) info=1;

  /* rho */
  dim = max_rho_index+1;  //WARNING: 1) can be exceeding large! 2) one black matrix at the end

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->rho ) ) info=1;

  /* mu00 */
  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, electrons_p->mu00 ) ) info=1;

  /* mu01 */
  dim = N_coor;

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->mu01 ) ) info=1;

  /* mu10 */
  dim = N_coor;

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->mu10 ) ) info=1;

  /* mu02 */
  dim = N_coor *N_coor; //WARNING: can be exceeding large!

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->mu02 ) ) info=1;

  /* mu11 */
  dim = N_coor *N_coor; //WARNING: can be exceeding large!

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->mu11 ) ) info=1;

  /* mu20 */
  dim = N_coor *N_coor; //WARNING: can be exceeding large!

  if( MATRIX_ARRAY_ALLOCATE( dim, N_levels_many, N_levels_many, electrons_p->mu20 ) ) info=1;

  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_electrons_free( electrons_p electrons_p ){

  /* dummies */
  int  info=0;


  /* mu20 */
  if( MATRIX_ARRAY_FREE( electrons_p->mu20 ) ) info=1;

  /* mu11 */
  if( MATRIX_ARRAY_FREE( electrons_p->mu11 ) ) info=1;

  /* mu02 */
  if( MATRIX_ARRAY_FREE( electrons_p->mu02 ) ) info=1;

  /* mu10 */
  if( MATRIX_ARRAY_FREE( electrons_p->mu10 ) ) info=1;

  /* mu01 */
  if( MATRIX_ARRAY_FREE( electrons_p->mu01 ) ) info=1;

  /* mu00 */
  if( MATRIX_FREE( electrons_p->mu00 ) ) info=1;

  /* rho */
  if( MATRIX_ARRAY_FREE( electrons_p->rho ) ) info=1;

  /* Ehrenfest_frame */
  if( MATRIX_FREE( electrons_p->Ehrenfest_frame ) ) info=1;

  /* dipole_many */
  if( MATRIX_FREE( electrons_p->dipole_many ) ) info=1;

  /* K_matrix */
  if( MATRIX_ARRAY_FREE( electrons_p->K_matrix ) ) info=1;

  /* delta_F_matrix */
  if( MATRIX_ARRAY_FREE( electrons_p->delta_F_matrix ) ) info=1;

  /* F_matrix */
  if( MATRIX_ARRAY_FREE( electrons_p->F_matrix ) ) info=1;

  /* H_matrix */
  if( MATRIX_FREE( electrons_p->H_matrix ) ) info=1;


  return info;

}

//------------------------------------------

/* COPY [ e1 = e2 ] */

int PolyCEID_electrons_copy( electrons_p electrons_p, const electrons electrons ){

  /* dummies */
  int info=0;


  /* H_matrix */
  if( MATRIX_COPY( electrons_p->H_matrix, electrons.H_matrix ) ) info=1;

  /* F_matrix */
  if( MATRIX_ARRAY_COPY( electrons_p->F_matrix, electrons.F_matrix ) ) info=1;

  /* delta_F_matrix */
  if( MATRIX_ARRAY_COPY( electrons_p->delta_F_matrix, electrons.delta_F_matrix ) ) info=1;

  /* K_matrix */
  if( MATRIX_ARRAY_COPY( electrons_p->K_matrix, electrons.K_matrix ) ) info=1;

  /* dipole_many */
  if( MATRIX_COPY( electrons_p->dipole_many, electrons.dipole_many ) ) info=1;

  /* Ehrenfest_frame */
  if( MATRIX_COPY( electrons_p->Ehrenfest_frame, electrons.Ehrenfest_frame ) ) info=1;

  /* rho */
  if( MATRIX_ARRAY_COPY( electrons_p->rho, electrons.rho ) ) info=1;

  /* mu00 */
  if( MATRIX_COPY( electrons_p->mu00, electrons.mu00 ) ) info=1;

  /* mu01 */
  if( MATRIX_ARRAY_COPY( electrons_p->mu01, electrons.mu01 ) ) info=1;

  /* mu10 */
  if( MATRIX_ARRAY_COPY( electrons_p->mu10, electrons.mu10 ) ) info=1;

  /* mu02 */
  if( MATRIX_ARRAY_COPY( electrons_p->mu02, electrons.mu02 ) ) info=1;

  /* mu11 */
  if( MATRIX_ARRAY_COPY( electrons_p->mu11, electrons.mu11 ) ) info=1;

  /* mu20 */
  if( MATRIX_ARRAY_COPY( electrons_p->mu20, electrons.mu20 ) ) info=1;


  return info;

}

//------------------------------------------

/* COMPARE [ e1 = e2 ] */

int PolyCEID_electrons_compare( const electrons electrons1, const electrons electrons2 ){

  /* dummies */
  int info=0;


  /* H_matrix */
  if( MATRIX_COMPARE( electrons1.H_matrix, electrons2.H_matrix ) ){

    fprintf( stderr, "ERROR: H_matrix are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* F_matrix */
  if( MATRIX_ARRAY_COMPARE( electrons1.F_matrix, electrons2.F_matrix ) ){

    fprintf( stderr, "ERROR: F_matrix are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* delta_F_matrix */
  if( MATRIX_ARRAY_COMPARE( electrons1.delta_F_matrix, electrons2.delta_F_matrix ) ){

    fprintf( stderr, "ERROR: delta_F_matrix are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* K_matrix */
  if( MATRIX_ARRAY_COMPARE( electrons1.K_matrix, electrons2.K_matrix ) ){

    fprintf( stderr, "ERROR: K_matrix are not equals\n" );
    fflush( stderr );

    info=1;

  }


  /* dipole_many */
  if( MATRIX_COMPARE( electrons1.dipole_many, electrons2.dipole_many ) ){

    fprintf( stderr, "ERROR: dipole_many are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* Ehrenfest_frame */
  if( MATRIX_COMPARE( electrons1.Ehrenfest_frame, electrons2.Ehrenfest_frame ) ){

    fprintf( stderr, "ERROR: Ehrenfest_frame are not equals\n" );
    fflush( stderr );

    info=1;

  }


  /* rho */
  if( MATRIX_ARRAY_COMPARE( electrons1.rho, electrons2.rho ) ){

    fprintf( stderr, "ERROR: rho are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* mu00 */
  if( MATRIX_COMPARE( electrons1.mu00, electrons2.mu00 ) ){

    fprintf( stderr, "ERROR: mu00 are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* mu01 */
  if( MATRIX_ARRAY_COMPARE( electrons1.mu01, electrons2.mu01 ) ){

    fprintf( stderr, "ERROR: mu01 are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* mu10 */
  if( MATRIX_ARRAY_COMPARE( electrons1.mu10, electrons2.mu10 ) ){

    fprintf( stderr, "ERROR: mu10 are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* mu02 */
  if( MATRIX_ARRAY_COMPARE( electrons1.mu02, electrons2.mu02 ) ){

    fprintf( stderr, "ERROR: mu02 are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* mu11 */
  if( MATRIX_ARRAY_COMPARE( electrons1.mu11, electrons2.mu11 ) ){

    fprintf( stderr, "ERROR: mu11 are not equals\n" );
    fflush( stderr );

    info=1;

  }

  /* mu20 */
  if( MATRIX_ARRAY_COMPARE( electrons1.mu20, electrons2.mu20 ) ){

    fprintf( stderr, "ERROR: mu20 are not equals\n" );
    fflush( stderr );

    info=1;

  }


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_electrons_read( FILE* fp, electrons_p electrons_p ){

  /* dummies */
  int info=0;


  /* H_matrix */
  if( MATRIX_READ( fp, electrons_p->H_matrix ) ) info=1;

  /* F_matrix */
  if( MATRIX_ARRAY_READ( fp, electrons_p->F_matrix ) ) info=1;

  /* delta_F_matrix */
  if( MATRIX_ARRAY_READ( fp, electrons_p->delta_F_matrix ) ) info=1;

  /* K_matrix */
  if( MATRIX_ARRAY_READ( fp, electrons_p->K_matrix ) ) info=1;

  /* dipole_many */
  if( MATRIX_READ( fp, electrons_p->dipole_many ) ) info=1;

  /* Ehrenfest_frame */
  if( MATRIX_READ( fp, electrons_p->Ehrenfest_frame ) ) info=1;

  /* rho */
  if( MATRIX_ARRAY_READ( fp, electrons_p->rho ) ) info=1;

  /* mu00 */
  if( MATRIX_READ( fp, electrons_p->mu00 ) ) info=1;

  /* mu01 */
  if( MATRIX_ARRAY_READ( fp, electrons_p->mu01 ) ) info=1;

  /* mu10 */
  if( MATRIX_ARRAY_READ( fp, electrons_p->mu10 ) ) info=1;

  /* mu02 */
  if( MATRIX_ARRAY_READ( fp, electrons_p->mu02 ) ) info=1;

  /* mu11 */
  if( MATRIX_ARRAY_READ( fp, electrons_p->mu11 ) ) info=1;

  /* mu20 */
  if( MATRIX_ARRAY_READ( fp, electrons_p->mu20 ) ) info=1;


  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_electrons_print( FILE* fp, const electrons electrons ){

  /* dummies */
  int info=0;


  /* H_matrix */
  if( MATRIX_PRINT( fp, electrons.H_matrix ) ) info=1;

  /* F_matrix */
  if( MATRIX_ARRAY_PRINT( fp, electrons.F_matrix ) ) info=1;

  /* delta_F_matrix */
  if( MATRIX_ARRAY_PRINT( fp, electrons.delta_F_matrix ) ) info=1;

  /* K_matrix */
  if( MATRIX_ARRAY_PRINT( fp, electrons.K_matrix ) ) info=1;

  /* dipole_many */
  if( MATRIX_PRINT( fp, electrons.dipole_many ) ) info=1;

  /* Ehrenfest_frame */
  if( MATRIX_PRINT( fp, electrons.Ehrenfest_frame ) ) info=1;

  /* rho */
  if( MATRIX_ARRAY_PRINT( fp, electrons.rho ) ) info=1;

  /* mu00 */
  if( MATRIX_PRINT( fp, electrons.mu00 ) ) info=1;

  /* mu01 */
  if( MATRIX_ARRAY_PRINT( fp, electrons.mu01 ) ) info=1;

  /* mu10 */
  if( MATRIX_ARRAY_PRINT( fp, electrons.mu10 ) ) info=1;

  /* mu02 */
  if( MATRIX_ARRAY_PRINT( fp, electrons.mu02 ) ) info=1;

  /* mu11 */
  if( MATRIX_ARRAY_PRINT( fp, electrons.mu11 ) ) info=1;

  /* mu20 */
  if( MATRIX_ARRAY_PRINT( fp, electrons.mu20 ) ) info=1;


  fflush( fp );


  return info;

}

//------------------------------------------

/* PRINT VERBOSE */

int PolyCEID_electrons_verbose_print( FILE* fp, const electrons electrons ){

  /* dummies */
  int info=0;


  /* H_matrix */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# H_matrix:\n" );

  if( MATRIX_PRINT_PLUS( fp, electrons.H_matrix ) ) info=1;

  /* F_matrix */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# F_matrix:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.F_matrix ) ) info=1;

  /* delta_F_matrix */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# delta_F_matrix:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.delta_F_matrix ) ) info=1;

  /* K_matrix */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# K_matrix:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.K_matrix ) ) info=1;

  /* dipole_many */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# dipole_many:\n" );

  if( MATRIX_PRINT_PLUS( fp, electrons.dipole_many ) ) info=1;

  /* Ehrenfest_frame */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# Ehrenfest_frame:\n" );

  if( MATRIX_PRINT_PLUS( fp, electrons.Ehrenfest_frame ) ) info=1;

#ifdef __DEBUG_PLUS__

  /* rho */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# rho:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.rho ) ) info=1;

#endif /* __DEBUG_PLUS__ */

  /* mu00 */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# mu00:\n" );

  if( MATRIX_PRINT_PLUS( fp, electrons.mu00 ) ) info=1;

#ifdef __DEBUG_PLUS__

  /* mu01 */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# mu01:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.mu01 ) ) info=1;

  /* mu10 */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# mu10:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.mu10 ) ) info=1;

  /* mu02 */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# mu02:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.mu02 ) ) info=1;

  /* mu11 */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# mu11:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.mu11 ) ) info=1;

  /* mu20 */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# mu20:\n" );

  if( MATRIX_ARRAY_PRINT_PLUS( fp, electrons.mu20 ) ) info=1;

#endif /* __DEBUG_PLUS__ */


  fflush( fp );


  return info;

}

//------------------------------------------

