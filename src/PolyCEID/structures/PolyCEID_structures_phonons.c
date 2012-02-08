
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
#include "PolyCEID_structures_phonons.h"


/* ALLOCATE */

int PolyCEID_phonons_allocate( int np, int* pp, phonons_p phonons_p ){

  /* phonons */
  int N_coor;
  int sqrt_max_rho_index;
  /* dummies */
  int info=0;


  if( np != NP_PHONONS ) info=1;


  N_coor              = pp[0];
  sqrt_max_rho_index  = pp[1];


  /* initial_condition_atoms */
  if( VECTOR_ALLOCATE( sqrt_max_rho_index, phonons_p->initial_condition_atoms ) ) info=1;

  /* original_hessian */
  if( MATRIX_ALLOCATE( N_coor, N_coor, phonons_p->original_hessian ) ) info=1;

  /* canonical_transform */
  if( MATRIX_ALLOCATE( N_coor, N_coor, phonons_p->canonical_transform ) ) info=1;

  /* canonical frequencies */
  if( RVECTOR_ALLOCATE( N_coor, phonons_p->canonical_frequencies ) ) info=1;

  /* omega */
  if( RVECTOR_ALLOCATE( N_coor, phonons_p->omega ) ) info=1;

  /* ar */
  if( RVECTOR_ALLOCATE( N_coor, phonons_p->ar ) ) info=1;

  /* ap */
  if( RVECTOR_ALLOCATE( N_coor, phonons_p->ap ) ) info=1;

  /* phononic_dos */
  if( RVECTOR_ALLOCATE( N_coor, phonons_p->phononic_dos ) ) info=1;

  /* symplectic transform R */
  if( MATRIX_ALLOCATE( N_coor, N_coor, phonons_p->symplectic_transform_R ) ) info=1;

  /* conjugate symplectic transform R */
  if( MATRIX_ALLOCATE( N_coor, N_coor, phonons_p->conjugate_symplectic_transform_R ) ) info=1;

  /* symplectic transform P */
  if( MATRIX_ALLOCATE( N_coor, N_coor, phonons_p->symplectic_transform_P ) ) info=1;

  /* conjugate symplectic transform P */
  if( MATRIX_ALLOCATE( N_coor, N_coor, phonons_p->conjugate_symplectic_transform_P ) ) info=1;


  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_phonons_free( phonons_p phonons_p ){

  /* dummies */
  int info=0;


  /* symplectic transform P */
  if( MATRIX_FREE( phonons_p->conjugate_symplectic_transform_P ) ) info=1;

  /* conjugate symplectic transform P */
  if( MATRIX_FREE( phonons_p->symplectic_transform_P ) ) info=1;

  /* symplectic transform R */
  if( MATRIX_FREE( phonons_p->conjugate_symplectic_transform_R ) ) info=1;

  /* conjugate symplectic transform R */
  if( MATRIX_FREE( phonons_p->symplectic_transform_R ) ) info=1;

  /* phononic_dos */
  if( RVECTOR_FREE( phonons_p->phononic_dos ) ) info=1;

  /* ap */
  if( RVECTOR_FREE( phonons_p->ap ) ) info=1;

  /* ar */
  if( RVECTOR_FREE( phonons_p->ar ) ) info=1;

  /* omega */
  if( RVECTOR_FREE( phonons_p->omega ) ) info=1;

  /* canonical frequencies */
  if( RVECTOR_FREE( phonons_p->canonical_frequencies ) ) info=1;

  /* canonical_transform */
  if( MATRIX_FREE( phonons_p->canonical_transform ) ) info=1;

  /* original_hessian */
  if( MATRIX_FREE( phonons_p->original_hessian ) ) info=1;

  /* initial_condition_atoms */
  if( VECTOR_FREE( phonons_p->initial_condition_atoms ) ) info=1;


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_phonons_read( FILE* fp, phonons_p phonons_p ){

  /* dummies */
  int info=0;


  if( VECTOR_READ( fp, phonons_p->initial_condition_atoms ) ) info=1;

  if( MATRIX_READ( fp, phonons_p->original_hessian ) ) info=1;

  if( MATRIX_READ( fp, phonons_p->canonical_transform ) ) info=1;

  if( RVECTOR_READ( fp, phonons_p->canonical_frequencies ) ) info=1;

  if( RVECTOR_READ( fp, phonons_p->omega ) ) info=1;

  if( RVECTOR_READ( fp, phonons_p->ar ) ) info=1;

  if( RVECTOR_READ( fp, phonons_p->ap ) ) info=1;

  if( RVECTOR_READ( fp, phonons_p->phononic_dos ) ) info=1;

  if( MATRIX_READ( fp, phonons_p->symplectic_transform_R ) ) info=1;

  if( MATRIX_READ( fp, phonons_p->conjugate_symplectic_transform_R ) ) info=1;

  if( MATRIX_READ( fp, phonons_p->symplectic_transform_P ) ) info=1;

  if( MATRIX_READ( fp, phonons_p->conjugate_symplectic_transform_P ) ) info=1;


  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_phonons_print( FILE* fp, const phonons phonons ){

  /* dummies */
  int info=0;


  /* phonons */
  if( VECTOR_PRINT( fp, phonons.initial_condition_atoms ) ) info=1;

  if( MATRIX_PRINT( fp, phonons.original_hessian ) ) info=1;

  if( MATRIX_PRINT( fp, phonons.canonical_transform ) ) info=1;

  if( RVECTOR_PRINT( fp, phonons.canonical_frequencies ) ) info=1;

  if( RVECTOR_PRINT( fp, phonons.omega ) ) info=1;

  if( RVECTOR_PRINT( fp, phonons.ar ) ) info=1;

  if( RVECTOR_PRINT( fp, phonons.ap ) ) info=1;

  if( RVECTOR_PRINT( fp, phonons.phononic_dos ) ) info=1;

  if( MATRIX_PRINT( fp, phonons.symplectic_transform_R ) ) info=1;

  if( MATRIX_PRINT( fp, phonons.conjugate_symplectic_transform_R ) ) info=1;

  if( MATRIX_PRINT( fp, phonons.symplectic_transform_P ) ) info=1;

  if( MATRIX_PRINT( fp, phonons.conjugate_symplectic_transform_P ) ) info=1;


  fflush( fp );


  return info;

}

//------------------------------------------

/* PRINT VERBOSE */

int PolyCEID_phonons_verbose_print( FILE* fp, const phonons phonons ){

  /* dummies */
  int info=0;


  /* initial_condition_atoms  */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# initial_condition_atoms:\n" );

  if( VECTOR_PRINT_PLUS( fp, phonons.initial_condition_atoms ) ) info=1;

  /* original_hessian */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# original_hessian:\n" );

  if( MATRIX_PRINT_PLUS( fp, phonons.original_hessian ) ) info=1;

  /* canonical_transform */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# canonical_transform:\n" );

  if( MATRIX_PRINT_PLUS( fp, phonons.canonical_transform ) ) info=1;

  /* canonical_frequencies */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# canonical_frequencies:\n" );

  if( RVECTOR_PRINT_PLUS( fp, phonons.canonical_frequencies ) ) info=1;

  /* omega */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# omega:\n" );

  if( RVECTOR_PRINT_PLUS( fp, phonons.omega ) ) info=1;

  /* ar */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# ar:\n" );

  if( RVECTOR_PRINT_PLUS( fp, phonons.ar ) ) info=1;

  /* ap */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# ap:\n" );

  if( RVECTOR_PRINT_PLUS( fp, phonons.ap ) ) info=1;

  /* phononic_dos */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# phononic_dos:\n" );

  if( RVECTOR_PRINT_PLUS( fp, phonons.phononic_dos ) ) info=1;

  /* symplectic_transform_R */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# symplectic_transform_R:\n" );

  if( MATRIX_PRINT_PLUS( fp, phonons.symplectic_transform_R ) ) info=1;

  /* conjugate_symplectic_transform R */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# conjugate_symplectic_transform_R:\n" );

  if( MATRIX_PRINT_PLUS( fp, phonons.conjugate_symplectic_transform_R ) ) info=1;

  /* symplectic_transform P */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# symplectic_transform_P:\n" );

  if( MATRIX_PRINT_PLUS( fp, phonons.symplectic_transform_P ) ) info=1;

  /* conjugate_symplectic_transform P */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# conjugate_symplectic_transform_P:\n" );

  if( MATRIX_PRINT_PLUS( fp, phonons.conjugate_symplectic_transform_P ) ) info=1;


  fflush( fp );


  return info;

}

//------------------------------------------

/* COPY [ c1 = c2 ] */

int PolyCEID_phonons_copy( phonons_p phonons_p, const phonons phonons ){

  /* dummies */
  int info=0;


  if( VECTOR_COPY( phonons_p->initial_condition_atoms, phonons.initial_condition_atoms ) ) info=1;

  if( MATRIX_COPY( phonons_p->original_hessian, phonons.original_hessian ) ) info=1;

  if( MATRIX_COPY( phonons_p->canonical_transform, phonons.canonical_transform ) ) info=1;

  if( RVECTOR_COPY( phonons_p->canonical_frequencies, phonons.canonical_frequencies ) ) info=1;

  if( RVECTOR_COPY( phonons_p->omega, phonons.omega ) ) info=1;

  if( RVECTOR_COPY( phonons_p->ar, phonons.ar ) ) info=1;

  if( RVECTOR_COPY( phonons_p->ap, phonons.ap ) ) info=1;

  if( RVECTOR_COPY( phonons_p->phononic_dos, phonons.phononic_dos ) ) info=1;

  if( MATRIX_COPY( phonons_p->symplectic_transform_R, phonons.symplectic_transform_R ) ) info=1;

  if( MATRIX_COPY( phonons_p->conjugate_symplectic_transform_R, phonons.conjugate_symplectic_transform_R ) ) info=1;

  if( MATRIX_COPY( phonons_p->symplectic_transform_P, phonons.symplectic_transform_P ) ) info=1;

  if( MATRIX_COPY( phonons_p->conjugate_symplectic_transform_P, phonons.conjugate_symplectic_transform_P ) ) info=1;


  return info;

}

//------------------------------------------
