
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

#ifndef __PolyCEID_STRUCTURES_CONSTANTS__
#define __PolyCEID_STRUCTURES_CONSTANTS__


/* -------------------------------- */
/* __DEBUG_PLUS__ implies __DEBUG__ */
#ifdef __DEBUG_PLUS__

#ifndef __DEBUG__
#define __DEBUG__
#endif /* __DEBUG__ */

#endif /* __DEBUG_PLUS__ */
/* -------------------------------- */


#include "../algebra/PolyCEID_constants.h"
#include "../structures/PolyCEID_structures_atoms.h"
#include "../structures/PolyCEID_structures_hamiltonian.h"
#include "../input/PolyCEID_parsing.h"
#include "../algebra/PolyCEID_linear_algebra.h"
#include <limits.h>
#include <float.h>
#include <string.h>


/* Global parameters */

#define BOLTZMANN_K        8.61734215e-5
// eV / K

#define HBAR               0.65821188926e0
#define IHBAR              CMPLX_INIT( 0.0e0, HBAR )
// eV *fs

#define HBARC              1.97326968e3
// eV *A

#define AMU_TO_INTERNAL    103.6426867e0
// eV *fs^2 /A^2
// i.e. 1 amu = AMU_TO_INTERNAL eV * fs^2 / A^2

#define E_MASS_OVER_AMU    0.54857991e-3
// dimesionless

#define HARTREE            27.2113845e0
// eV

#define BOHR               0.5291772108e0
// A

#define COULOMB_CONST      HARTREE *BOHR
// e^2/( 4 \pi \epsilon_0 )
// eV *A

#define FINE_STR_CONST     7.297352568e-3
// dimensionless

#define SPIN_DEG           2

#define N_SYMMETRIES       1

#define NP_CONSTANTS       4


/*********************
  CONSTANTS
*********************/

typedef struct{

  /* input constants */
  char                 output_label[MAX_STRING_LENGTH];
  unsigned short int   flag_restart;                //
  unsigned short int   flag_relaxation;             //
  unsigned short int   flag_no_spin;                //
  unsigned short int   flag_no_spin_flip;           //
  unsigned short int   flag_normal_mode_expansion;  //
  unsigned short int   flag_no_Ehrenfest_frame;     //
  unsigned short int   flag_periodic_boundary_condition; //
  unsigned short int   flag_Ehrenfest;              //
  rvector              cell_dim;                    //
  int                  N_atoms;                     // Number of atoms
  int                  N_levels_single;             // Number of 1-electronic levels
  int                  N_electrons;                 // filling
  int                  N_electrons_CAS;             // filling
  int                  spacial_dimension;           //
  int                  N_coor_red;                  //
  int                  N_chain;                     //
  rvector              thermostat_masses;           //
  double               temperature;                 //
  int                  N_chain_steps;               //
  ivector              relevant_modes;              // [N_coor_red] vector with indeces of the effective coordinates
  int                  CEID_order;                  //
  double               dt;                          // integrators time_step (input)
  double               time_length;                 // simulation time
  int                  skip_write;                  // steps skipped before recording on file
  int                  skip_save;                   // steps skipped before saving the reference position
  char                 hamiltonian_label[MAX_STRING_LENGTH];
  int                  initial_ionic_state;
  atoms                initial_atoms;
  double               cutoff_energy;
  char                 initial_condition_type[MAX_STRING_LENGTH];
  char                 initial_many_body_occup[MAX_STRING_LENGTH];
  char                 excited_many_body_occup[MAX_STRING_LENGTH];
  unsigned long  int   seed;                        // random number generator seed
  hamiltonian          hamiltonian;                 // the hamiltonian
  unsigned short int   flag_observable_all;
  unsigned short int   flag_observable_geometry;
  unsigned short int   flag_observable_positions;
  unsigned short int   flag_observable_momenta;
  unsigned short int   flag_observable_forces;
  unsigned short int   flag_observable_populations;
  unsigned short int   flag_observable_mus;
  unsigned short int   flag_observable_purity;
  unsigned short int   flag_observable_energies;
  unsigned short int   flag_observable_adiabatic_populations;
  unsigned short int   flag_observable_projections;
  unsigned short int   flag_observable_adiabatic_pes_many;
  unsigned short int   flag_observable_adiabatic_pes_single;
  unsigned short int   flag_observable_single_level_populations;
  unsigned short int   flag_observable_nonadiabatic_couplings;
  unsigned short int   flag_observable_density_matrix;
  unsigned short int   flag_observable_transition_matrices;
  unsigned short int   flag_observable_adiabatic_states;
  unsigned short int   flag_observable_electronic_density_states;
  unsigned short int   flag_observable_ionic_density_states;
  unsigned short int   flag_observable_dipoles_many;
  unsigned short int   flag_observable_dipoles_single;

  /* derived constants */
  unsigned short int   flag_pruning;                //
  int                  N_coor;                      // N_atoms *sdim or effective number of atomic degrees of freedom
  int                  N_levels_single_CAS;              // Number of 1-electronic levels
  int                  N_levels_many;                    // Number of many-electronic levels
  imatrix              single_level_occupation;          // [N_levels_many][N_levels_single*SPIN_DEG]
  imatrix              single_level_occupation_reduced;  // [N_levels_many][N_levels_single]
  imatrix              many_body_hybridisation_table;    // [N_levels_many*(N_levels_many-1)/2][6]
  imatrix              symmetry_multiplicity;            // [N_levels_many][2 *N_SYMMETRIES +1]
  vector               symmetry_coefficient;             // [N_levels_many*(N_levels_many-1)/2]
  int                  max_rho_index;                    //
  int                  sqrt_max_rho_index;               //
  int                  rho_index_border_length;
  int                  rho_index_next_to_border_length;
  rvector              initial_many_body_state;
  rvector              excited_many_body_state;

} PolyCEID_constants;

/*********************
  POINTERS
*********************/

typedef PolyCEID_constants constants;

typedef PolyCEID_constants* constants_p;

/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_constants_allocate( int, int*, constants_p, unsigned short int );

#define CONSTANTS_ALLOCATE( np, pp, c, flag )  FUNCTION_CHECK( PolyCEID_constants_allocate( (np), (pp), &(c), (flag) ), PolyCEID_constants_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_constants_free( constants_p );

#define CONSTANTS_FREE( c )  FUNCTION_CHECK( PolyCEID_constants_free( &(c) ), PolyCEID_constants_free )

//------------------------------------------

/* READ */

int PolyCEID_constants_read( FILE*, constants_p );

#define CONSTANTS_READ( fp, c )  FUNCTION_CHECK( PolyCEID_constants_read( fp, &(c) ), PolyCEID_constants_read )

//------------------------------------------

/* PRINT */

int PolyCEID_constants_print( FILE*, const constants );

#define CONSTANTS_PRINT( fp, c )  FUNCTION_CHECK( PolyCEID_constants_print( (fp), (c) ), PolyCEID_constants_print )

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_constants_verbose_print( FILE*, const constants );

#define CONSTANTS_VERBOSE_PRINT( fp, c )  FUNCTION_CHECK( PolyCEID_constants_verbose_print( (fp), (c) ), PolyCEID_constants_verbose_print )

//------------------------------------------


/* COPY [ c1 = c2 ] */

int PolyCEID_constants_copy( constants_p, const constants );

#define CONSTANTS_COPY( c1, c2 )  FUNCTION_CHECK( PolyCEID_constants_copy( &(c1), (c2) ), PolyCEID_constants_copy )

//------------------------------------------


#endif /* __PolyCEID_STRUCTURES_CONSTANTS__ */

