
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

#ifndef __PolyCEID_STRUCTURES_STATE__
#define __PolyCEID_STRUCTURES_STATE__


#include <stdio.h>
#include "../structures/PolyCEID_structures_phonons.h"
#include "../structures/PolyCEID_structures_observables.h"


/* parameters */
#define NP_STATE  8


/*********************
  STATE
*********************/

typedef struct{

  long unsigned int     step_counter;                                     // integration step counter
  double                time;                                             // internal time
  //
  matrix                initial_rho_electron;                             //[N_levels_many, N_levels_many]
  matrix_array          rho_dot;                                          //[max_rho_index]
  matrix                Ehrenfest_frame_dot;                              //[N_levels_single,N_levels_single] 
  //
  double                rho_norm;
  double                rho_dot_norm;
  double                rho_projection;
  //
  rvector               adiabatic_populations;                            //[N_levels_many]
  matrix                adiabatic_projection;                             //[N_levels_many,N_levels_many]
  matrix                adiabatic_density_matrix;                         //[N_levels_many,N_levels_many]
  //
  rvector               single_level_populations_Ehrenfest;               //[N_levels_single]
  rvector               single_level_populations_adiabatic;               //[N_levels_single]
  //
  matrix_array          initial_state_recorded;                           //[max_rho_index]
  matrix_array          excited_state_recorded;                           //[max_rho_index]
  //
  matrix                one_body_electronic_density_matrix;               //[N_levels_single,N_levels_single]
  matrix                one_body_electronic_density_matrix_Ehrenfest;     //[N_levels_single,N_levels_single]
  //
  matrix                one_body_electronic_hole_matrix;                  //[N_levels_single,N_levels_single]
  matrix                one_body_electronic_hole_matrix_Ehrenfest;        //[N_levels_single,N_levels_single]
  matrix                one_body_electronic_particle_matrix;              //[N_levels_single,N_levels_single]
  matrix                one_body_electronic_particle_matrix_Ehrenfest;    //[N_levels_single,N_levels_single]
  //
  rvector               hole_energies;                                    //[N_level_single]
  rvector               particle_energies;                                //[N_level_single]
  //
  matrix                hole_orbitals;                                    //[N_level_single,N_levels_single]
  matrix                particle_orbitals;                                //[N_level_single,N_levels_single]
  //        
  rvector               hole_populations;                                 //[N_level_single]
  rvector               particle_populations;                             //[N_level_single]
  //
  rvector               adiabatic_PES_many;                               //[N_levels_many]
  matrix                adiabatic_states_many;                            //[N_levels_many,N_levels_many]
  //
  rvector               adiabatic_PES_single;                             //[N_levels_single]
  matrix                adiabatic_states_single;                          //[N_levels_single,N_levels_single]
  //
  rvector               electronic_density_eigenvalues;                   //[N_levels_many]
  matrix                electronic_density_eigenvectors;                  //[N_levels_many,N_levels_many]
  //
  matrix                ionic_density;                                    //[sqrt_max_rho_index,sqrt_max_rho_index]
  rvector               ionic_density_eigenvalues;                        //[sqrt_max_rho_index]
  matrix                ionic_density_eigenvectors;                       //[sqrt_max_rho_index,sqrt_max_rho_index]
  //
  rvector               nonadiabatic_coupling;                            //[N_coor]
  rvector               nonadiabatic_rate;                                //[N_coor]
  //
  phonons               phonons;
  //
  observables           observables;
  // Why here?
#ifdef __DEBUG__
  rvector               dummy_positions;      //[N_coor]
#endif /* __DEBUG__*/
  rvector               positions_saved;      //[N_coor]
  vector                delta_positions;      //[N_coor]
  rvector               momenta_saved;        //[N_coor]
  rvector               forces_saved;         //[N_coor]
  //
  matrix                H_matrix_saved;       //[N_levels_many,N_levels_many]
  matrix_array          F_matrix_saved;       //[N_coor,N_levels_many,N_levels_many]
  matrix_array          delta_F_matrix_saved; //[N_coor,N_levels_many,N_levels_many]
  matrix_array          K_matrix_saved;       //[N_coor*N_coor,N_levels_many,N_levels_many]
  //
  vector_array          F_matrix_saved_var;   //[N_levels_many*N_levels_many,N_coor]
  matrix_array          K_matrix_saved_var;   //[N_levels_many*N_levels_many,N_coor,N_coor]
  //
  rvector               dummy_rvector;        //[N_levels_many]
  vector                dummy_vector1;        //[N_levels_many]
  vector                dummy_vector2;        //[N_levels_many]
  vector                dummy_vector3;        //[N_levels_many]
  matrix                dummy_matrix1;        //[N_levels_many, N_levels_many]
  matrix                dummy_matrix2;        //[N_levels_many, N_levels_many]
  matrix                dummy_matrix3;        //[N_levels_many, N_levels_many]
  matrix                dummy_matrix4;        //[N_levels_many, N_levels_many]
  //
  vector                dummy_vector_single1; //[N_levels_single]
  vector                dummy_vector_single2; //[N_levels_single]
  vector                dummy_vector_single3; //[N_levels_single]
  rvector               dummy_rvector_single; //[N_levels_single]
  matrix                dummy_matrix_single1; //[N_levels_single, N_levels_single]
  matrix                dummy_matrix_single2; //[N_levels_single, N_levels_single]
  matrix                dummy_matrix_single3; //[N_levels_single, N_levels_single]
  matrix                dummy_matrix_single4; //[N_levels_single, N_levels_single]
  //
  matrix                dummy_Ehrenfest_frame;//[N_levels_single, N_levels_single]
  //
  vector                dummy_vector_coor1;   //[N_coor]
  vector                dummy_vector_coor2;   //[N_coor]
  matrix                dummy_matrix_coor1;   //[N_coor,N_coor]
  matrix                dummy_matrix_coor2;   //[N_coor,N_coor]

} PolyCEID_state;


/*********************
  POINTERS
*********************/

typedef PolyCEID_state  state;

typedef PolyCEID_state* state_p;


/*********************
  FUNCTIONS & MACROS
*********************/

/* ALLOCATE */

int PolyCEID_state_allocate( int, int*, state_p );

#define STATE_ALLOCATE( np, pp, s )  FUNCTION_CHECK( PolyCEID_state_allocate( (np), (pp), &(s) ), PolyCEID_state_allocate )

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_state_free( state_p );

#define STATE_FREE( s )  FUNCTION_CHECK( PolyCEID_state_free( &(s) ), PolyCEID_state_free )

//------------------------------------------

/* READ */

int PolyCEID_state_read( FILE*, state_p );

#define STATE_READ( fp, s )  FUNCTION_CHECK( PolyCEID_state_read( fp, &(s) ), PolyCEID_state_read )

//------------------------------------------

/* PRINT */

int PolyCEID_state_print( FILE*, const state );

#define STATE_PRINT( fp, s )  FUNCTION_CHECK( PolyCEID_state_print( (fp), (s) ), PolyCEID_state_print )

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_state_verbose_print( FILE*, const state );

#define STATE_VERBOSE_PRINT( fp, s )  FUNCTION_CHECK( PolyCEID_state_verbose_print( (fp), (s) ), PolyCEID_state_verbose_print )

//------------------------------------------


/* COPY [ s1 = s2 ] */

int PolyCEID_state_copy( state_p, const state );

#define STATE_COPY( s1, s2 )  FUNCTION_CHECK( PolyCEID_state_copy( &(s1), (s2) ), PolyCEID_state_copy )

//------------------------------------------


#endif /* __PolyCEID_STRUCTURES_STATE__ */

