
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
#include "PolyCEID_structures_state.h"

/* global */
int is_thermostat_state=0;


/* ALLOCATE */

int PolyCEID_state_allocate( int np, int* pp, state_p state_p ){

  /* state */
  int sdim;
  int N_coor;
  int N_levels_single;
  int N_levels_many;
  int max_rho_index;
  int sqrt_max_rho_index;
  int N_atoms;
  int N_chain;
  /* dummies */
  int pp_phonons[ NP_PHONONS ];
  int pp_observables[ NP_OBSERVABLES ];
  int info=0;


  if( np != NP_STATE ) info=1;

  /* read parameters */
  N_coor             = pp[0];
  N_levels_single    = pp[1];
  N_levels_many      = pp[2];
  max_rho_index      = pp[3];
  sqrt_max_rho_index = pp[4];
  N_atoms            = pp[5];
  N_chain            = pp[6];
  sdim               = pp[7];


  /* initial_rho_electrons */
  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->initial_rho_electron ) ) info=1;

  /* Ehrenfest_frame_dot */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->Ehrenfest_frame_dot ) ) info=1;

  /* rho_dot */
  if( MATRIX_ARRAY_ALLOCATE( max_rho_index, N_levels_many, N_levels_many, state_p->rho_dot ) ) info=1;

  /* adiabatic_populations */
  if( RVECTOR_ALLOCATE( N_levels_many, state_p->adiabatic_populations ) ) info=1;

  /* adiabatic_projection */
  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->adiabatic_projection ) ) info=1;

  /* adiabatic_density_matrix */
  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->adiabatic_density_matrix ) ) info=1;

  /* single_level_populations_Ehrenfest */
  if( RVECTOR_ALLOCATE( N_levels_single, state_p->single_level_populations_Ehrenfest ) ) info=1;

  /* single_level_populations_adiabatic */
  if( RVECTOR_ALLOCATE( N_levels_single, state_p->single_level_populations_adiabatic ) ) info=1;

  /* initial_state_recorded */
  if( MATRIX_ARRAY_ALLOCATE( max_rho_index +1, N_levels_many, N_levels_many, state_p->initial_state_recorded ) ) info=1;

  /* excited_state_recorded */
  if( MATRIX_ARRAY_ALLOCATE( max_rho_index +1, N_levels_many, N_levels_many, state_p->excited_state_recorded ) ) info=1;

  /* one_body_electronic_density_matrix */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->one_body_electronic_density_matrix ) ) info=1;

  /* one_body_electronic_density_matrix_Ehrenfest */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->one_body_electronic_density_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_hole_matrix */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->one_body_electronic_hole_matrix ) ) info=1;

  /* one_body_electronic_hole_matrix_Ehrenfest */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->one_body_electronic_hole_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_particle_matrix */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->one_body_electronic_particle_matrix ) ) info=1;

  /* one_body_electronic_particle_matrix_Ehrenfest */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->one_body_electronic_particle_matrix_Ehrenfest ) ) info=1;

  /* hole_energies */
  if( RVECTOR_ALLOCATE( N_levels_single, state_p->hole_energies ) )      info=1;

  /* particle_energies */
  if( RVECTOR_ALLOCATE( N_levels_single, state_p->particle_energies ) )  info=1;

  /* hole_orbitals */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->hole_orbitals ) )      info=1;

  /* particle_orbitals */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->particle_orbitals ) )  info=1;

  /* hole_populations */
  if( RVECTOR_ALLOCATE( N_levels_single, state_p->hole_populations ) )      info=1;

  /* particle_populations */
  if( RVECTOR_ALLOCATE( N_levels_single, state_p->particle_populations ) )  info=1;

  /* adiabatic_PES_many */
  if( RVECTOR_ALLOCATE( N_levels_many, state_p->adiabatic_PES_many ) )                            info=1;

  /* adiabatic_states_many */
  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->adiabatic_states_many ) )           info=1;

  /* adiabatic_PES_single */
  if( RVECTOR_ALLOCATE( N_levels_single, state_p->adiabatic_PES_single ) )                        info=1;

  /* adiabatic_states_single */
  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->adiabatic_states_single ) )     info=1;

  /* electronic_density_eigenvalues */
  if( RVECTOR_ALLOCATE( N_levels_many, state_p->electronic_density_eigenvalues ) )                info=1;

  /* electronic_density_eigenvectors */
  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->electronic_density_eigenvectors ) ) info=1;

  /* ionic_density */
  if( MATRIX_ALLOCATE( sqrt_max_rho_index, sqrt_max_rho_index, state_p->ionic_density ) )         info=1;

  /* ionic_density_eigenvalues */
  if( RVECTOR_ALLOCATE( sqrt_max_rho_index, state_p->ionic_density_eigenvalues ) )                info=1;

  /* ionic_density_eigenvectors */
  if( MATRIX_ALLOCATE( sqrt_max_rho_index, sqrt_max_rho_index, state_p->ionic_density_eigenvectors ) ) info=1;

  /* nonadiabaticity */
  if( RVECTOR_ALLOCATE( N_coor, state_p->nonadiabaticity ) ) info=1;

  /* nonadiabatic_coupling */
  if( RVECTOR_ALLOCATE( N_coor *N_levels_many *N_levels_many, state_p->nonadiabatic_coupling ) ) info=1;

  /* phonons */
  pp_phonons[0] = N_coor;
  pp_phonons[1] = sqrt_max_rho_index;

  if( PHONONS_ALLOCATE( NP_PHONONS, pp_phonons, state_p->phonons ) ) info=1;

  /* observables */
  pp_observables[0] = N_atoms;
  pp_observables[1] = N_levels_single;
  pp_observables[2] = N_levels_many;
  pp_observables[3] = N_chain;
  pp_observables[4] = sdim;;

  if( OBSERVABLES_ALLOCATE( NP_OBSERVABLES, pp_observables, state_p->observables ) ) info=1;

  /* dummies */
#ifdef __DEBUG__

  if( RVECTOR_ALLOCATE( N_coor, state_p->dummy_positions ) ) info=1;

#endif /* __DEBUG__ */

  if( RVECTOR_ALLOCATE( N_coor, state_p->positions_saved ) ) info=1;

  if( VECTOR_ALLOCATE( N_coor, state_p->delta_positions ) ) info=1;

  if( RVECTOR_ALLOCATE( N_coor, state_p->momenta_saved ) ) info=1;

  if( RVECTOR_ALLOCATE( N_coor, state_p->forces_saved ) ) info=1;

  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->H_matrix_saved ) )                info=1;

  if( MATRIX_ARRAY_ALLOCATE( N_coor, N_levels_many, N_levels_many, state_p->F_matrix_saved ) )  info=1;

  if( MATRIX_ARRAY_ALLOCATE( N_coor, N_levels_many, N_levels_many, state_p->delta_F_matrix_saved ) )  info=1;

  if( MATRIX_ARRAY_ALLOCATE( N_coor *N_coor, N_levels_many, N_levels_many, state_p->K_matrix_saved ) )  info=1;

  if( VECTOR_ARRAY_ALLOCATE( N_levels_many *N_levels_many, N_coor, state_p->F_matrix_saved_var ) )  info=1;

  if( MATRIX_ARRAY_ALLOCATE( N_levels_many *N_levels_many, N_coor, N_coor, state_p->K_matrix_saved_var ) )  info=1;

  if( RVECTOR_ALLOCATE( N_levels_many, state_p->dummy_rvector ) )                      info=1;

  if( VECTOR_ALLOCATE( N_levels_many, state_p->dummy_vector1 ) )                       info=1;

  if( VECTOR_ALLOCATE( N_levels_many, state_p->dummy_vector2 ) )                       info=1;

  if( VECTOR_ALLOCATE( N_levels_many, state_p->dummy_vector3 ) )                       info=1;

  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->dummy_matrix1 ) )        info=1;

  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->dummy_matrix2 ) )        info=1;

  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->dummy_matrix3 ) )        info=1;

  if( MATRIX_ALLOCATE( N_levels_many, N_levels_many, state_p->dummy_matrix4 ) )        info=1;

  if( VECTOR_ALLOCATE( N_levels_single, state_p->dummy_vector_single1 ) )              info=1;

  if( VECTOR_ALLOCATE( N_levels_single, state_p->dummy_vector_single2 ) )              info=1;

  if( VECTOR_ALLOCATE( N_levels_single, state_p->dummy_vector_single3 ) )              info=1;

  if( RVECTOR_ALLOCATE( N_levels_single, state_p->dummy_rvector_single ) )             info=1;

  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->dummy_matrix_single1 ) )  info=1;

  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->dummy_matrix_single2 ) )  info=1;

  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->dummy_matrix_single3 ) )  info=1;

  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->dummy_matrix_single4 ) )  info=1;

  if( MATRIX_ALLOCATE( N_levels_single, N_levels_single, state_p->dummy_Ehrenfest_frame ) ) info=1;


  if( VECTOR_ALLOCATE( N_coor, state_p->dummy_vector_coor1 ) )         info=1;

  if( VECTOR_ALLOCATE( N_coor, state_p->dummy_vector_coor2 ) )         info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, state_p->dummy_matrix_coor1 ) ) info=1;

  if( MATRIX_ALLOCATE( N_coor, N_coor, state_p->dummy_matrix_coor2 ) ) info=1;


  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_state_free( state_p state_p ){

  /* dummies */
  int info=0;


  /* dummies */

  if( MATRIX_FREE( state_p->dummy_matrix_coor2 ) ) info=1;

  if( MATRIX_FREE( state_p->dummy_matrix_coor1 ) ) info=1;

  if( VECTOR_FREE( state_p->dummy_vector_coor2 ) ) info=1;

  if( VECTOR_FREE( state_p->dummy_vector_coor1 ) ) info=1;

  if( MATRIX_FREE( state_p->dummy_Ehrenfest_frame ) )   info=1;

  if( MATRIX_FREE( state_p->dummy_matrix_single4 ) )    info=1;

  if( MATRIX_FREE( state_p->dummy_matrix_single3 ) )    info=1;

  if( MATRIX_FREE( state_p->dummy_matrix_single2 ) )    info=1;

  if( MATRIX_FREE( state_p->dummy_matrix_single1 ) )    info=1;

  if( RVECTOR_FREE( state_p->dummy_rvector_single ) )  info=1;

  if( VECTOR_FREE( state_p->dummy_vector_single3 ) )   info=1;

  if( VECTOR_FREE( state_p->dummy_vector_single2 ) )   info=1;

  if( VECTOR_FREE( state_p->dummy_vector_single1 ) )   info=1;

  if( MATRIX_FREE( state_p->dummy_matrix4 ) )          info=1;

  if( MATRIX_FREE( state_p->dummy_matrix3 ) )          info=1;

  if( MATRIX_FREE( state_p->dummy_matrix2 ) )          info=1;

  if( MATRIX_FREE( state_p->dummy_matrix1 ) )          info=1;

  if( VECTOR_FREE( state_p->dummy_vector3 ) )          info=1;

  if( VECTOR_FREE( state_p->dummy_vector2 ) )          info=1;

  if( VECTOR_FREE( state_p->dummy_vector1 ) )          info=1;

  if( RVECTOR_FREE( state_p->dummy_rvector ) )         info=1;

  if( MATRIX_ARRAY_FREE( state_p->K_matrix_saved_var ) )  info=1;

  if( VECTOR_ARRAY_FREE( state_p->F_matrix_saved_var ) )  info=1;

  if( MATRIX_ARRAY_FREE( state_p->K_matrix_saved ) )      info=1;

  if( MATRIX_ARRAY_FREE( state_p->delta_F_matrix_saved ) )  info=1;

  if( MATRIX_ARRAY_FREE( state_p->F_matrix_saved ) )      info=1;

  if( MATRIX_FREE( state_p->H_matrix_saved ) )            info=1;

  if( RVECTOR_FREE( state_p->forces_saved ) ) info=1;

  if( RVECTOR_FREE( state_p->momenta_saved ) ) info=1;

  if( VECTOR_FREE( state_p->delta_positions ) ) info=1;

  if( RVECTOR_FREE( state_p->positions_saved ) ) info=1;

#ifdef __DEBUG__

  if( RVECTOR_FREE( state_p->dummy_positions ) ) info=1;

#endif /* __DEBUG__ */

  /* observables */
  if( OBSERVABLES_FREE( state_p->observables ) ) info=1;

  /* phonons */
  if( PHONONS_FREE( state_p->phonons ) ) info=1;

  /* nonadiabatic_coupling */
  if( RVECTOR_FREE( state_p->nonadiabatic_coupling ) ) info=1;

  /* nonadiabaticity */
  if( RVECTOR_FREE( state_p->nonadiabaticity ) ) info=1;

  /* ionic_density_eigenvectors */
  if( MATRIX_FREE( state_p->ionic_density_eigenvectors ) )      info=1;

  /* ionic_density_eigenvalues */
  if( RVECTOR_FREE( state_p->ionic_density_eigenvalues ) )      info=1;

  /* ionic_density */
  if( MATRIX_FREE( state_p->ionic_density ) )                   info=1;

  /* electronic_density_eigenvectors */
  if( MATRIX_FREE( state_p->electronic_density_eigenvectors ) ) info=1;

  /* electronic_density_eigenvalues */
  if( RVECTOR_FREE( state_p->electronic_density_eigenvalues ) ) info=1;

  /* adiabatic_states_single */
  if( MATRIX_FREE( state_p->adiabatic_states_single ) )         info=1;

  /* adiabatic_PES_single */
  if( RVECTOR_FREE( state_p->adiabatic_PES_single ) )           info=1;

  /* adiabatic_states_many */
  if( MATRIX_FREE( state_p->adiabatic_states_many ) )           info=1;

  /* adiabatic_PES_many */
  if( RVECTOR_FREE( state_p->adiabatic_PES_many ) )             info=1;

  /* particle_populations */
  if( RVECTOR_FREE( state_p->particle_populations ) )              info=1;

  /* hole_populations */
  if( RVECTOR_FREE( state_p->hole_populations ) )                  info=1;

  /* particle_orbitals */
  if( MATRIX_FREE( state_p->particle_orbitals ) )              info=1;

  /* hole_orbitals */
  if( MATRIX_FREE( state_p->hole_orbitals ) )                  info=1;

  /* particle_energies */
  if( RVECTOR_FREE( state_p->particle_energies ) )              info=1;

  /* hole_energies */
  if( RVECTOR_FREE( state_p->hole_energies ) )                  info=1;

  /* one_body_electronic_particle_matrix_Ehrenfest */
  if( MATRIX_FREE( state_p->one_body_electronic_particle_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_particle_matrix */
  if( MATRIX_FREE( state_p->one_body_electronic_particle_matrix ) ) info=1;

  /* one_body_electronic_hole_matrix_Ehrenfest */
  if( MATRIX_FREE( state_p->one_body_electronic_hole_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_hole_matrix */
  if( MATRIX_FREE( state_p->one_body_electronic_hole_matrix ) ) info=1;

  /* one_body_electronic_density_matrix_Ehrenfest */
  if( MATRIX_FREE( state_p->one_body_electronic_density_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_density_matrix */
  if( MATRIX_FREE( state_p->one_body_electronic_density_matrix ) ) info=1;

  /* excited_state_recorded */
  if( MATRIX_ARRAY_FREE( state_p->excited_state_recorded ) ) info=1;

  /* initial_state_recorded */
  if( MATRIX_ARRAY_FREE( state_p->initial_state_recorded ) ) info=1;

  /* single_level_populations_adiabatic */
  if( RVECTOR_FREE( state_p->single_level_populations_adiabatic ) ) info=1;

  /* single_level_populations_Ehrenfest */
  if( RVECTOR_FREE( state_p->single_level_populations_Ehrenfest ) ) info=1;

  /* Ehrenfest_frame_dot */
  if( MATRIX_FREE( state_p->Ehrenfest_frame_dot ) ) info=1;

  /* adiabatic_density_matrix */
  if( MATRIX_FREE( state_p->adiabatic_density_matrix ) ) info=1;

  /* adiabatic_projection */
  if( MATRIX_FREE( state_p->adiabatic_projection ) ) info=1;

  /* adiabatic_populations */
  if( RVECTOR_FREE( state_p->adiabatic_populations ) ) info=1;

  /* rho_dot */
  if( MATRIX_ARRAY_FREE( state_p->rho_dot ) ) info=1;

  /* initial_rho_electron */
  if( MATRIX_FREE( state_p->initial_rho_electron ) ) info=1;


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_state_read( FILE* fp, state_p state_p ){

  /* dummies */
  int    info=0;


  /* step counter */
  if( fread( ( void* ) &state_p->step_counter, sizeof( unsigned long ), 1, fp ) < 1 ) info=1;

  /* initial_rho_electron */
  if( MATRIX_READ( fp, state_p->initial_rho_electron ) ) info=1;
  
  /* rho_dot */
  if( MATRIX_ARRAY_READ( fp, state_p->rho_dot ) ) info=1;

  /* Ehrenfest_frame_dot */
  if( MATRIX_READ( fp, state_p->Ehrenfest_frame_dot ) ) info=1;

  /* rho_norm */
  if( fread( ( void* ) &state_p->rho_norm, sizeof( double ), 1, fp ) < 1 ) info=1;

  /* rho_dot_norm */
  if( fread( ( void* ) &state_p->rho_dot_norm, sizeof( double ), 1, fp ) < 1 ) info=1;

  /* rho_projection */
  if( fread( ( void* ) &state_p->rho_projection, sizeof( double ), 1, fp ) < 1 ) info=1;

  /* adiabatic_populations */
  if( RVECTOR_READ( fp, state_p->adiabatic_populations ) ) info=1;

  /* adiabatic_projection */
  if( MATRIX_READ( fp, state_p->adiabatic_projection ) ) info=1;

  /* adiabatic_density_matrix */
  if( MATRIX_READ( fp, state_p->adiabatic_density_matrix ) ) info=1;

  /* single_level_populations_Ehrenfest */
  if( RVECTOR_READ( fp, state_p->single_level_populations_Ehrenfest ) ) info=1;

  /* single_level_populations_adiabatic */
  if( RVECTOR_READ( fp, state_p->single_level_populations_adiabatic ) ) info=1;
  
  /* initial_state_recorded */
  if( MATRIX_ARRAY_READ( fp, state_p->initial_state_recorded ) ) info=1;

  /* excited_state_recorded */
  if( MATRIX_ARRAY_READ( fp, state_p->excited_state_recorded ) ) info=1;

  /* one_body_electronic_density_matrix */
  if( MATRIX_READ( fp, state_p->one_body_electronic_density_matrix ) ) info=1;

  /* one_body_electronic_density_matrix_Ehrenfest */
  if( MATRIX_READ( fp, state_p->one_body_electronic_density_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_hole_matrix */
  if( MATRIX_READ( fp, state_p->one_body_electronic_hole_matrix ) ) info=1;

  /* one_body_electronic_hole_matrix_Ehrenfest */
  if( MATRIX_READ( fp, state_p->one_body_electronic_hole_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_particle_matrix */
  if( MATRIX_READ( fp, state_p->one_body_electronic_particle_matrix ) ) info=1;

  /* one_body_electronic_particle_matrix_Ehrenfest */
  if( MATRIX_READ( fp, state_p->one_body_electronic_particle_matrix_Ehrenfest ) ) info=1;

  /* hole_energies */
  if( RVECTOR_READ( fp, state_p->hole_energies ) )     info=1;

  /* particle_energies */
  if( RVECTOR_READ( fp, state_p->particle_energies ) ) info=1;

  /* hole_orbitals */
  if( MATRIX_READ( fp, state_p->hole_orbitals ) )     info=1;

  /* particle_orbitals */
  if( MATRIX_READ( fp, state_p->particle_orbitals ) ) info=1;

  /* hole_populations */
  if( RVECTOR_READ( fp, state_p->hole_populations ) )     info=1;

  /* particle_populations */
  if( RVECTOR_READ( fp, state_p->particle_populations ) ) info=1;

  /* adiabatic_PES_many */
  if( RVECTOR_READ( fp, state_p->adiabatic_PES_many ) )             info=1;

  /* adiabatic_states_many */
  if( MATRIX_READ( fp, state_p->adiabatic_states_many ) )           info=1;

  /* adiabatic_PES_single */
  if( RVECTOR_READ( fp, state_p->adiabatic_PES_single ) )           info=1;

  /* adiabatic_states_single */
  if( MATRIX_READ( fp, state_p->adiabatic_states_single ) )         info=1;

  /* electronic_density_eigenvalues */
  if( RVECTOR_READ( fp, state_p->electronic_density_eigenvalues ) ) info=1;

  /* electronic_density_eigenvectors */
  if( MATRIX_READ( fp, state_p->electronic_density_eigenvectors ) ) info=1;

  /* ionic_density */
  if( MATRIX_READ( fp, state_p->ionic_density ) )                   info=1;

  /* ionic_density_eigenvalues */
  if( RVECTOR_READ( fp, state_p->ionic_density_eigenvalues ) )      info=1;

  /* ionic_density_eigenvectors */
  if( MATRIX_READ( fp, state_p->ionic_density_eigenvectors ) )      info=1;

  /* nonadiabaticity */
  if( RVECTOR_READ( fp, state_p->nonadiabaticity ) )                info=1;

  /* nonadiabatic_coupling */
  if( RVECTOR_READ( fp, state_p->nonadiabatic_coupling ) )          info=1;

  /* phonons */
  if( PHONONS_READ( fp, state_p->phonons ) ) info=1;

  /* observables */
  if( OBSERVABLES_READ( fp, state_p->observables ) ) info=1;

  /* WARNING! dummy matrices are not recorded */


  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_state_print( FILE* fp, const state state ){

  /* dummies */
  int    info=0;


  /* step_counter */
  if( fwrite( ( const void* ) &state.step_counter, sizeof( unsigned long ), 1, fp ) < 1 ) info=1;	

  /* initial_rho_electron */
  if( MATRIX_PRINT( fp, state.initial_rho_electron ) ) info=1;

  /* rho_dot */
  if( MATRIX_ARRAY_PRINT( fp, state.rho_dot ) ) info=1;

  /* Ehrenfest_frame_dot */
  if( MATRIX_PRINT( fp, state.Ehrenfest_frame_dot ) ) info=1;

  /* rho_norm */
  if( fwrite( ( const void* ) &state.rho_norm, sizeof( double ), 1, fp ) < 1 ) info=1;

  /* rho_dot_norm */
  if( fwrite( ( const void* ) &state.rho_dot_norm, sizeof( double ), 1, fp ) < 1 ) info=1;

  /* rho_projection */
  if( fwrite( ( const void* ) &state.rho_projection, sizeof( double ), 1, fp ) < 1 ) info=1;

  /* adiabatic_populations */
  if( RVECTOR_PRINT( fp, state.adiabatic_populations ) ) info=1;

  /* adiabatic_projection */
  if( MATRIX_PRINT( fp, state.adiabatic_projection ) ) info=1;

  /* adiabatic_density_matrix */
  if( MATRIX_PRINT( fp, state.adiabatic_density_matrix ) ) info=1;

  /* single_level_populations_Ehrenfest */
  if( RVECTOR_PRINT( fp, state.single_level_populations_Ehrenfest ) ) info=1;

  /* single_level_populations_adiabatic */
  if( RVECTOR_PRINT( fp, state.single_level_populations_adiabatic ) ) info=1;

  /* initial_state_recorded */
  if( MATRIX_ARRAY_PRINT( fp, state.initial_state_recorded ) ) info=1;

  /* excited_state_recorded */
  if( MATRIX_ARRAY_PRINT( fp, state.excited_state_recorded ) ) info=1;

  /* one_body_electronic_density_matrix */
  if( MATRIX_PRINT( fp, state.one_body_electronic_density_matrix ) ) info=1;

  /* one_body_electronic_density_matrix_Ehrenfest */
  if( MATRIX_PRINT( fp, state.one_body_electronic_density_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_hole_matrix */
  if( MATRIX_PRINT( fp, state.one_body_electronic_hole_matrix ) ) info=1;

  /* one_body_electronic_hole_matrix_Ehrenfest */
  if( MATRIX_PRINT( fp, state.one_body_electronic_hole_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_particle_matrix */
  if( MATRIX_PRINT( fp, state.one_body_electronic_particle_matrix ) ) info=1;

  /* one_body_electronic_particle_matrix_Ehrenfest */
  if( MATRIX_PRINT( fp, state.one_body_electronic_particle_matrix_Ehrenfest ) ) info=1;

  /* hole_energies */
  if( RVECTOR_PRINT( fp, state.hole_energies ) )     info=1;

  /* particle_energies */
  if( RVECTOR_PRINT( fp, state.particle_energies ) ) info=1;

  /* hole_orbitals */
  if( MATRIX_PRINT( fp, state.hole_orbitals ) )     info=1;

  /* particle_orbitals */
  if( MATRIX_PRINT( fp, state.particle_orbitals ) ) info=1;

  /* hole_populations */
  if( RVECTOR_PRINT( fp, state.hole_populations ) )     info=1;

  /* particle_populations */
  if( RVECTOR_PRINT( fp, state.particle_populations ) ) info=1;

  /* adiabatic_PES_many */
  if( RVECTOR_PRINT( fp, state.adiabatic_PES_many ) )             info=1;

  /* adiabatic_states_many */
  if( MATRIX_PRINT( fp, state.adiabatic_states_many ) )           info=1;

  /* adiabatic_PES_single */
  if( RVECTOR_PRINT( fp, state.adiabatic_PES_single ) )           info=1;

  /* adiabatic_states_single */
  if( MATRIX_PRINT( fp, state.adiabatic_states_single ) )         info=1;

  /* electronic_density_eigenvalues */
  if( RVECTOR_PRINT( fp, state.electronic_density_eigenvalues ) ) info=1;

  /* electronic_density_eigenvectors */
  if( MATRIX_PRINT( fp, state.electronic_density_eigenvectors ) ) info=1;

  /* ionic_density */
  if( MATRIX_PRINT( fp, state.ionic_density ) )                   info=1;

  /* ionic_density_eigenvalues */
  if( RVECTOR_PRINT( fp, state.ionic_density_eigenvalues ) )      info=1;

  /* ionic_density_eigenvectors */
  if( MATRIX_PRINT( fp, state.ionic_density_eigenvectors ) )      info=1;

  /* nonadiabaticity */
  if( RVECTOR_PRINT( fp, state.nonadiabaticity ) )                info=1;

  /* nonadiabatic_coupling */
  if( RVECTOR_PRINT( fp, state.nonadiabatic_coupling ) )          info=1;

  /* phonons */
  if( PHONONS_PRINT( fp, state.phonons ) ) info=1;

  /* observables */
  if( OBSERVABLES_PRINT( fp, state.observables ) ) info=1;

  /* WARNING! dummy matrices are not recorded */

  fflush( fp );


  return info;

}

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_state_verbose_print( FILE* fp, const state state ){

  /* dummies */
  long unsigned int step_counter;
  int    info=0;


  step_counter = state.step_counter;

  /* step_counter */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# step_counter:\n");

  if( fprintf( fp, "# %lu\n", step_counter ) < 1 ) info=1;

  /* initial_rho_electron */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# initial_rho_electron:\n");

  if( MATRIX_PRINT_PLUS( fp, state.initial_rho_electron ) ) info=1;

  /* Ehrenfest_frame_dot */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# Ehrenfest_frame_dot:\n");

  if( MATRIX_PRINT_PLUS( fp, state.Ehrenfest_frame_dot ) ) info=1;

#ifdef __DEBUG_PLUS__

  /* rho_dot */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# rho_dot:\n");

  if( MATRIX_ARRAY_PRINT_PLUS( fp, state.rho_dot ) ) info=1;

#endif /* __DEBUG_PLUS__ */

  /* adiabatic_populations */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# adiabatic_populations:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.adiabatic_populations ) ) info=1;

  /* adiabatic_projection */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# adiabatic_projection:\n");

  if( MATRIX_PRINT_PLUS( fp, state.adiabatic_projection ) ) info=1;

  /* adiabatic_density_matrix */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# adiabatic_density_matrix:\n");

  if( MATRIX_PRINT_PLUS( fp, state.adiabatic_density_matrix ) ) info=1;

  /* single_level_populations_Ehrenfest */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# single_level_populations_Ehrenfest:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.single_level_populations_Ehrenfest ) ) info=1;

  /* single_level_populations_adiabatic */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# single_level_populations_adiabatic:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.single_level_populations_adiabatic ) ) info=1;

#ifdef __DEBUG_PLUS__

  /* initial_state_recorded */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# initial_state_recorded:\n");

  if( MATRIX_ARRAY_PRINT_PLUS( fp, state.initial_state_recorded ) ) info=1;

  /* excited_state_recorded */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# excited_state_recorded:\n");

  if( MATRIX_ARRAY_PRINT_PLUS( fp, state.excited_state_recorded ) ) info=1;

#endif /* __DEBUG_PLUS__ */

  /* one_body_electronic_density_matrix */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# one_body_electronic_density_matrix:\n");

  if( MATRIX_PRINT_PLUS( fp, state.one_body_electronic_density_matrix ) ) info=1;

  /* one_body_electronic_density_matrix_Ehrenfest */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# one_body_electronic_density_matrix_Ehrenfest:\n");

  if( MATRIX_PRINT_PLUS( fp, state.one_body_electronic_density_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_hole_matrix */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# one_body_electronic_hole_matrix:\n");

  if( MATRIX_PRINT_PLUS( fp, state.one_body_electronic_hole_matrix ) ) info=1;

  /* one_body_electronic_hole_matrix_Ehrenfest */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# one_body_electronic_hole_matrix_Ehrenfest:\n");

  if( MATRIX_PRINT_PLUS( fp, state.one_body_electronic_hole_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_particle_matrix */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# one_body_electronic_particle_matrix:\n");

  if( MATRIX_PRINT_PLUS( fp, state.one_body_electronic_particle_matrix ) ) info=1;

  /* one_body_electronic_particle_matrix_Ehrenfest */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# one_body_electronic_particle_matrix_Ehrenfest:\n");

  if( MATRIX_PRINT_PLUS( fp, state.one_body_electronic_particle_matrix_Ehrenfest ) ) info=1;

  /* hole_energies */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# hole_energies:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.hole_energies ) ) info=1;

  /* particle_energies */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# particle_energies:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.particle_energies ) ) info=1;

  /* hole_orbitals */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# hole_orbitals:\n");

  if( MATRIX_PRINT_PLUS( fp, state.hole_orbitals ) ) info=1;

  /* particle_orbitals */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# particle_orbitals:\n");

  if( MATRIX_PRINT_PLUS( fp, state.particle_orbitals ) ) info=1;

  /* hole_populations */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# hole_populations:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.hole_populations ) ) info=1;

  /* particle_populations */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# particle_populations:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.particle_populations ) ) info=1;

  /* adiabatic_PES_many */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# adiabatic_PES_many:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.adiabatic_PES_many ) ) info=1;

  /* adiabatic_states_many */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# adiabatic_states_many:\n");

  if( MATRIX_PRINT_PLUS( fp, state.adiabatic_states_many ) ) info=1;

  /* adiabatic_PES_single */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# adiabatic_PES_single:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.adiabatic_PES_single ) ) info=1;

  /* adiabatic_states_single */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# adiabatic_states_single:\n");

  if( MATRIX_PRINT_PLUS( fp, state.adiabatic_states_single ) ) info=1;

  /* electronic_density_eigenvalues */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# electronic_density_eigenvalues:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.electronic_density_eigenvalues ) ) info=1;

  /* electronic_density_eigenvectors */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# electronic_density_eigenvectors:\n");

  if( MATRIX_PRINT_PLUS( fp, state.electronic_density_eigenvectors ) ) info=1;

  /* ionic_density */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# ionic_density:\n");

  if( MATRIX_PRINT_PLUS( fp, state.ionic_density ) )                   info=1;

  /* ionic_density_eigenvalues */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# ionic_density_eigenvalues:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.ionic_density_eigenvalues ) )      info=1;

  /* ionic_density_eigenvectors */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# ionic_density_eigenvectors:\n");

  if( MATRIX_PRINT_PLUS( fp, state.ionic_density_eigenvectors ) )      info=1;

  /* nonadiabaticity */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# nonadiabaticity:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.nonadiabaticity ) )       info=1;

  /* nonadiabatic_coupling */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# nonadiabatic_coupling:\n");

  if( RVECTOR_PRINT_PLUS( fp, state.nonadiabatic_coupling ) ) info=1;

  /* rho_norm */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# rho_norm:\n");

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", state.rho_norm ) <1 ) info=1;

  /* rho_dot_norm */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# rho_dot_norm:\n");

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", state.rho_dot_norm ) <1 ) info=1;

  /* rho_projection */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# rho_projection:\n");

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", state.rho_projection ) <1 ) info=1;

  /* phonons */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# phonons:\n");

  if( PHONONS_VERBOSE_PRINT( fp, state.phonons ) ) info=1;

  /* observables */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# observables:\n");

  if( OBSERVABLES_VERBOSE_PRINT( fp, state.observables ) ) info=1;

  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "\n");


  /* WARNING! dummy matrices are not recorded */

  fflush( fp );


  return info;

}

//------------------------------------------

/* COPY [ s1 = s2 ] */

int PolyCEID_state_copy( state_p state_p, const state state ){

  /* dummies */
  int info=0;


  state_p->step_counter = state.step_counter;

  /* initial_rho_electron */
  if( MATRIX_COPY( state_p->initial_rho_electron, state.initial_rho_electron ) ) info=1;

  /* Ehrenfest_frame_dot */
  if( MATRIX_COPY( state_p->Ehrenfest_frame_dot, state.Ehrenfest_frame_dot ) ) info=1;

  /* rho_dot */
  if( MATRIX_ARRAY_COPY( state_p->rho_dot, state.rho_dot ) ) info=1;

  /* adiabatic_populations */
  if( RVECTOR_COPY( state_p->adiabatic_populations, state.adiabatic_populations ) ) info=1;

  /* adiabatic_projection */
  if( MATRIX_COPY( state_p->adiabatic_projection, state.adiabatic_projection ) ) info=1;

  /* adiabatic_density_matrix */
  if( MATRIX_COPY( state_p->adiabatic_density_matrix, state.adiabatic_density_matrix ) ) info=1;

  /* single_level_populations_Ehrenfest */
  if( RVECTOR_COPY( state_p->single_level_populations_Ehrenfest, state.single_level_populations_Ehrenfest ) ) info=1;

  /* single_level_populations_adiabatic */
  if( RVECTOR_COPY( state_p->single_level_populations_adiabatic, state.single_level_populations_adiabatic ) ) info=1;

  /* initial_state_recorded */
  if( MATRIX_ARRAY_COPY( state_p->initial_state_recorded, state.initial_state_recorded ) ) info=1;

  /* excited_state_recorded */
  if( MATRIX_ARRAY_COPY( state_p->excited_state_recorded, state.excited_state_recorded ) ) info=1;

  /* one_body_electronic_density_matrix */
  if( MATRIX_COPY( state_p->one_body_electronic_density_matrix, state.one_body_electronic_density_matrix ) ) info=1;

  /* one_body_electronic_density_matrix_Ehrenfest */
  if( MATRIX_COPY( state_p->one_body_electronic_density_matrix_Ehrenfest, state.one_body_electronic_density_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_hole_matrix */
  if( MATRIX_COPY( state_p->one_body_electronic_hole_matrix, state.one_body_electronic_hole_matrix ) ) info=1;

  /* one_body_electronic_hole_matrix_Ehrenfest */
  if( MATRIX_COPY( state_p->one_body_electronic_hole_matrix_Ehrenfest, state.one_body_electronic_hole_matrix_Ehrenfest ) ) info=1;

  /* one_body_electronic_particle_matrix */
  if( MATRIX_COPY( state_p->one_body_electronic_particle_matrix, state.one_body_electronic_particle_matrix ) ) info=1;

  /* one_body_electronic_particle_matrix_Ehrenfest */
  if( MATRIX_COPY( state_p->one_body_electronic_particle_matrix_Ehrenfest, state.one_body_electronic_particle_matrix_Ehrenfest ) ) info=1;

  /* hole_energies */
  if( RVECTOR_COPY( state_p->hole_energies, state.hole_energies ) ) info=1;

  /* particle_energies */
  if( RVECTOR_COPY( state_p->particle_energies, state.particle_energies ) ) info=1;

  /* hole_orbitals */
  if( MATRIX_COPY( state_p->hole_orbitals, state.hole_orbitals ) ) info=1;

  /* particle_orbitals */
  if( MATRIX_COPY( state_p->particle_orbitals, state.particle_orbitals ) ) info=1;

  /* hole_populations */
  if( RVECTOR_COPY( state_p->hole_populations, state.hole_populations ) ) info=1;

  /* particle_populations */
  if( RVECTOR_COPY( state_p->particle_populations, state.particle_populations ) ) info=1;

  /* adiabatic_PES_many */
  if( RVECTOR_COPY( state_p->adiabatic_PES_many, state.adiabatic_PES_many ) )                           info=1;

  /* adiabatic_states_many */
  if( MATRIX_COPY( state_p->adiabatic_states_many, state.adiabatic_states_many ) )                      info=1;

  /* adiabatic_PES_single */
  if( RVECTOR_COPY( state_p->adiabatic_PES_single, state.adiabatic_PES_single ) )                       info=1;

  /* adiabatic_states_single */
  if( MATRIX_COPY( state_p->adiabatic_states_single, state.adiabatic_states_single ) )                  info=1;

  /* electronic_density_eigenvalues */
  if( RVECTOR_COPY( state_p->electronic_density_eigenvalues, state.electronic_density_eigenvalues ) )   info=1;

  /* electronic_density_eigenvectors */
  if( MATRIX_COPY( state_p->electronic_density_eigenvectors, state.electronic_density_eigenvectors ) )  info=1;

  /* ionic_density  */
  if( MATRIX_COPY( state_p->ionic_density, state.ionic_density ) )                                      info=1;

  /* ionic_density_eigenvalues */
  if( RVECTOR_COPY( state_p->ionic_density_eigenvalues, state.ionic_density_eigenvalues ) )             info=1;

  /* ionic_density_eigenvectors */
  if( MATRIX_COPY( state_p->ionic_density_eigenvectors, state.ionic_density_eigenvectors ) )            info=1;

  /* nonadiabaticity */
  state_p->nonadiabaticity = state.nonadiabaticity;

  /* nonadiabatic_coupling */
  state_p->nonadiabatic_coupling = state.nonadiabatic_coupling;

  /* rho_norm */
  state_p->rho_norm = state.rho_norm;

  /* rho_dot_norm */
  state_p->rho_dot_norm = state.rho_dot_norm;

  /* rho_projection */
  state_p->rho_projection = state.rho_projection;

  /* phonons */
  if( PHONONS_COPY( state_p->phonons, state.phonons ) ) info=1;

  /* observables */
  if( OBSERVABLES_COPY( state_p->observables, state.observables ) ) info=1;

  /* WARNING! dummy matrices are not recorded */


  return info;

}

//------------------------------------------
