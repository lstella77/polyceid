
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

#include "config.h"
#include "PolyCEID_structures_constants.h"


/* global */
int is_thermostat_constants=0;


/* ALLOCATE */

int PolyCEID_constants_allocate( int np, int* pp, constants_p constants_p, unsigned short int flag ){

  /* constants */
  int N_levels_single;
  int N_levels_many;
  int N_chain;
  /* dummies */
  int info=0;


  if( np != NP_CONSTANTS ) info=1;

  //N_atoms             = pp[0];
  N_levels_single     = pp[1];
  N_levels_many       = pp[2];
  N_chain             = pp[3];


  // is_thermostat_constants
  if( N_chain ){ 

    is_thermostat_constants=1;

  }
  else{

    is_thermostat_constants=0;

  }

  
  if( !flag ){

    /* single_level_occupation */
    if( IMATRIX_ALLOCATE( N_levels_many, (N_levels_single *SPIN_DEG), constants_p->single_level_occupation ) ) info=1;
    
    /* single_level_occupation_reduced */
    if( IMATRIX_ALLOCATE( N_levels_many, N_levels_single, constants_p->single_level_occupation_reduced ) ) info=1;

  }

  if( flag && N_levels_many > 1 ){

    /* many_body_table_hybridisation_table */
    if( IMATRIX_ALLOCATE( (N_levels_many*(N_levels_many -1)/2), 6, constants_p->many_body_hybridisation_table ) ) info=1;

    /* symmetry_coefficient */
    if( VECTOR_ALLOCATE( (N_levels_many*(N_levels_many -1)/2), constants_p->symmetry_coefficient ) ) info=1;

  }

  if( !flag ){

    /* symmetry_multiplicity  */
    if( IMATRIX_ALLOCATE( N_levels_many, 2*N_SYMMETRIES +1, constants_p->symmetry_multiplicity ) ) info=1;
    
    /* initial_many_body_state */
    if( RVECTOR_ALLOCATE( N_levels_many, constants_p->initial_many_body_state ) ) info=1;
    
    /* excited_many_body_state */
    if( RVECTOR_ALLOCATE( N_levels_many, constants_p->excited_many_body_state ) ) info=1;

  }


  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_constants_free( constants_p constants_p ){

  /* dummies */
  int info=0;


  /* initial_atoms */
  if( ATOMS_FREE( constants_p->initial_atoms ) ) info=1;

  /* excited_many_body_state */
  if( RVECTOR_FREE( constants_p->excited_many_body_state ) ) info=1;

  /* initial_many_body_state */
  if( RVECTOR_FREE( constants_p->initial_many_body_state ) ) info=1;
  
  /* excited_many_body_occup */
  if( IVECTOR_FREE( constants_p->excited_many_body_occup ) ) info=1;

  /* initial_many_body_occup */
  if( IVECTOR_FREE( constants_p->initial_many_body_occup ) ) info=1;
  
  /* symmetry_multiplicity  */
  if( IMATRIX_FREE( constants_p->symmetry_multiplicity ) ) info=1;

  if( constants_p->N_levels_many > 1 ){

    /* symmetry_coefficient */
    if( VECTOR_FREE( constants_p->symmetry_coefficient ) ) info=1;

    /* many_body_table_hybridisation_table */
    if( IMATRIX_FREE( constants_p->many_body_hybridisation_table ) ) info=1;

  }

  /* single_level_occupation_reduced */
  if( IMATRIX_FREE( constants_p->single_level_occupation_reduced ) ) info=1;

  /* single_level_occupation */
  if( IMATRIX_FREE( constants_p->single_level_occupation ) ) info=1;

  /* thermostat masses */
  if( is_thermostat_constants ){
	      
    /* masses */
    if( RVECTOR_FREE( constants_p->thermostat_masses ) ) info=1;

  }
			      
  /* cell_dim */
  if( RVECTOR_FREE( constants_p->cell_dim ) ) info=1;
  
  /* hamiltonian */
  if( HAMILTONIAN_FREE( constants_p->hamiltonian ) ) info=1;

  /* relevant_modes */
  if( IVECTOR_FREE( constants_p->relevant_modes ) ) info=1;


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_constants_read( FILE* fp, constants_p constants_p ){

  /* dummies */
  int info=0;


  /* general */
  if( fscanf( fp, "%s", constants_p->output_label ) < 1 ) info=1;

  if( fscanf( fp, "%lu", &constants_p->seed ) < 1 ) info=1;

  if( fscanf( fp, "%hu", &constants_p->flag_restart ) < 1 ) info=1;

  if( fscanf( fp, "%hu", &constants_p->flag_relaxation ) < 1 ) info=1;

  if( fscanf( fp, "%hu", &constants_p->flag_no_spin ) < 1 ) info=1;

  if( fscanf( fp, "%hu", &constants_p->flag_no_spin_flip ) < 1 ) info=1;

  if( fscanf( fp, "%hu", &constants_p->flag_normal_mode_expansion ) < 1 ) info=1;

  if( fscanf( fp, "%hu", &constants_p->flag_no_Ehrenfest_frame ) < 1 ) info=1;

  if( fscanf( fp, "%hu", &constants_p->flag_periodic_boundary_condition ) < 1 ) info=1;

  if( fscanf( fp, "%hu", &constants_p->flag_pruning ) < 1 ) info=1;

  /* atoms */
  if( fscanf( fp, "%d", &constants_p->N_atoms ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->N_coor ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->N_coor_red ) < 1 ) info=1;

  if( IVECTOR_READ( fp, constants_p->relevant_modes ) ) info=1;

  if( fscanf( fp, "%d", &constants_p->CEID_order ) ) info=1;

  if( HAMILTONIAN_READ( fp, constants_p->hamiltonian ) ) info=1;

  if( RVECTOR_READ( fp, constants_p->cell_dim ) ) info=1;

  /* electrons */
  if( fscanf( fp, "%d", &constants_p->N_levels_single ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->N_levels_single_CAS ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->N_levels_many ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->N_electrons ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->N_electrons_CAS ) < 1 ) info=1;

  if( IMATRIX_READ( fp, constants_p->single_level_occupation ) ) info=1;

  if( IMATRIX_READ( fp, constants_p->single_level_occupation_reduced ) ) info=1;

  if( constants_p->N_levels_many > 1 ){

    if( IMATRIX_READ( fp, constants_p->many_body_hybridisation_table ) ) info=1;

    if( VECTOR_READ( fp, constants_p->symmetry_coefficient ) ) info=1;

  }

  if( IMATRIX_READ( fp, constants_p->symmetry_multiplicity ) ) info=1;

  if( fscanf( fp, "%d", &constants_p->max_rho_index ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->sqrt_max_rho_index ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->rho_index_border_length ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->rho_index_next_to_border_length ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->initial_ionic_state ) < 1 ) info=1;

  if( fscanf( fp, "%s", constants_p->initial_condition_type ) < 1 ) info=1;

  if( !strcmp( constants_p->initial_condition_type, "pure" ) ){

    if( IVECTOR_READ( fp, constants_p->initial_many_body_occup ) ) info=1;

    if( IVECTOR_READ( fp, constants_p->excited_many_body_occup ) ) info=1;

    if( RVECTOR_READ( fp, constants_p->initial_many_body_state ) ) info=1;

    if( RVECTOR_READ( fp, constants_p->excited_many_body_state ) ) info=1;

  }  

  if( fscanf( fp, "%le", &constants_p->cutoff_energy ) < 1 ) info=1;

  /* time evolution */
  if( fscanf( fp, "%le", &constants_p->dt ) < 1 ) info=1;

  if( fscanf( fp, "%le", &constants_p->time_length ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->skip_write ) < 1 ) info=1;

  if( fscanf( fp, "%d", &constants_p->skip_save ) < 1 ) info=1;

  /* thermostat  */
  if( fscanf( fp, "%d",  &constants_p->N_chain         ) < 1 ) info=1;

  if( is_thermostat_constants ){
	  
    if( RVECTOR_READ( fp, constants_p->thermostat_masses ) ) info=1;

    if( fscanf( fp, "%le", &constants_p->temperature     ) < 1 ) info=1;

    if( fscanf( fp, "%d",  &constants_p->N_chain_steps ) < 1 ) info=1;

  }

  
  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_constants_print( FILE* fp, const constants constants ){

  /* dummies */
  int info=0;


  /* general */
  if( fprintf( fp, "%s\n", constants.output_label ) < 1 ) info=1;

  if( fprintf( fp, "%lu\n", constants.seed ) < 1 ) info=1;

  if( fprintf( fp, "%hu\n", constants.flag_restart ) < 1 ) info=1;

  if( fprintf( fp, "%hu\n", constants.flag_relaxation ) < 1 ) info=1;

  if( fprintf( fp, "%hu\n", constants.flag_no_spin ) < 1 ) info=1;

  if( fprintf( fp, "%hu\n", constants.flag_no_spin_flip ) < 1 ) info=1;

  if( fprintf( fp, "%hu\n", constants.flag_normal_mode_expansion ) < 1 ) info=1;

  if( fprintf( fp, "%hu\n", constants.flag_no_Ehrenfest_frame ) < 1 ) info=1;

  if( fprintf( fp, "%hu\n", constants.flag_periodic_boundary_condition ) < 1 ) info=1;

  if( fprintf( fp, "%hu\n", constants.flag_pruning ) < 1 ) info=1;

  /* atoms */
  if( fprintf( fp, "%d\n", constants.N_atoms ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.N_coor ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.N_coor_red ) < 1 ) info=1;

  if( IVECTOR_PRINT( fp, constants.relevant_modes ) ) info=1;

  if( fprintf( fp, "%d\n", constants.CEID_order ) < 1 ) info=1;
 
  if( HAMILTONIAN_PRINT( fp, constants.hamiltonian ) ) info=1;

  if( RVECTOR_PRINT( fp, constants.cell_dim ) ) info=1;

  /* electrons */
  if( fprintf( fp, "%d\n", constants.N_levels_single ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.N_levels_single_CAS ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.N_levels_many ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.N_electrons ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.N_electrons_CAS ) < 1 ) info=1;

  if( IMATRIX_PRINT( fp, constants.single_level_occupation ) < 1 ) info=1;

  if( IMATRIX_PRINT( fp, constants.single_level_occupation_reduced ) < 1 ) info=1;

  if( constants.N_levels_many > 1 ){

    if( IMATRIX_PRINT( fp, constants.many_body_hybridisation_table ) ) info=1;

    if( VECTOR_PRINT( fp, constants.symmetry_coefficient ) ) info=1;

  }

  if( IMATRIX_PRINT( fp, constants.symmetry_multiplicity ) ) info=1;

  if( fprintf( fp, "%d\n", constants.max_rho_index ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.sqrt_max_rho_index ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.rho_index_border_length ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.rho_index_next_to_border_length ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.initial_ionic_state ) ) info=1;

  if( fprintf( fp, "%s\n", constants.initial_condition_type ) < 1 ) info=1;

  if( !strcmp( constants.initial_condition_type, "pure" ) ){

    if( IVECTOR_PRINT( fp, constants.initial_many_body_occup ) ) info=1;

    if( IVECTOR_PRINT( fp, constants.excited_many_body_occup ) ) info=1;

    if( RVECTOR_PRINT( fp, constants.initial_many_body_state ) ) info=1;

    if( RVECTOR_PRINT( fp, constants.excited_many_body_state ) ) info=1;

  }  

  if( fprintf( fp, DOUBLE_FORMAT"\n", constants.cutoff_energy ) < 1 ) info=1;

  /* time evolution */
  if( fprintf( fp, DOUBLE_FORMAT"\n", constants.dt ) < 1 ) info=1;

  if( fprintf( fp, DOUBLE_FORMAT"\n", constants.time_length ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.skip_write ) < 1 ) info=1;

  if( fprintf( fp, "%d\n", constants.skip_save ) < 1 ) info=1;

  /* thermostat  */
  if( fprintf( fp, "%d\n", constants.N_chain          ) < 1 ) info=1;

  if( is_thermostat_constants ){
	  
    if( RVECTOR_PRINT( fp, constants.thermostat_masses ) ) info=1;

    if( fprintf( fp, DOUBLE_FORMAT"\n", constants.temperature     ) < 1 ) info=1;

    if( fprintf( fp, "%d\n",  constants.N_chain_steps ) < 1 ) info=1;

  }


  fflush( fp );


  return info;

}

//------------------------------------------

/* PRINT VERBOSE */

int PolyCEID_constants_verbose_print( FILE* fp, const constants constants ){

  /* dummies */
  int info=0;



#ifdef __DEBUG_PLUS__

  fprintf( fp, "# Preprocessor: using __DEBUG_PLUS__\n" );

#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG__

  fprintf( fp, "# Preprocessor: using __DEBUG__\n" );

#endif /* __DEBUG__ */


#ifdef __RK4__

  fprintf( fp, "# Preprocessor: using __RK4__\n" );

#endif /* __RK4__ */


#ifdef __NO_VERLET__

  fprintf( fp, "# Preprocessor: using __NO_VERLET__\n" );

#endif /* __NO_VERLET__ */


#ifdef __ENERGY_MOD__

  fprintf( fp, "# Preprocessor: using __ENERGY_MOD__\n" );

#endif /* __ENERGY_MOD__ */


#ifdef __NO_HESSIAN_CORRECTIONS__

  fprintf( fp, "# Preprocessor: using __NO_HESSIAN_CORRECTIONS__\n" );

#endif /* __NO_HESSIAN_CORRECTIONS__ */


#ifdef __NO_RELAXATION__

  fprintf( fp, "# Preprocessor: using __NO_RELAXATION__\n" );

#endif /* __NO_RELAXATION__ */


#ifdef __HAMILTONIAN_CONSISTENCY__

  fprintf( fp, "# Preprocessor: using __HAMILTONIAN_CONSISTENCY__\n" );

#endif /* __HAMILTONIAN_CONSISTENCY__ */


#ifdef __NO_ENERGY_CORRECTION__

  fprintf( fp, "# Preprocessor: using __NO_ENERGY_CORRECTION__\n" );

#endif /* __NO_ENERGY_CORRECTION__ */


#ifdef __CHECK_SYMMETRIES__

  fprintf( fp, "# Preprocessor: using __CHECK_SYMMETRIES__\n" );

#endif /* __CHECK_SYMMETRIES__ */


#ifdef __RHO_DOT_CHECK__

  fprintf( fp, "# Preprocessor: using __RHO_DOT_CHECK__\n" );

#endif /* __RHO_DOT_CHECK__ */


#ifdef __K_MATRIX_UPDATE__

  fprintf( fp, "# Preprocessor: using __K_MATRIX_UPDATE__\n" );

#endif /* __K_MATRIX_UPDATE__ */


#ifdef __NO_CEID__

  fprintf( fp, "# Preprocessor: using __NO_CEID__\n" );

#endif /* __NO_CEID__ */

#ifdef __ORTHONORMALISE__

  fprintf( fp, "# Preprocessor: using __ORTHONORMALISE__\n" );

#endif /* __ORTHONORMALISE__ */


#ifdef __CHECK_PARTICLE_HOLE__

  fprintf( fp, "# Preprocessor: using __CHECK_PARTICLE_HOLE__\n" );

#endif /*  __CHECK_PARTICLE_HOLE__ */


#ifdef __NO_DIPOLE__

  fprintf( fp, "# Preprocessor: using __NO_DIPOLE__\n" );

#endif /*  __NO_DIPOLE__ */


#ifdef __TIME_SCALE__

   fprintf( fp, "# Preprocessor: using __TIME_SCALE__\n" );

#endif /* __TIME_SCALE__ */

#ifdef __FORCES_NON_CONS__

   fprintf( fp, "# Preprocessor: using __FORCES_NON_CONS__\n" );

#endif /* __FORCES_NON_CONS__ */


  /* spin */
  fprintf( fp, "# SPIN_DEG = %d\n", SPIN_DEG );

  /* symmetries */
  fprintf( fp, "# N_SYMMETRIES = %d\n", N_SYMMETRIES );


  /* general */

  /* output_label */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# output_label:\n" );

  if( fprintf( fp, "# %s\n", constants.output_label ) < 1 ) info=1;

  /* seed */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# seed:\n" );

  if( fprintf( fp, "# %lu\n", constants.seed ) < 1 ) info=1;

  /* flag_restart */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# flag_restart:\n" );

  if( fprintf( fp, "# %hu\n", constants.flag_restart ) < 1 ) info=1;

  /* flag_relaxation */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# flag_relaxation:\n" );

  if( fprintf( fp, "# %hu\n", constants.flag_relaxation ) < 1 ) info=1;

  /* flag_no_spin */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# flag_no_spin:\n" );

  if( fprintf( fp, "# %hu\n", constants.flag_no_spin ) < 1 ) info=1;

  /* flag_no_spin_flip */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# flag_no_spin_flip:\n" );

  if( fprintf( fp, "# %hu\n", constants.flag_no_spin_flip ) < 1 ) info=1;

  /* flag_normal_mode_expansion */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# flag_normal_mode_expansion:\n" );

  if( fprintf( fp, "# %hu\n", constants.flag_normal_mode_expansion ) < 1 ) info=1;

  /* flag_no_Ehrenfest_frame */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# flag_no_Ehrenfest_frame:\n" );

  if( fprintf( fp, "# %hu\n", constants.flag_no_Ehrenfest_frame ) < 1 ) info=1;

  /* flag_periodic_boundary_condition */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# flag_periodic_boundary_condition:\n" );

  if( fprintf( fp, "# %hu\n", constants.flag_periodic_boundary_condition) < 1 ) info=1;

  /* flag_pruning */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# flag_pruning:\n" );

  if( fprintf( fp, "# %hu\n", constants.flag_pruning ) < 1 ) info=1;


  /* atoms */

  /* N_atoms */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_atoms:\n" );

  if( fprintf( fp, "# %d\n", constants.N_atoms ) < 1 ) info=1;

  /* N_coor */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_coor:\n" );

  if( fprintf( fp, "# %d\n", constants.N_coor ) < 1 ) info=1;

  /* N_coor_red */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_coor_red:\n" );

  if( fprintf( fp, "# %d\n", constants.N_coor_red ) < 1 ) info=1;

  /* relevant_modes */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# relevant_modes\n" );

  if( IVECTOR_PRINT_PLUS( fp, constants.relevant_modes ) ) info=1;

  /* CEID_order  */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# CEID_order:\n" );

  if( fprintf( fp, "# %d\n", constants.CEID_order ) < 1 ) info=1;

  /* hamiltonian */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# hamiltonian:\n" );

  if( HAMILTONIAN_VERBOSE_PRINT( fp, constants.hamiltonian ) ) info=1;

  /* cell_dim */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# cell_dim:\n" );

  if( RVECTOR_PRINT_PLUS( fp, constants.cell_dim ) ) info=1;

  /* electrons */

  /* N_levels_single */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_levels_single:\n" );

  if( fprintf( fp, "# %d\n", constants.N_levels_single ) < 1 ) info=1;

  /* N_levels_single_CAS */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_levels_single_CAS:\n" );

  if( fprintf( fp, "# %d\n", constants.N_levels_single_CAS ) < 1 ) info=1;

  /* N_levels_many */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_levels_many:\n" );

  if( fprintf( fp, "# %d\n", constants.N_levels_many ) < 1 ) info=1;

  /* N_electrons */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_electrons:\n" );

  if( fprintf( fp, "# %d\n", constants.N_electrons ) < 1 ) info=1;

  /* N_electrons_CAS */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_electrons_CAS:\n" );

  if( fprintf( fp, "# %d\n", constants.N_electrons_CAS ) < 1 ) info=1;

  /* single_level_occupation */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# single_level_occupation:\n" );

  if( IMATRIX_PRINT_PLUS( fp, constants.single_level_occupation ) ) info=1;

  /* single_level_occupation_reduced */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# single_level_occupation_reduced:\n" );

  if( IMATRIX_PRINT_PLUS( fp, constants.single_level_occupation_reduced ) ) info=1;

  if( constants.N_levels_many > 1 ){

    /* many_body_hybridisation_table */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# many_body_hybridisation_table:\n" );

    if( IMATRIX_PRINT_PLUS( fp, constants.many_body_hybridisation_table ) ) info=1;

    /* symmetry_coefficient */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# symmetry_coefficient:\n" );

    if( VECTOR_PRINT_PLUS( fp, constants.symmetry_coefficient ) ) info=1;

  }

  /* symmetry_multiplicity */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# symmetry_multiplicity:\n" );

  if( IMATRIX_PRINT_PLUS( fp, constants.symmetry_multiplicity ) ) info=1;

  /* max_rho_index */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# max_rho_index:\n" );

  if( fprintf( fp, "# %d\n", constants.max_rho_index ) < 1 ) info=1;

  /* sqrt_max_rho_index */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# sqrt_max_rho_index:\n" );

  if( fprintf( fp, "# %d\n", constants.sqrt_max_rho_index ) < 1 ) info=1;

  /* rho_index_border_length */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# rho_index_border_length:\n" );

  if( fprintf( fp, "# %d\n", constants.rho_index_border_length ) < 1 ) info=1;

  /* rho_index_next_to_border_length */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# rho_index_next_to_border_length:\n" );

  if( fprintf( fp, "# %d\n", constants.rho_index_next_to_border_length ) < 1 ) info=1;

  /* initial_ionic_state */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# initial_ionic_state:\n" );

  if( fprintf( fp, "# %d\n", constants.initial_ionic_state ) < 1 ) info=1;

  /* initial_condition_type */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# initial_condition_type:\n" );

  if( fprintf( fp, "# %s\n", constants.initial_condition_type ) < 1 ) info=1;

  if( !strcmp( constants.initial_condition_type, "pure" ) ){

    /* initial_many_body_occup */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# initial_many_body_occup:\n" );

    if( IVECTOR_PRINT_PLUS( fp, constants.initial_many_body_occup ) ) info=1;

    /* excited_many_body_occup */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# excited_many_body_occup:\n" );

    if( IVECTOR_PRINT_PLUS( fp, constants.excited_many_body_occup ) ) info=1;

    /* initial_many_body_state */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# initial_many_body_state:\n" );

    if( RVECTOR_PRINT_PLUS( fp, constants.initial_many_body_state ) ) info=1;

    /* excited_many_body_state */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# excited_many_body_state:\n" );

    if( RVECTOR_PRINT_PLUS( fp, constants.excited_many_body_state ) ) info=1;

  }

  /* cutoff_energy */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# cutoff_energy:\n" );

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", constants.cutoff_energy ) < 1 ) info=1;

  /* time evolution */

  /* dt */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# dt:\n" );

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", constants.dt ) < 1 ) info=1;

  /* time_length */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# time_length:\n" );

  if( fprintf( fp, "# "DOUBLE_FORMAT"\n", constants.time_length ) < 1 ) info=1;

  /* skip_write */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# skip write:\n" );

  if( fprintf( fp, "# %d\n", constants.skip_write ) < 1 ) info=1;

  /* skip_save */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# skip save:\n" );

  if( fprintf( fp, "# %d\n", constants.skip_save ) < 1 ) info=1;

  /* thermostat  */

  /* N_chain */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_chain:\n" );

  if( fprintf( fp, "# %d\n", constants.N_chain ) < 1 ) info=1;

  if( is_thermostat_constants ){

    /* thermostat masses */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# thermostat_masses:\n" );

    if( RVECTOR_PRINT_PLUS( fp, constants.thermostat_masses ) ) info=1;

    /* temperature */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# temperature:\n" );

    if( fprintf( fp, "# "DOUBLE_FORMAT"\n", constants.temperature ) < 1 ) info=1;

    /* N_chain_steps */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# N_chain_steps:\n" );

    if( fprintf( fp, "# %d\n", constants.N_chain_steps ) < 1 ) info=1;

  }  

  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "\n" );


  fflush( fp );


  return info;

}

//------------------------------------------

/* COPY [ c1 = c2 ] */

int PolyCEID_constants_copy( constants_p constants_p, const constants constants ){

  /* dummies */
  int info=0;


  /* general */
  if( !strcpy( constants_p->output_label, constants.output_label ) ) info=1;

  constants_p->seed                       = constants.seed;

  constants_p->flag_restart               = constants.flag_restart;

  constants_p->flag_relaxation            = constants.flag_relaxation;

  constants_p->flag_no_spin               = constants.flag_no_spin;

  constants_p->flag_no_spin_flip          = constants.flag_no_spin_flip;

  constants_p->flag_normal_mode_expansion = constants.flag_normal_mode_expansion;

  constants_p->flag_no_Ehrenfest_frame    = constants.flag_no_Ehrenfest_frame;

  constants_p->flag_periodic_boundary_condition = constants.flag_periodic_boundary_condition;

  constants_p->flag_pruning               = constants.flag_pruning;

  /* atoms */
  constants_p->N_atoms                    = constants.N_atoms;

  constants_p->N_coor                     = constants.N_coor;

  constants_p->N_coor_red                 = constants.N_coor_red;

  if( IVECTOR_COPY( constants_p->relevant_modes, constants.relevant_modes ) ) info=1;

  constants_p->CEID_order                 = constants.CEID_order;

  if( HAMILTONIAN_COPY( constants_p->hamiltonian, constants.hamiltonian ) ) info=1;

  if( RVECTOR_COPY( constants_p->cell_dim, constants.cell_dim ) ) info=1;

  /* electrons */
  constants_p->N_levels_single            = constants.N_levels_single;

  constants_p->N_levels_single_CAS        = constants.N_levels_single_CAS;

  constants_p->N_levels_many              = constants.N_levels_many;

  constants_p->N_electrons                = constants.N_electrons;

  constants_p->N_electrons_CAS            = constants.N_electrons_CAS;

  if( IMATRIX_COPY( constants_p->single_level_occupation, constants.single_level_occupation ) ) info=1;

  if( IMATRIX_COPY( constants_p->single_level_occupation_reduced, constants.single_level_occupation_reduced ) ) info=1;

  if( constants.N_levels_many > 1 ){

    if( IMATRIX_COPY( constants_p->many_body_hybridisation_table, constants.many_body_hybridisation_table ) ) info=1;

    if( VECTOR_COPY( constants_p->symmetry_coefficient, constants.symmetry_coefficient ) ) info=1;

  }

  if( IMATRIX_COPY( constants_p->symmetry_multiplicity, constants.symmetry_multiplicity ) ) info=1;

  constants_p->max_rho_index                   = constants.max_rho_index;

  constants_p->sqrt_max_rho_index              = constants.sqrt_max_rho_index;

  constants_p->rho_index_border_length         = constants.rho_index_border_length;

  constants_p->rho_index_next_to_border_length = constants.rho_index_next_to_border_length;

  constants_p->initial_ionic_state             = constants.initial_ionic_state;

  if( !strcpy( constants_p->initial_condition_type, constants.initial_condition_type ) ) info=1;

  if( !strcmp( constants.initial_condition_type, "pure" ) ){

    if( IVECTOR_COPY( constants_p->initial_many_body_occup, constants.initial_many_body_occup ) ) info=1;

    if( IVECTOR_COPY( constants_p->excited_many_body_occup, constants.excited_many_body_occup ) ) info=1;

    if( RVECTOR_COPY( constants_p->initial_many_body_state, constants.initial_many_body_state ) ) info=1;

    if( RVECTOR_COPY( constants_p->excited_many_body_state, constants.excited_many_body_state ) ) info=1;

  }  

  constants_p->cutoff_energy   = constants.cutoff_energy;

  /* time evolution */
  constants_p->dt              = constants.dt;

  constants_p->time_length     = constants.time_length;

  constants_p->skip_write      = constants.skip_write;

  constants_p->skip_save       = constants.skip_save;

  /* thermostat  */
  constants_p->N_chain         = constants.N_chain;

  if( is_thermostat_constants ){
	  
    if( RVECTOR_COPY( constants_p->thermostat_masses, constants.thermostat_masses ) ) info=1;

    constants_p->temperature     = constants.temperature;

    constants_p->N_chain_steps   = constants.N_chain_steps;

  }  


  return info;

}

//------------------------------------------
