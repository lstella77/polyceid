
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
#include "PolyCEID_input_parsing.h"


#define MAX_N_FIELDS   10000


/*********************
  FUNCTIONS & MACROS
*********************/


/* INPUT PARSING */

int input_parsing( FILE* fp, constants_p constants_p, int counter, rvector initial_positions_saved ){

  /* dummies */
  int      pp[ NP_CONSTANTS ];
  int      i;
  list     input_list;
  double   dummy_seed;
  int      info=0;


  // reading list 
  LIST_READ( fp, &input_list );

  /*
    S_ECHO( "list before assigning variable", input_list.title.name );
    LIST_PRINT( stdout, &input_list );
    fflush( stdout );
  */
  

  // setting default for "output_label"
  strcpy( constants_p->output_label, "noname" );

  // extracting "output_label"
  if( ASSING_STRING_VARIABLE( &input_list, "output_label", constants_p->output_label, constants_p->output_label ) ) info=1;


  // extracting "flag_restart"
  constants_p->flag_restart = ASSING_FLAG_VARIABLE( &input_list, "restart" );


  // extracting "flag_relaxation"
  constants_p->flag_relaxation = ASSING_FLAG_VARIABLE( &input_list, "relaxation" );


  // extracting "flag_no_spin"
  constants_p->flag_no_spin = ASSING_FLAG_VARIABLE( &input_list, "no_spin" );


  // extracting "flag_no_spin_flip"
  constants_p->flag_no_spin_flip = ASSING_FLAG_VARIABLE( &input_list, "no_spin_flip" );


  // extracting "flag_normal_mode_expansion"
  constants_p->flag_normal_mode_expansion = ASSING_FLAG_VARIABLE( &input_list, "normal_mode_expansion" );


  // extracting "flag_no_Ehrenfest_frame"
  constants_p->flag_no_Ehrenfest_frame = ASSING_FLAG_VARIABLE( &input_list, "no_Ehrenfest_frame" );


  // extracting "periodic_boundary_condition"
  constants_p->flag_periodic_boundary_condition = ASSING_FLAG_VARIABLE( &input_list, "periodic_boundary_condition" );


  // extracting "N_atoms"
  if( ASSING_INT_VARIABLE( &input_list, "N_atoms", 1, &constants_p->N_atoms, NULL ) ) info=1;


  // extracting "N_levels_single"
  if( ASSING_INT_VARIABLE( &input_list, "N_levels_single", 1, &constants_p->N_levels_single, NULL ) ) info=1;


  // extract "N_electrons"
  if( ASSING_INT_VARIABLE( &input_list, "N_electrons", 1, &constants_p->N_electrons, NULL ) ) info=1;


  // setting default for "N_electrons_CAS"
  constants_p->N_electrons_CAS = constants_p->N_electrons;

  // extracting "N_electrons_CAS"
  if( ASSING_INT_VARIABLE( &input_list, "N_electrons_CAS", 1, &constants_p->N_electrons_CAS, &constants_p->N_electrons_CAS ) ) info=1;


  // setting default for "spacial_dimension"
  constants_p->spacial_dimension = 1;

  // extracting "spacial_dimension"
  if( ASSING_INT_VARIABLE( &input_list, "spacial_dimension", 1, &constants_p->spacial_dimension, &constants_p->spacial_dimension ) ) info=1;


  // setting default for "cell_dim"
  if( RVECTOR_ALLOCATE( constants_p->spacial_dimension, constants_p->cell_dim ) ) info=1;

  // extracting "cell_dim"
  if( ASSING_DOUBLE_VARIABLE( &input_list, "cell_dim", constants_p->spacial_dimension, constants_p->cell_dim.rvector, constants_p->cell_dim.rvector ) ) info=1;


  // setting default "N_coor_red"
  constants_p->N_coor_red = ( constants_p->spacial_dimension ) *( constants_p->N_atoms );

  // extracting "N_coor_red"
  if( ASSING_INT_VARIABLE( &input_list, "N_coor_red", 1, &constants_p->N_coor_red, &constants_p->N_coor_red ) ) info=1;


  // setting default for  "N_chain"
  constants_p->N_chain = 0;

  // extracting "N_chain"
  if( ASSING_INT_VARIABLE( &input_list, "N_chain", 1, &constants_p->N_chain,  &constants_p->N_chain ) ) info=1;


  if( constants_p->N_chain ){

    // setting default for "thermostat_masses"
    if( RVECTOR_ALLOCATE( constants_p->N_chain, constants_p->thermostat_masses ) ) info=1;

    // extracting "thermostat_masses"
    if( ASSING_DOUBLE_VARIABLE( &input_list, "thermostat_masses",  constants_p->N_chain, constants_p->thermostat_masses.rvector, constants_p->thermostat_masses.rvector ) ) info=1;

    // mass rescaling
    // BUGFIX: this is not very transparent...
    for( i=0; i< constants_p->N_chain; i++ ){

      constants_p->thermostat_masses.rvector[ i ] *= AMU_TO_INTERNAL;

    }       


    // setting default for "N_chain_steps"
    constants_p->N_chain_steps = 1;

    // extracting "N_chain_steps"
    if( ASSING_INT_VARIABLE( &input_list, "N_chain_steps", 1, &constants_p->N_chain_steps,  &constants_p->N_chain_steps ) ) info=1;


    // setting default for "temperature"
    constants_p->temperature=ZERO;

    // extracting "temperature"
    if( ASSING_DOUBLE_VARIABLE( &input_list, "temperature", 1, &constants_p->temperature, &constants_p->temperature ) ) info=1;


  } /* end N_chain */


  // setting default for "relevant_modes"
  if( IVECTOR_ALLOCATE( constants_p->N_coor_red, constants_p->relevant_modes ) ) info=1;

  for( i=0; i<( constants_p->N_coor_red ); i++ ){

    constants_p->relevant_modes.ivector[ i ] = (constants_p->N_coor_red) -1 -i;

  }        

  // extracting "relevant_modes"
  if( ASSING_INT_VARIABLE( &input_list, "relevant_modes", constants_p->N_coor_red, constants_p->relevant_modes.ivector, constants_p->relevant_modes.ivector ) ) info=1;


  // setting default for "CEID_order"
  constants_p->CEID_order = 0;

  // extracting "CEID_order"
  if( ASSING_INT_VARIABLE( &input_list, "CEID_order", 1, &constants_p->CEID_order,  &constants_p->CEID_order ) ) info=1;


  // setting default for "dt"
  constants_p->dt = ZERO;

  // extracting "dt"
  if( ASSING_DOUBLE_VARIABLE( &input_list, "dt", 1, &constants_p->dt, &constants_p->dt ) ) info=1;

    
  // setting default for "time_length"
  constants_p->time_length = ZERO;

  // extracting "time length"
  if( ASSING_DOUBLE_VARIABLE( &input_list, "time_length", 1, &constants_p->time_length, &constants_p->time_length ) ) info=1;

    
  // setting default for "skip_write"
  constants_p->skip_write = 1;

  // extracting "skip_write"
  if( ASSING_INT_VARIABLE( &input_list, "skip_write", 1, &constants_p->skip_write,  &constants_p->skip_write ) ) info=1;


  // setting default for "skip_save"
  constants_p->skip_save = 1;

  // extracting "skip_save"
  if( ASSING_INT_VARIABLE( &input_list, "skip_save", 1, &constants_p->skip_save,  &constants_p->skip_save ) ) info=1;


  // setting default for "initial_ionic_state"
  constants_p->initial_ionic_state = 0;

  // extracting "initial_ionic_state"
  if( ASSING_INT_VARIABLE( &input_list, "initial_ionic_state", 1, &constants_p->initial_ionic_state, &constants_p->initial_ionic_state ) ) info=1;

  // initial condition atoms
  if( ASSIGN_INITIAL_CONDITION_ATOMS( *constants_p, &input_list ) ) info=1;

  // setting default for "cutoff_energy"
  constants_p->cutoff_energy = 10.0;
  // BUGFIX: this value is totally arbitrary...

  // extracting "cutoff_energy"
  if( ASSING_DOUBLE_VARIABLE( &input_list, "cutoff_energy", 1, &constants_p->cutoff_energy, &constants_p->cutoff_energy ) ) info=1;


  // setting default for "initial_condition_type"
  strcpy( constants_p->initial_condition_type, "pure" );

  // extracting "initial_condition_type"
  if( ASSING_STRING_VARIABLE( &input_list, "initial_condition_type", constants_p->initial_condition_type, constants_p->initial_condition_type ) ) info=1;


  // setting default for "initial_many_body_occup"
  if( IVECTOR_ALLOCATE( 2, constants_p->initial_many_body_occup ) ) info=1;

  // extracting "initial_many_body_occup"
  if( ASSING_INT_VARIABLE( &input_list, "initial_many_body_occup", 2, constants_p->initial_many_body_occup.ivector, constants_p->initial_many_body_occup.ivector ) ) info=1;


  // setting default for "excited_many_body_occup"
  if( IVECTOR_ALLOCATE( 2, constants_p->excited_many_body_occup ) ) info=1;

  // extracting "excited_many_body_occup"
  if( ASSING_INT_VARIABLE( &input_list, "excited_many_body_occup", 2, constants_p->excited_many_body_occup.ivector, constants_p->excited_many_body_occup.ivector ) ) info=1;


  // setting default for "seed"
  dummy_seed = 12345;
  // BUGFIX: OK, this is a silly default value...

  // extracting "seed"
  if( ASSING_DOUBLE_VARIABLE( &input_list, "seed", 1, &dummy_seed, &dummy_seed ) ) info=1;

  // casting
  constants_p->seed = (unsigned long int) dummy_seed;

  // parsing hamiltonian
  if( HAMILTONIAN_PARSING( &input_list, *constants_p ) )info=1;

  /*
    S_ECHO( "list after assigning variable", input_list.title.name );
    LIST_PRINT( stdout, &input_list );
    fflush( stdout );
  */

  // free the input list
  LIST_CUT( &input_list );


  /*****************************
   * REMAINING INITIALISATIONS *
   *****************************/

  /* flag_pruning default */
  constants_p->flag_pruning=0;


  /* N_coor */
  constants_p->N_coor = ( constants_p->N_atoms ) *( constants_p->spacial_dimension );


  /* check */
  if(  constants_p->N_electrons_CAS > constants_p->N_electrons || ( constants_p->N_electrons -constants_p->N_electrons_CAS )%2 ){

    fprintf( stderr, "ERROR: Invalid value of N_electron_CAS [%d]\n", constants_p->N_electrons_CAS );
    fflush( stderr );
    
    info=1;

  }


  /* N_levels_single_CAS
   * This is twice the number of levels occupied by the CAS electrons
   *
   */
  constants_p->N_levels_single_CAS = MIN( constants_p->N_levels_single,
                                            2 *( ( constants_p->N_electrons_CAS +SPIN_DEG -1 ) /SPIN_DEG ) );


  /* N_levels_many */
  if( !info ){

    if( !constants_p->flag_no_spin ){
    // all the spin states possible

      /* BUGFIX: open shell just for 1 electron */
      if( 1 != constants_p->N_electrons){
	    
        fprintf( stderr, "ERROR: CI_table_spin for more than 1 electron not implemented yet!\n" );
        fflush( stderr );
  
        info=1;

      }
      else{

        if( constants_p->flag_no_spin_flip ){

          constants_p->N_levels_many = constants_p->N_levels_single;

        }
        else{

          constants_p->N_levels_many = SPIN_DEG *constants_p->N_levels_single;

        }

      }
   
    }
    else{
    // just singlet states possible

      /* BUGFIX: closed shell only */
      if( 2 != SPIN_DEG ){
    
        fprintf( stderr, "ERROR: CI_table_no_spin with SPIN_DEG != 2 not implemented yet!\n" );
        fflush( stderr );
    
        info=1;
	
      }
  
      if( ( constants_p->N_electrons_CAS )%2 ){
    
        fprintf( stderr, "ERROR: only closed shells with flag_no_spin = %u \n", constants_p->flag_no_spin );
        fflush( stderr );
    
        info=1;
	
      }

      i = CHOOSE( constants_p->N_levels_single_CAS, ( constants_p->N_electrons_CAS )/2 );


      if( constants_p->flag_no_spin_flip ){

        constants_p->N_levels_many = ( i *(i+1) ) /2;

      }
      else{

        constants_p->N_levels_many = i *i;

      }

    }

  }


  /* constants allocation */
  if( !info ){

    pp[0] = constants_p->N_atoms;
    pp[1] = constants_p->N_levels_single;
    pp[2] = constants_p->N_levels_many;
    pp[3] = constants_p->N_chain;

    if( CONSTANTS_ALLOCATE( NP_CONSTANTS, pp, *constants_p, 0 ) ) info=1; //WARNING: initial allocation without many_body_hybridisation_table and symmetry_coefficients

  }


#ifdef __NO_CEID__

  if( constants_p->CEID_order > 0 ){


    fprintf( stderr, "ERROR: CEID_order must be equal to 0 if __NO_CEID__ is in use\n" );

    info=1;

  }

#endif /* __NO_CEID__ */


  /* compute_rho_dimension */
  if( COMPUTE_RHO_DIMENSIONS( *constants_p, constants_p->CEID_order, constants_p->N_coor_red ) ) info=1;

  /* initial_condition_atoms_parsing */
  if( INITIAL_CONDITION_ATOMS_PARSING( *constants_p, counter, initial_positions_saved ) ) info=1;

  /* initial_condition_electrons_parsing */
  if( INITIAL_CONDITION_ELECTRONS_PARSING( *constants_p ) ) info=1;


  return info;

}

//------------------------------------------

int hamiltonian_parsing( list_p list_p, constants_p constants_p ){

  /* dummies */
  entry_p  entry_found_p=NULL;
  entry_p  dummy_entry_p=NULL;
  int      N_atoms=0;
  int      N_nodes=0;
  int      N_links=0;
  int      dim=0;
  double** list_nodes;
  double** list_links;
  imatrix  link_table;
  int      i=0;
  int      j=0;
  int      count=0;
  int      in_node, out_node;
  int      info=0;

  
  // extracting "class"
  if( ASSING_STRING_VARIABLE( list_p, "class", constants_p->hamiltonian.class, NULL ) ) info=1;


  // setting default for "N_par_node"
  constants_p->hamiltonian.N_par_node = 0;

  // extracting "N_par_node"
  if( ASSING_INT_VARIABLE( list_p, "N_par_node", 1, &constants_p->hamiltonian.N_par_node, &constants_p->hamiltonian.N_par_node ) ) info=1;


  // setting default for "N_par_link"
  constants_p->hamiltonian.N_par_link = 0;

  // extracting "N_par_link"
  if( ASSING_INT_VARIABLE( list_p, "N_par_link", 1, &constants_p->hamiltonian.N_par_link, &constants_p->hamiltonian.N_par_link ) ) info=1;


  // setting default for "N_par_extra"
  constants_p->hamiltonian.N_par_extra = 0;

  // extracting "N_par_extra"
  if( ASSING_INT_VARIABLE( list_p, "N_par_extra", 1, &constants_p->hamiltonian.N_par_extra, &constants_p->hamiltonian.N_par_extra ) ) info=1;


  N_atoms = constants_p->N_atoms;
  
  // list_nodes allocation
  dim = N_atoms;

  list_nodes = (double**) calloc( (size_t) dim, sizeof(double*) );

  for( i=0; i<dim; i++ ){

    list_nodes[ i ] = (double*) calloc( (size_t) constants_p->hamiltonian.N_par_node, sizeof(double) );

  }

  // list_links allocation
  dim = N_atoms *(N_atoms -1)/2;

  list_links = (double**) calloc( (size_t) dim, sizeof(double*) );

  for( i=0; i<dim; i++ ){

    list_links[ i ] = (double*) calloc( (size_t) constants_p->hamiltonian.N_par_link, sizeof(double) );

  }

  // link_table allocation
  if( IMATRIX_ALLOCATE( N_atoms, N_atoms, link_table ) ) info=1;


  // parse nodes
  N_nodes=0; 
  while( 1 ){

    entry_found_p = LIST_SEARCH( list_p, "node" );      

    if( !entry_found_p ) break;


    if( N_nodes > (N_atoms -1) ){

      fprintf( PolyCEID_STDERR, "ERROR: the number of nodes is larger than N_atoms [%d]\n", N_atoms );

      fflush( PolyCEID_STDERR );

      info=1;

      break;

    }        

    if( !entry_found_p->sublist_p ){

      fprintf( PolyCEID_STDERR, "ERROR: no value to extract from variable %s\n", entry_found_p->field.name );

      fflush( PolyCEID_STDERR );

      info=1;

    }
    else{

      if( entry_found_p->sublist_p->N_entries != constants_p->hamiltonian.N_par_node ){

        fprintf( PolyCEID_STDERR, "ERROR: node %d contains an inconsistent number of parameters: ", N_nodes );
        fprintf( PolyCEID_STDERR, "[ %d != %d]\n", entry_found_p->sublist_p->N_entries, constants_p->hamiltonian.N_par_node );

        fflush( PolyCEID_STDERR );

        info=1;

        break;

      }
      else{

        dummy_entry_p = entry_found_p->sublist_p->first_entry_p;

        // setting parameters
        for( i=0; i<(constants_p->hamiltonian.N_par_node); i++ ){

          list_nodes[ N_nodes ][ i ] = atof( dummy_entry_p->field.name );
	  
	  dummy_entry_p = dummy_entry_p->next_entry_p;

        }

        // cut the entry
        ENTRY_CUT( entry_found_p );

      }

      N_nodes++;

    }

  }

  // consistency checks
  if( !strcmp( constants_p->hamiltonian.class, "SSH" ) ){

    
    if( N_nodes != N_atoms ){

      fprintf( PolyCEID_STDERR, "ERROR: using class \"%s\", but N_nodes != N_atoms ", constants_p->hamiltonian.class );
      fprintf( PolyCEID_STDERR, "[ %d != %d]\n", N_nodes, N_atoms );

      fflush( PolyCEID_STDERR );

      info=1;

    }
    

  }        

  // set N_nodes
  constants_p->hamiltonian.N_nodes = N_nodes;

  if( constants_p->hamiltonian.N_nodes > 0 ){

    // allocate par_node: [N_nodes,N_par_node]
    constants_p->hamiltonian.par_node = (double**) calloc( (size_t) N_nodes, sizeof(double*) );

    for( i=0; i<N_nodes; i++ ){

      constants_p->hamiltonian.par_node[ i ] = (double*) calloc( (size_t) constants_p->hamiltonian.N_par_node, sizeof(double) );

    }        

    // setting par_node
    for( i=0; i<N_nodes; i++ ){

      for( j=0; j<(constants_p->hamiltonian.N_par_node); j++ ){

        constants_p->hamiltonian.par_node[ i ][ j ] = list_nodes[ i ][ j ];

      }

    }        

  }

  // list_nodes
  dim = N_atoms;

  for( i=0; i<dim; i++ ){

    free( list_nodes[ i ] );

  }

  free( list_nodes );


  // parse links
  N_links=0; 
  while( 1 ){

    entry_found_p = LIST_SEARCH( list_p, "link" );      

    if( !entry_found_p ) break;


    if( N_links > (N_atoms *(N_atoms -1)/2 ) ){

      fprintf( PolyCEID_STDERR, "ERROR: the number of links is larger than  N_atoms*(N_atoms-1)/2 [%d]\n", N_atoms *(N_atoms -1)/2  );

      fflush( PolyCEID_STDERR );

      info=1;

      break;

    }        

    if( !entry_found_p->sublist_p ){

      fprintf( PolyCEID_STDERR, "ERROR: no value to extract from variable %s\n", entry_found_p->field.name );

      fflush( PolyCEID_STDERR );

      info=1;

    }
    else{

      if( entry_found_p->sublist_p->N_entries != constants_p->hamiltonian.N_par_link +2 ){

        fprintf( PolyCEID_STDERR, "ERROR: link %d contains an inconsistent number of parameters: ", N_links );
        fprintf( PolyCEID_STDERR, "[ %d != %d]\n", entry_found_p->sublist_p->N_entries -2, constants_p->hamiltonian.N_par_link );

        fflush( PolyCEID_STDERR );

        info=1;

        break;

      }
      else{

        dummy_entry_p = entry_found_p->sublist_p->first_entry_p;

	// setting link_table
	in_node = atoi( dummy_entry_p->field.name );

	dummy_entry_p = dummy_entry_p->next_entry_p;

	out_node = atoi( dummy_entry_p->field.name );

	dummy_entry_p = dummy_entry_p->next_entry_p;

	//BUGFIX: I have to put N_links+1, otherwise it can be zero
        link_table.imatrix[ in_node ][ out_node ] = link_table.imatrix[ out_node ][ in_node ] = N_links+1;

        // setting parameters
        for( i=0; i<(constants_p->hamiltonian.N_par_link); i++ ){

          list_links[ N_links ][ i ] = atof( dummy_entry_p->field.name );
	  
	  dummy_entry_p = dummy_entry_p->next_entry_p;

        }

        // cut the entry
        ENTRY_CUT( entry_found_p );

      }

      N_links++;

    }

  }        

  // set N_links
  constants_p->hamiltonian.N_links = N_links;

  if( constants_p->hamiltonian.N_links > 0){

    // allocate par_link: [N_links,N_par_link]
    constants_p->hamiltonian.par_link = (double**) calloc( (size_t) N_links, sizeof(double*) );

    for( i=0; i<N_links; i++ ){

      constants_p->hamiltonian.par_link[ i ] = (double*) calloc( (size_t) constants_p->hamiltonian.N_par_link, sizeof(double) );

    }        

    // setting par_link
    for( i=0; i<N_links; i++ ){

      for( j=0; j<(constants_p->hamiltonian.N_par_link); j++ ){

        constants_p->hamiltonian.par_link[ i ][ j ] = list_links[ i ][ j ];

      }

    }        

  }

  // list_links
  dim = N_atoms *(N_atoms -1)/2;

  for( i=0; i<dim; i++ ){

    free( list_links[ i ] );

  }

  free( list_links);

  if( constants_p->hamiltonian.N_links > 0 ){

    // neighbour_index allocate
    constants_p->hamiltonian.neighbour_index = (int*) calloc( (size_t) (N_nodes+1), sizeof(int) );
  
    // neighbour_list allocate
    constants_p->hamiltonian.neighbour_list  = (int*) calloc( (size_t) (2 *N_links),   sizeof(int) );

    // neighbour_link allocate
    constants_p->hamiltonian.neighbour_link  = (int*) calloc( (size_t) (2 *N_links),   sizeof(int) );

    // setting neighbours
    constants_p->hamiltonian.neighbour_index[ 0 ] = count = 0;

    for( i=0; i<(constants_p->hamiltonian.N_nodes); i++ ){

      for( j=0; j<(constants_p->hamiltonian.N_nodes); j++ ){

        if( link_table.imatrix[ i ][ j ] ){
	
          constants_p->hamiltonian.neighbour_list[ count ] = j; 
	      
	  //BUGFIX: I have to put a "-1", because "0" means no link. See othe note above
          constants_p->hamiltonian.neighbour_link[ count ] = link_table.imatrix[ i ][ j ] -1;

          count++;

        }

      }

      constants_p->hamiltonian.neighbour_index[ i+1 ] = count; 

    }        

    /*
      fprintf( stdout, "# link_table\n");
      if( IMATRIX_PRINT( stdout, link_table ) ) info=1;
      fprintf( stdout, "\n");
  
      fprintf( stdout, "# hamiltonian\n");
      if( HAMILTONIAN_VERBOSE_PRINT( stdout, constants_p->hamiltonian) ) info=1;
      fprintf( stdout, "\n");
    */
  
  }

  // link_table deallocation
  if( IMATRIX_FREE( link_table ) ) info=1;


  // par_extra allocate
  if( constants_p->hamiltonian.N_par_extra > 0 ){

    constants_p->hamiltonian.par_extra = (double*) calloc( (size_t) constants_p->hamiltonian.N_par_extra, sizeof(double) );

    // extracting "par_extra"
    if( ASSING_DOUBLE_VARIABLE( list_p, "extra", constants_p->hamiltonian.N_par_extra, constants_p->hamiltonian.par_extra, constants_p->hamiltonian.par_extra ) ) info=1;

  }        


  return info;

}

//------------------------------------------

int assign_initial_condition_atoms( constants_p constants_p, list_p list_p ){

  /* constants */
  int      N_atoms;
  int      sdim;
  /* dummies */
  unsigned short int flag_fixed=0;
  entry_p  entry_found_p;
  int      pp_atoms[ NP_ATOMS ]; 
  atoms_p  atoms_p;
  int      comp;
  int      N_atoms_check;
  double*  dummy_position;
  double*  dummy_momentum;
  double   mass;
  int      info=0;


  N_atoms                     =  constants_p->N_atoms;
  sdim                        =  constants_p->spacial_dimension;
  atoms_p                     = &constants_p->initial_atoms;

  // setting default for "initial_atoms"
  pp_atoms[0] = N_atoms;
  pp_atoms[1] = sdim;

  if( ATOMS_ALLOCATE( NP_ATOMS, pp_atoms, *atoms_p ) );

  // dummy allocations
  dummy_position= ( double* ) calloc( (size_t) sdim, sizeof( double ) );

  dummy_momentum= ( double* ) calloc( (size_t) sdim, sizeof( double ) );

  // read atoms
  N_atoms_check=0;
  while( 1 ){

    entry_found_p = LIST_SEARCH( list_p, "atom" );      

    if( !entry_found_p || N_atoms_check > N_atoms ) break;

    if( !entry_found_p->sublist_p ){

      fprintf( PolyCEID_STDERR, "ERROR: no value to extract from variable %s\n", entry_found_p->field.name );

      fflush( PolyCEID_STDERR );

      info=1;
      
      break;

    }

    // setting default for "name"
    strcpy( atoms_p->names[ N_atoms_check ], "noname" );

    // extracting "name"
    if( ASSING_STRING_VARIABLE( entry_found_p->sublist_p, "name", atoms_p->names[ N_atoms_check ], atoms_p->names[ N_atoms_check ] ) ) info=1;

    // extracting "mass"    
    if( ASSING_DOUBLE_VARIABLE( entry_found_p->sublist_p, "mass", 1, &mass, NULL ) ) info=1;

    for( comp=0; comp<sdim; comp++ ){

      atoms_p->masses.rvector[ sdim *N_atoms_check +comp ] = mass;

    }        

    // mass rescaling
    // BUGFIX: this is not very transparent...
    atoms_p->masses.rvector[ N_atoms_check  ] *= AMU_TO_INTERNAL;


    // extracting "flag_fixed"
    flag_fixed = ASSING_FLAG_VARIABLE( entry_found_p->sublist_p, "fixed" );

    if( !flag_fixed ){

      for( comp=0; comp<sdim; comp++ ){ 

        atoms_p->mask.ivector[ N_atoms_check *sdim +comp ] = 1;

      }  

    }        

    // extracting "position"    
    if( ASSING_DOUBLE_VARIABLE( entry_found_p->sublist_p, "position", sdim, dummy_position, NULL ) ) info=1;

    // extracting "momentum"    
    if( ASSING_DOUBLE_VARIABLE( entry_found_p->sublist_p, "momentum", sdim, dummy_momentum, NULL ) ) info=1;

    // copy
    for( comp=0; comp<sdim; comp++ ){
      
      atoms_p->positions.rvector[ N_atoms_check *sdim +comp ] =  dummy_position[ comp ];

      atoms_p->momenta.rvector[ N_atoms_check *sdim +comp ]   =  dummy_momentum[ comp ];

    }        

    ENTRY_CUT( entry_found_p );    

    N_atoms_check++;

  } /* end while */        
  
  // dummy deallocate
  free( dummy_position );
 
  free( dummy_momentum );
  
  // consistency check
  if( N_atoms != N_atoms_check ){

    fprintf( stderr, "ERROR: %d atoms found instead of %d\n", N_atoms_check, N_atoms );	
    fflush( stderr );

    info=1;

  }        
  

  return info;

}

//------------------------------------------

int initial_condition_atoms_parsing( constants_p constants_p, int counter, rvector initial_positions_saved ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  int      N_atoms;
  int      N_coor;
  int      sdim;
  /* dummies */
  atoms_p  atoms_p;
  int      i, comp;
  int      index;
  double   dummy;
  int      info=0;


  flag_normal_mode_expansion  =  constants_p->flag_normal_mode_expansion;
  N_coor                      =  constants_p->N_coor;
  N_atoms                     =  constants_p->N_atoms;
  sdim                        =  constants_p->spacial_dimension;
  atoms_p                     = &constants_p->initial_atoms;


  // In the case, enforcing PBC
  if( constants_p->flag_periodic_boundary_condition ){

    for( i=0; i<N_atoms; i++ ){

      for( comp=0; comp<sdim; comp++ ){
        
        index = comp +i *sdim;
        
	atoms_p->positions.rvector[ index ] = 
	PBC_DIFF( atoms_p->positions.rvector[ index ], constants_p->cell_dim.rvector[ comp ] );

      }

    }	

  }

  /* initial positions */
  if( counter ){

    RVECTOR_COPY( atoms_p->positions, initial_positions_saved );

  }

  /* initial momenta -- consistency check */
  for( i=0; i<N_atoms; i++ ){

    for( comp=0; comp<sdim; comp++ ){
        
      index = comp +i *sdim;
   
      if( !atoms_p->mask.ivector[ index ] && fabs( atoms_p->momenta.rvector[ index ] ) > EPS ){
        
        fprintf( stderr, "ERROR: coor %d of atom %d is fixed, while its monetum is non-zero [%12.5le]\n", comp, i, atoms_p->momenta.rvector[ index ] );
	fflush( stderr );

	info=1;

      }        

    }

  }  


  // BUGFIX: the following is highly error prone...

  /* masses */
  counter=0;
  atoms_p->mass_tot=0.0e0;
  for( i=0; i<N_coor; i++ ){

    if( !atoms_p->mask.ivector[ i ] ){

      atoms_p->masses.rvector[ i ] = 0.0e0;

    }
    else{

      atoms_p->mass_tot += atoms_p->masses.rvector[ i ];

      counter++;

    }        

  }  
  
  atoms_p->mass_ave = atoms_p->mass_tot / (double) counter;

  atoms_p->mass_tot /= (double) sdim;
      

  /* mass_ave */
  dummy=0.0e0;
  for( i=0; i<N_coor; i++ ){
    
    if( atoms_p->mask.ivector[ i ] ){

      dummy += ( atoms_p->mass_ave ) /( atoms_p->masses.rvector[ i ] );

    }  
	
  }

  atoms_p->mass_ave = (double) counter /dummy *( atoms_p->mass_ave );
      

  /* masses_aux */
  for( i=0; i<N_coor; i++ ){
	
    if( flag_normal_mode_expansion ){

      atoms_p->masses_aux.rvector[ i ] = atoms_p->mass_ave;
	
    }
    else{

      atoms_p->masses_aux.rvector[ i ] = atoms_p->masses.rvector[ i ];

    }

  }
      

  /* compute centre of mass */
  if( COMPUTE_CENTRE_OF_MASS( *constants_p, *atoms_p ) ) info=1;


#ifdef __DEBUG__

  fprintf( stdout, "# initial_atoms\n" );
  if( ATOMS_VERBOSE_PRINT( stdout, *atoms_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG__*/


  return info;

}

//------------------------------------------
//------------------------------------------

int initial_condition_electrons_parsing( constants_p constants_p ){

  /* constants */
  int      N_levels_single;
  int*     N_levels_many_p;
  /* dummies */
  int      i, j, k;
  int      info=0;


  N_levels_single               =   constants_p->N_levels_single;
  N_levels_many_p               = &(constants_p->N_levels_many);


  /* create CI_table */
  if( !info ){
    
    if( !constants_p->flag_no_spin ){
    // all the spin states are allowed
    
      /* tmp */
      if( 1 != constants_p->N_electrons){
	
	fprintf( stderr, "ERROR: CI_table_spin for more than 1 electron not implemented yet!\n" );
	fflush( stderr );
	  
	info=1;
	  
      }  
      else{
	
	/* create CI_table_spin */
	if( CI_TABLE_SPIN( *constants_p ) ) info=1;
	
      }  
      
    }
    else{
    // just singlet states are allowed
      
      /* tmp */
      if( 2 != SPIN_DEG ){

	fprintf( stderr, "ERROR: CI_table_no_spin with SPIN_DEG != 2 not implemented yet!\n" );
	fflush( stderr );
	
	info=1;

      }
      else{

	/* create CI_table_no_spin */
	if( CI_TABLE_NO_SPIN( *constants_p ) ) info=1;

      }   
      
    }

  }


  if( IMATRIX_ZERO( constants_p->single_level_occupation_reduced ) ) info=1;
  
  if( !info ){
    
    for( i=0; i<*N_levels_many_p; i++ ){
      
      for( j=0; j<N_levels_single; j++ ){
	
	for( k=0; k<SPIN_DEG; k++ ){

	  constants_p->single_level_occupation_reduced.imatrix[ i ][ j ] += constants_p->single_level_occupation.imatrix[ i ][ k +j *SPIN_DEG ];

	} /* end k loop */

      } /* end j loop */

    } /* end i loop */

  }


  /* allocate_many_body_hybridisation_table */
  if( !info && constants_p->N_levels_many > 1 ){

    if( ALLOCATE_SYMMETRY_MULTIPLICITY( *constants_p ) ) info=1;

    if( ALLOCATE_MANY_BODY_HYBRIDISATION_TABLE(*constants_p ) ) info=1;

  }

      
  if( CONSTRUCT_TRANSITION( *constants_p, constants_p->initial_many_body_occup, constants_p->initial_many_body_state ) ) info=1;

  
  if( CONSTRUCT_TRANSITION( *constants_p, constants_p->excited_many_body_occup, constants_p->excited_many_body_state ) ) info=1;


  return info;

}

//------------------------------------------
/* 
 * BUGFIX: this function assume the ground-state is a single Slater
 *         determinant, filled according to the aufbau principle, e.g.,
 *         |1,1,1,1,0,0,0,0>
 *         SPIN_DEG=1,2 assumed
 */
int construct_transition( const constants constants, ivector occup, rvector_p transition_p ){

  /* constants */
  int      N_levels_single;
  int      N_levels_many;
  int      N_electrons;
  imatrix  single_level_occupation;
  imatrix  symmetry_multiplicity;
  /* dummy */
  ivector  ivec1;
  ivector  ivec2;
  int      hole, electron;
  int      i, h, k;
  int      level_hole, level_electron;
  int      ph_sign=1;
  unsigned short int flag_first=0;
  double   norm;
  int      info=0;


  N_levels_single          = constants.N_levels_single;
  N_levels_many            = constants.N_levels_many;
  N_electrons              = constants.N_electrons;
  single_level_occupation  = constants.single_level_occupation;
  symmetry_multiplicity    = constants.symmetry_multiplicity;


  hole     = occup.ivector[ 0 ];
  electron = occup.ivector[ 1 ];

  if( hole < 0 || electron < 0 ){

    fprintf( stderr, "ERROR: hole and electron must non-negative\n" );
    fflush( stderr );

    info=1;

  }


  if( !info ){

    if( IVECTOR_ALLOCATE( SPIN_DEG *N_levels_single, ivec1 ) ) info=1;

    if( IVECTOR_ALLOCATE( SPIN_DEG *N_levels_single, ivec2 ) ) info=1;

  }


  /* set to zero */
  if( RVECTOR_ZERO( *transition_p ) ) info=1;


  /* initial state */
  if( IVECTOR_ZERO( ivec1 ) ) info=1;


  /********************
   * the ground-state 
   ********************/
  for( i=0; i<N_electrons; i++ ){

    ivec1.ivector[ i ] = 1;

  }


  /*
    fprintf( stdout, "hole     = %d\n", hole );
    fprintf( stdout, "electron = %d\n", electron );
  */


  if( hole != 0 && electron != 0 ){

    level_hole = SPIN_DEG *( (N_electrons+1)/SPIN_DEG - hole );

    if( 1 == ivec1.ivector[ level_hole ] ){
      
      ivec1.ivector[ level_hole ] = 0; // the hole

    }
    else{
      
      fprintf( stderr, "ERROR: invalid hole [%d]\n", hole );
      fflush( stderr );
      
      info=1;

    }

    
    level_electron = SPIN_DEG *( (N_electrons+1)/2 +electron -1);

    if( 0 ==  ivec1.ivector[ level_electron] ){

      ivec1.ivector[ level_electron] = 1; // the electron

    }
    else{

      fprintf( stderr, "ERROR: invalid electron [%d]\n", hole );
      fflush( stderr );

      info=1;

    }

  }
  else if( ( hole == 0 && electron != 0 ) || ( hole != 0 && electron == 0 ) ){

    fprintf( stderr, "ERROR: invalid transition [%d -> %d]\n", hole, electron );
    fflush( stderr );

    info=1;

  }

  /*
    fprintf( stdout, "ivec1:\n" );
    if( IVECTOR_PRINT_PLUS( stdout, ivec1 ) ) info=1;
    fprintf( stdout, "\n" );
  */
  

  /* the search */
  for( i=0; i<N_levels_many; i++ ){

    /* ivector copy */
    for( k=0; k<(SPIN_DEG *N_levels_single); k++ ){
      
      ivec2.ivector[ k ] = single_level_occupation.imatrix[ i ][ k ];
      
    }

    if( !IVECTOR_COMPARE_MUTE( ivec1, ivec2 ) ){

      transition_p->rvector[ i ] = 1.0e0;
     
      flag_first=1;
     
    }

  } /* end i loop */

  /* 
    fprintf( stdout, "transition:\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *transition_p ) ) info=1;
    fprintf( stdout, "\n" );
  
    fprintf( stdout, "ivec2:\n" );
    if( IVECTOR_PRINT_PLUS( stdout, ivec2 ) ) info=1;
    fprintf( stdout, "\n" );
  */


  /************* 
   * spin flip * 
   *************/
  if( IVECTOR_COPY( ivec2, ivec1 ) ) info=1;

  for( k=0; k<N_levels_single; k++ ){

    for( h=0; h<SPIN_DEG; h++ ){

      ivec2.ivector[ k *SPIN_DEG +h ] = ivec1.ivector[ k *SPIN_DEG +SPIN_DEG -1 -h ];

    } /* end h loop */

  } /* end k loop */

  if( IVECTOR_COPY( ivec1, ivec2 ) ) info=1;

  /*
    fprintf( stdout, "transition:\n" );
    if( RVECTOR_PRINT_PLUS( stdout, *transition_p ) ) info=1;
    fprintf( stdout, "\n" );
  
    fprintf( stdout, "ivec1 [flipped]:\n" );
    if( IVECTOR_PRINT_PLUS( stdout, ivec1 ) ) info=1;
    fprintf( stdout, "\n" );
  */

  /* the search */
  for( i=0; i<N_levels_many; i++ ){

    /* ivector copy */
    for( k=0; k<(SPIN_DEG *N_levels_single); k++ ){
      
      ivec2.ivector[ k ] = single_level_occupation.imatrix[ i ][ k ];
      
    }

    if( !IVECTOR_COMPARE_MUTE( ivec1, ivec2 ) ){

      if( !flag_first ){
     
        transition_p->rvector[ i ] = 1.0e0;
     
        ph_sign = symmetry_multiplicity.imatrix[ i ][ 2 *N_SYMMETRIES ]; 

        flag_first=1;
        
      }
      else{

        transition_p->rvector[ i ] = (double) symmetry_multiplicity.imatrix[ i ][ 1 ];

      }        
      
    }
    
  } /* end i loop */


  /***************** 
   * particle-hole *
   *****************/
  if( hole != electron ){

    // p-h and inversion symmetries are ralated
    //if( (abs(hole -electron))%2 ){

      // is odd: even photon transition
      //ph_sign = +ph_sign;

    //}
    //else{

      // is even: odd photon transition
      //ph_sign = -ph_sign;

    //}  


    // assuming singlets have always -1 particle-hole character
    ph_sign = -ph_sign;


    if( IVECTOR_COPY( ivec2, ivec1 ) ) info=1;

      for( k=0; k<(SPIN_DEG *N_levels_single); k++ ){

        ivec2.ivector[ k ] = !ivec1.ivector[ (SPIN_DEG *N_levels_single -1) -k ];

      } /* end k loop */

      if( IVECTOR_COPY( ivec1, ivec2 ) ) info=1;

    /*
      fprintf( stdout, "transition:\n" );
      if( RVECTOR_PRINT_PLUS( stdout, *transition_p ) ) info=1;
      fprintf( stdout, "\n" );
  
      fprintf( stdout, "ivec1 [flipped]:\n" );
      if( IVECTOR_PRINT_PLUS( stdout, ivec1 ) ) info=1;
      fprintf( stdout, "\n" );
    */

      /* the search */
      for( i=0; i<N_levels_many; i++ ){

      /* ivector copy */
      for( k=0; k<(SPIN_DEG *N_levels_single); k++ ){
      
        ivec2.ivector[ k ] = single_level_occupation.imatrix[ i ][ k ];
      
      }

      if( !IVECTOR_COMPARE_MUTE( ivec1, ivec2 ) ){

        transition_p->rvector[ i ] = (double) ph_sign;
     
      }
    
    } /* end i loop */


    /******************* 
     * spin flip again * 
     *******************/
    if( IVECTOR_COPY( ivec2, ivec1 ) ) info=1;

    for( k=0; k<N_levels_single; k++ ){

      for( h=0; h<SPIN_DEG; h++ ){

        ivec2.ivector[ k *SPIN_DEG +h ] = ivec1.ivector[ k *SPIN_DEG +SPIN_DEG -1 -h ];

      } /* end h loop */

    } /* end k loop */

    if( IVECTOR_COPY( ivec1, ivec2 ) ) info=1;

    /*
      fprintf( stdout, "transition:\n" );
      if( RVECTOR_PRINT_PLUS( stdout, *transition_p ) ) info=1;
      fprintf( stdout, "\n" );
  
      fprintf( stdout, "ivec1 [flipped]:\n" );
      if( IVECTOR_PRINT_PLUS( stdout, ivec1 ) ) info=1;
      fprintf( stdout, "\n" );
    */

    /* the search */
    for( i=0; i<N_levels_many; i++ ){

      /* ivector copy */
      
      for( k=0; k<(SPIN_DEG *N_levels_single); k++ ){
      
        ivec2.ivector[ k ] = single_level_occupation.imatrix[ i ][ k ];
      
      }

      if( !IVECTOR_COMPARE_MUTE( ivec1, ivec2 ) ){

        transition_p->rvector[ i ] = (double) ph_sign *symmetry_multiplicity.imatrix[ i ][ 1 ];

      }
  
    } /* end particle-hole conditional */

  } /* end i loop */


  /* the normalisation */
  norm = 0.0e0;

  for( i=0; i<N_levels_many; i++ ){

    norm += ( transition_p->rvector[ i ] )*( transition_p->rvector[ i ] ); 
    
  }
  
  norm = sqrt( norm );
	
  if( norm > EPS ){
    
    for( i=0; i<N_levels_many; i++ ){

      transition_p->rvector[ i ] /= norm;
	      
    }
    
  }
  else{

    fprintf( stderr, "ERROR: norm of excited_many_body_state too small [%le]\n", norm );
    fflush( stderr );
    
    info=1;

  }


  if( IVECTOR_FREE( ivec2 ) ) info=1;

  if( IVECTOR_FREE( ivec1 ) ) info=1;


  return info;

}


//------------------------------------------

/* utilites */

//------------------------------------------


//------------------------------------------

