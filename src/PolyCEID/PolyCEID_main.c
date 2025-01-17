
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

#include "input/PolyCEID_input_parsing.h"
#include "input/PolyCEID_start_file.h"
#include "initial/PolyCEID_initialise_state.h"
#include "evolution/PolyCEID_integrator.h"
#include "output/PolyCEID_output.h" 
#include "initial/PolyCEID_rho_extra_indexing.h"
#include "initial/PolyCEID_CI_table.h"
#include "output/PolyCEID_output_general.h"
#include <time.h>


int main( int argc, char* argv[] ){

  constants          constants;
  state              state;
  config             config_def;
  config             config_tmp;
  /* dummies */
  FILE*              fp;
  time_t             time1, time2;
  struct tm*         timeinfo;
  int                days, hours, mins;
  double             time_length=0.0e0;
  double             secs;
  rvector            initial_positions_saved;
  unsigned short int flag_static=1; //WARNING: at the beginning a static calculation is assumed
  unsigned short int flag_restart=0;
  int                N_levels_many_saved=-1;
  int                counter=0;
  int                info=0;


  /* print header */
  PRINT_HEADER( stdout );

  /* start time */
  time1 = time( NULL );
  fprintf( stdout, "#   The run has started on %s\n", asctime( timeinfo = localtime( &time1 ) ) );

  /* line parsing */
  if( argc !=2 ){

    fprintf( stderr, "USAGE %s \"input file name\"\n", argv[ 0 ] );
    fflush( stderr );

    info=1;

  }

  if( !info ){

    while( !info && flag_static ){

      // set counter to zero
      state.step_counter = 0;
 
      
#ifdef __DEBUG__
      
      fprintf( stdout, "# static calculation: %d\n\n", counter );
      
      if( counter ){
	
	fprintf( stdout, "# initial_positions_saved:\n" );
	if( RVECTOR_PRINT_PLUS( stdout, initial_positions_saved ) ) info=1;
	fprintf( stdout, "\n" );
	
      }

      //      fprintf( stdout, "input file name: %s\n\n", argv[ 1 ] );
           
#endif /* __DEBUG__ */

      
      /* read the input file */
      if( !info ){

	fp = fopen( argv[ 1 ], "r" );

	if( INPUT_PARSING( fp, constants, counter, initial_positions_saved ) ) info=1; //WARNING: if conter>0 use initial_positions_saved

	fclose( fp );

      }


      /* allocate indexing */
      if( !info ){

	if( INDEXING_ALLOCATE( constants, 0 ) ) info=1; 
	
      }

      
      /* allocate rho_indexing_extra */
      if( !info ){

	if( RHO_EXTRA_INDEXING_ALLOCATE( constants ) ) info=1; 

      }

      
      /* check indexing */
#ifdef __DEBUG__
      
      if( !info ){

	if( CHECK_INDEXING() ) info=1;

      }

      if( !info ){
      
	if( CHECK_RHO_EXTRA_INDEXING() ) info=1;
      
      }

#endif /* __DEBUG__*/

      
      /* compute_index-aux_translate */
      if( !info ){

	if( COMPUTE_RHO_INDEX_AUX_TRANSLATE( constants ) ) info=1;
	
      }


      /* initialise state */
      if( !info ){
    
	if( INITIALISE_STATE( constants, state, config_tmp, config_def ) ) info=1;

      }

#ifdef __DEBUG__
      
      fprintf( stdout, "# N_levels_many_saved = %d\n", N_levels_many_saved );
      fprintf( stdout, "# N_levels_many       = %d\n", constants.N_levels_many );
      fprintf( stdout, "\n" );

      if( !info ){
      
        fprintf( stdout, "# constants [at the beginning]\n" );
       
      }
      else{

        fprintf( stdout, "# constants [at the beginning --- after an error occurred]\n" );
         
      }
     
      if( CONSTANTS_VERBOSE_PRINT( stdout, constants ) ) info=1;
       
      fprintf( stdout, "\n" );


      if( !info ){

        fprintf( stdout, "# config & state [at the beginning]\n" );
       
      }
      else{
       
        fprintf( stdout, "# config & state [at the beginning --- after an error occurred]\n" );
         
      }
     
     
      if( CONFIG_VERBOSE_PRINT( stdout, config_def ) ) info=1;
       
      if( STATE_VERBOSE_PRINT( stdout, state ) ) info=1;
       
      if( CONFIG_VERBOSE_PRINT( stdout, config_tmp ) ) info=1;
       
      fprintf( stdout, "\n" );
      fflush( stdout );
      
#endif /* __DEBUG__ */


      /* check_pruning */
      if( 0 == constants.flag_relaxation || N_levels_many_saved == constants.N_levels_many ){

        if( counter ){

          if( RVECTOR_FREE( initial_positions_saved ) ) info=1;

        }

        flag_static=0;

      }
      else{

        if( !counter ){

          if( RVECTOR_ALLOCATE( constants.N_atoms, initial_positions_saved ) ) info=1;

        }

        if( RVECTOR_COPY( initial_positions_saved, config_def.atoms.positions ) ) info=1;

        N_levels_many_saved = constants.N_levels_many;
	
        counter++;

      }

      if( !flag_static ){

        /* in the case, restart */
        if( constants.flag_restart ){

          time_length = constants.time_length;

          flag_restart = constants.flag_restart;

          if( READ_START_FILE( constants, state, config_def ) )  info=1;

          // state.step_counter = (int) ( config_def.time /constants.dt +EPS );
	   
          constants.time_length = time_length;

          constants.flag_restart = flag_restart;

#ifdef __DEBUG__
      
	  fprintf( stdout, "# ----> config after restart\n" );      
          if( CONFIG_VERBOSE_PRINT( stdout, config_def ) ) info=1;
          fflush( stdout );

#endif /* __DEBUG__ */

        }

	/* opening output files */
	if( !info ){

	  if( OUTPUT_OPENING( constants ) ) info=1;

	}


	/* integration of the CEID EOM */
	if( !info ){
      
	  if( GLOBAL_CEID_INTEGRATOR( constants, state, config_tmp, config_def ) ) info=1;

	}


	/* closing output files */
	if( !info ){
      
	  if( OUTPUT_CLOSING( constants ) ) info=1;
	  
	}

      } /* end conditional flag_static*/    


      /* deallocate state */
      if( !info ){

	if( STATE_FREE( state ) ) info=1;
	
      }


      /* deallocate config_tmp */
      if( !info ){

	if( CONFIG_FREE( config_tmp ) ) info=1;
	
      }


      /* deallocate config_def */
      if( !info ){

	if( CONFIG_FREE( config_def ) ) info=1;
	
      }


      /* deallocate rho_extra_indexing */
      if( !info ){

	if( RHO_EXTRA_INDEXING_FREE() ) info=1; 

      }

      
      /* deallocate indexing */
      if( !info ){
    
	if( INDEXING_FREE( 0 ) ) info=1; 

      }


      /* deallocate constants */
      if( !info ){
    
	if( CONSTANTS_FREE( constants ) ) info=1;

      }
      
    } /* end while*/

  } /* end info conditional */


  /* end time */
  time2 = time( NULL );
  fprintf( stdout, "#   The run has ended on %s\n", asctime( localtime( &time2 ) ) );

  //secs  = difftime( time2, time1 ) /CLOCKS_PER_SEC;
  secs  = difftime( time2, time1 );
  
  days  = ceil( secs ) /DAY2SEC;
  secs -= days *DAY2SEC;
  
  hours = ceil( secs ) /HOUR2SEC;
  secs -= hours *HOUR2SEC;
  
  mins  = ceil( secs ) /MIN2SEC;
  secs -= mins *MIN2SEC;
  
  fprintf( stdout, "#   The run has lasted for %d days, %d hours, %d mins, and %4.2f secs\n\n", days, hours, mins, secs ); 


  return info;

}

