
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
#include "PolyCEID_structures_hamiltonian.h"


/* ALLOCATE */

int PolyCEID_hamiltonian_allocate( int np, int* pp, hamiltonian_p hamiltonian_p ){

  /* dummies */
  int info=0;


  return info;

}

//------------------------------------------

/* DEALLOCATE */

int PolyCEID_hamiltonian_free( hamiltonian_p hamiltonian_p ){

  /* dummies */
  int i=0;
  int info=0;


  /* par_extra */
  if( hamiltonian_p->N_par_extra > 0 ){

    free( hamiltonian_p->par_extra );

  }  
    
  /* par_link */
  if( hamiltonian_p->N_links > 0 ){

    for( i=0; i<hamiltonian_p->N_links; i++ ){

      free( hamiltonian_p->par_link[ i ] );
    
    }        

    free( hamiltonian_p->par_link );

  }  

  /* par_node */
  if( hamiltonian_p->N_nodes > 0 ){

    for( i=0; i<hamiltonian_p->N_nodes; i++ ){

      free( hamiltonian_p->par_node[ i ] );
    
    }        

    free( hamiltonian_p->par_node );

  }

  if( hamiltonian_p->N_links > 0 ){

    /* neighbour_link */
    free( hamiltonian_p->neighbour_link );
    
    /* neighbour_list */
    free( hamiltonian_p->neighbour_list );
    
    /* neighbour_index */
    free( hamiltonian_p->neighbour_index );

  }  


  return info;

}

//------------------------------------------

/* COPY [ h1 = h2 ] */

int PolyCEID_hamiltonian_copy( hamiltonian_p hamiltonian_p, const hamiltonian hamiltonian ){

  /* dummies */
  int info=0;


  return info;

}

//------------------------------------------

/* COMPARE [ h1 = h2 ] */

int PolyCEID_hamiltonian_compare( const hamiltonian hamiltonian1, const hamiltonian hamiltonian2 ){

  /* dummies */
  int info=0;


  return info;

}

//------------------------------------------

/* READ */

int PolyCEID_hamiltonian_read( FILE* fp, hamiltonian_p hamiltonian_p ){

  /* dummies */
  int    i,j;
  int    info=0;


  /* class */
  if( fread( ( void* ) hamiltonian_p->class, sizeof( char ), MAX_STRING_LENGTH, fp ) < 1 ) info=1;

  /* N_nodes */
  if( fread( ( void* ) &hamiltonian_p->N_nodes, sizeof( int ), 1, fp ) < 1 ) info=1;

  /* N_links */
  if( fread( ( void* ) &hamiltonian_p->N_links, sizeof( int ), 1, fp ) < 1 ) info=1;

  if( hamiltonian_p->N_links > 0 ){
  
    if( fread( ( void* ) &hamiltonian_p->neighbour_index, sizeof( int ), hamiltonian_p->N_nodes +1, fp ) < 1 ) info=1;

    if( fread( ( void* ) &hamiltonian_p->neighbour_list, sizeof( int ), 2 *hamiltonian_p->N_links, fp ) < 1 ) info=1;

    if( fread( ( void* ) &hamiltonian_p->neighbour_link, sizeof( int ), 2 *hamiltonian_p->N_links, fp ) < 1 ) info=1;

  }

  /* N_par_node */
  if( fread( ( void* ) &hamiltonian_p->N_par_node, sizeof( int ), 1, fp ) < 1 ) info=1;

  /* N_par_link */
  if( fread( ( void* ) &hamiltonian_p->N_par_link, sizeof( int ), 1, fp ) < 1 ) info=1;

  /* N_par_extra */
  if( fread( ( void* ) &hamiltonian_p->N_par_extra, sizeof( int ), 1, fp ) < 1 ) info=1;

  /* par_node */
  if( hamiltonian_p->N_nodes > 0 ){

    for( i=0; i<hamiltonian_p->N_nodes; i++ ){

      for( j=0; j<hamiltonian_p->N_par_node; j++ ){

        if( fread( ( void* ) &hamiltonian_p->par_node[i][j], sizeof( double ), 1, fp ) < 1 ) info=1;

      }        

    }

  }        

  /* par_link */
  if( hamiltonian_p->N_links > 0 ){

    for( i=0; i<hamiltonian_p->N_links; i++ ){

      for( j=0; j<hamiltonian_p->N_par_link; j++ ){

        if( fread( ( void* ) &hamiltonian_p->par_link[i][j], sizeof( double ), 1, fp ) < 1 ) info=1;

      }        

    }

  }

  /* par_extra */
  if( hamiltonian_p->N_par_extra > 0 ){

    if( fread( ( void* ) &hamiltonian_p->par_extra, sizeof( double ), hamiltonian_p->N_par_extra, fp ) < 1 ) info=1;
    
  }


  return info;

}

//------------------------------------------

/* PRINT */

int PolyCEID_hamiltonian_print( FILE* fp, const hamiltonian hamiltonian ){

  /* dummies */
  int    i,j;
  int    info=0;


  /* class */
  if( fwrite( ( const void* ) hamiltonian.class, sizeof( char ), MAX_STRING_LENGTH, fp ) < 1 ) info=1;

  /* N_nodes */
  if( fwrite( ( const void* ) &hamiltonian.N_nodes, sizeof( int ), 1, fp ) < 1 ) info=1;

  /* N_links */
  if( fwrite( ( const void* ) &hamiltonian.N_links, sizeof( int ), 1, fp ) < 1 ) info=1;

  if( hamiltonian.N_links > 0 ){
  
    if( fwrite( ( const void* ) &hamiltonian.neighbour_index, sizeof( int ), hamiltonian.N_nodes +1, fp ) < 1 ) info=1;

    if( fwrite( ( const void* ) &hamiltonian.neighbour_list, sizeof( int ), 2 *hamiltonian.N_links, fp ) < 1 ) info=1;

    if( fwrite( ( const void* ) &hamiltonian.neighbour_link, sizeof( int ), 2 *hamiltonian.N_links, fp ) < 1 ) info=1;

  }

  /* N_par_node */
  if( fwrite( ( const void* ) &hamiltonian.N_par_node, sizeof( int ), 1, fp ) < 1 ) info=1;

  /* N_par_link */
  if( fwrite( ( const void* ) &hamiltonian.N_par_link, sizeof( int ), 1, fp ) < 1 ) info=1;

  /* N_par_extra */
  if( fwrite( ( const void* ) &hamiltonian.N_par_extra, sizeof( int ), 1, fp ) < 1 ) info=1;

  /* par_node */
  if( hamiltonian.N_nodes > 0 ){

    for( i=0; i<hamiltonian.N_nodes; i++ ){

      for( j=0; j<hamiltonian.N_par_node; j++ ){

        if( fwrite( ( const void* ) &hamiltonian.par_node[i][j], sizeof( double ), 1, fp ) < 1 ) info=1;

      }        

    }

  }        

  /* par_link */
  if( hamiltonian.N_links > 0 ){

    for( i=0; i<hamiltonian.N_links; i++ ){

      for( j=0; j<hamiltonian.N_par_link; j++ ){

        if( fwrite( ( const void* ) &hamiltonian.par_link[i][j], sizeof( double ), 1, fp ) < 1 ) info=1;

      }        

    }

  }

  /* par_extra */
  if( hamiltonian.N_par_extra > 0 ){

    if( fwrite( ( const void* ) &hamiltonian.par_extra, sizeof( double ), hamiltonian.N_par_extra, fp ) < 1 ) info=1;
    
  }


  return info;

}

//------------------------------------------

/* VERBOSE PRINT */

int PolyCEID_hamiltonian_verbose_print( FILE* fp, const hamiltonian hamiltonian ){

  /* dummies */
  int i,j;
  int info=0;


  /* class */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# class:\n" );

  if( fprintf( fp, "# %s\n", hamiltonian.class ) < 1 ) info=1;

  /* N_nodes */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_nodes:\n" );

  if( fprintf( fp, "# %d\n", hamiltonian.N_nodes ) < 1 ) info=1;

  /* N_links */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_links:\n" );

  if( fprintf( fp, "# %d\n", hamiltonian.N_links ) < 1 ) info=1;

  if( hamiltonian.N_links > 0 ){
  
    /* neighbour_index */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# neighbour_index:\n" );

    fprintf( fp, "#" );

    for( i=0; i<(hamiltonian.N_nodes +1); i++ ){

      fprintf( fp, " " );

      if( fprintf( fp, INT_FORMAT, hamiltonian.neighbour_index[ i ] ) < 1 );
    
    }        

    fprintf( fp, "\n" );

    /* neighbour_list */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# neighbour_list:\n" );

    fprintf( fp, "#" );

    for( i=0; i<(2 *hamiltonian.N_links); i++ ){

      fprintf( fp, " " );

      if( fprintf( fp, INT_FORMAT, hamiltonian.neighbour_list[ i ] ) < 1 );

    }        
    
    fprintf( fp, "\n" );

    /* neighbour_link */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# neighbour_link:\n" );

    fprintf( fp, "#" );

    for( i=0; i<(2 *hamiltonian.N_links); i++ ){

      fprintf( fp, " " );

      if( fprintf( fp, INT_FORMAT, hamiltonian.neighbour_link[ i ] ) < 1 );

    }        
    
    fprintf( fp, "\n" );
    
  }

  /* N_par_node */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_par_node:\n" );

  if( fprintf( fp, "# %d\n", hamiltonian.N_par_node ) < 1 ) info=1;

  /* N_par_link */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_par_link:\n" );

  if( fprintf( fp, "# %d\n", hamiltonian.N_par_link ) < 1 ) info=1;

  /* N_par_extra */
  fprintf( fp, "#------------------------------------------#\n" );
  fprintf( fp, "# N_par_extra:\n" );

  if( fprintf( fp, "# %d\n", hamiltonian.N_par_extra ) < 1 ) info=1;

  if( hamiltonian.N_nodes > 0 ){

    /* par_node */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# par_node:\n" );

    for( i=0; i<hamiltonian.N_nodes; i++ ){

      fprintf( fp, "#" );

      for( j=0; j<hamiltonian.N_par_node; j++ ){

        fprintf( fp, " " );

        if( fprintf( fp, DOUBLE_FORMAT, hamiltonian.par_node[i][j] ) < 1 );
    
      }        

      fprintf( fp, "\n" );

    }

  }        

  if( hamiltonian.N_links > 0 ){

    /* par_link */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# par_link:\n" );

    for( i=0; i<hamiltonian.N_links; i++ ){

      fprintf( fp, "#" );

      for( j=0; j<hamiltonian.N_par_link; j++ ){

        fprintf( fp, " " );

        if( fprintf( fp, DOUBLE_FORMAT, hamiltonian.par_link[i][j] ) < 1 );
    
      }        

      fprintf( fp, "\n" );

    }

  }

  if( hamiltonian.N_par_extra > 0 ){

    /* par_extra */
    fprintf( fp, "#------------------------------------------#\n" );
    fprintf( fp, "# par_extra:\n" );

    for( i=0; i<hamiltonian.N_par_extra; i++ ){

      if( fprintf( fp, "# "DOUBLE_FORMAT"\n", hamiltonian.par_extra[i] ) < 1 );
    
    }        

  }

  fflush( fp );


  return info;

}

//------------------------------------------
