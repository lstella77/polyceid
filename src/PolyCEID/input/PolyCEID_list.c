
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
/* FILE: PolyCEID_list.c
 *
 * contains basis list utilities 
 *
 */

// headers
#include "PolyCEID_list.h"

/*****************************
 **********FUNCTIONS**********
 *****************************/

/* FUNCTION: _list_read
 *
 * ARGUMENTS: FILE*    stream
 *
 *            list_p  list_p
 *
 * reads a list from the stream
 *
 */
int _list_read( FILE* fp, list_p list_p ){

  /* dummies */
  field                      dummy_field;
  entry_p                    dummy_entry_p=NULL;
  FILE*                      fp_tmp;
  static int                 level=0;
  int                        level_ini;
  int                        info=0;


  // pointer check
  if( list_p ){

    // initialisation
    if( !level ){

      LIST_INIT( list_p );

    }        

    // set initial level
    level_ini = level;
  
  
    // main loop
    while( FIELD_NEW_READ( fp, &dummy_field) ){

      S_ECHO( "field read", dummy_field.name );

      if( IS_FIELD_NEW( FILE_SEP, &dummy_field ) ){

        // read the file name
        FIELD_NEW_READ( fp, &dummy_field); 

        // add the entry to the list
        ENTRY_ADD( dummy_field.name, list_p );
        
        // update dummy_entry
        // if first entry
        if( 1 == list_p->N_entries ){
		
          dummy_entry_p = list_p->first_entry_p;

        }
        else{

          dummy_entry_p = dummy_entry_p->next_entry_p;

        } /* end if first entry */     
	
	// read from file
        S_ECHO( "reading from file", dummy_field.name );

        fp_tmp = fopen( dummy_field.name, "rt" );

        if( !fp ){

          fprintf( PolyCEID_STDERR, "ERROR: file \"%s\" not there\n", dummy_field.name );

          fflush( PolyCEID_STDERR ); 

          info=1;

	  break;

	}	
        else{

          // increase level --- a sublist starts
          level++;

          // add a sublist
          SUBLIST_ADD( dummy_entry_p );
        
          // recurrency
          LIST_READ( fp_tmp, dummy_entry_p->sublist_p );

          // close file
          S_ECHO( "file closing", dummy_field.name );
  
          fclose( fp_tmp );

          // decrease leve --- a sublist is closed
          level--;
	
	}  

      }        
      else if( IS_FIELD_NEW( OPEN_SEP, &dummy_field ) ){

        // increase level --- a sublist starts
        level++;
        D_ECHO( "open", level );

        // add a sublist
        SUBLIST_ADD( dummy_entry_p );
        
        // recurrency
        LIST_READ( fp, dummy_entry_p->sublist_p );
      
      }
      else if( IS_FIELD_NEW( CLOSE_SEP, &dummy_field ) ){

        // level consistency check --- the final level must be equal to the initial one
        if( level != level_ini ){

          fprintf( PolyCEID_STDERR, "ERROR: level inconsistency in reading list '%s': level %d found instead of %d\n", list_p->title.name, level, level_ini );

          fflush( PolyCEID_STDERR );

          info=1;

          break;
	  
        }

        D_ECHO( "close", level );

        // decrease leve --- a sublist is closed
        level--;
	
        // level consistency check --- the fianl level must be positive
        if( level < 0 ){

          fprintf( PolyCEID_STDERR, "ERROR: negative level %d found instead of a positive one\n", level );

          fflush( PolyCEID_STDERR );

          info=1;

          break;

        }

        // avoid infinite loop...
        break;

      }
      else{

        // add the entry to the list
        ENTRY_ADD( dummy_field.name, list_p );
        
        // update dummy_entry
        // if first entry
        if( 1 == list_p->N_entries ){
		
          dummy_entry_p = list_p->first_entry_p;

        }
        else{

          dummy_entry_p = dummy_entry_p->next_entry_p;

        } /* end if first entry */
	
      }        

    } /* end main loop */       
                               
    // level consistency check --- at the end of the file all the brackets must be closed
    if( feof( fp ) && level != level_ini ){

      fprintf( PolyCEID_STDERR, "ERROR: level %d instead of 0 found at the end of file '%p'\n", level, fp );

      fflush( PolyCEID_STDERR );

      info=1;

    }        
  
  }
  else{

    fprintf( PolyCEID_STDERR, "ERROR: list_p %p is pointing to NULL\n", list_p );

    fflush( PolyCEID_STDERR );

    info=1;

    
  } /* end poiner check */        


  return info;

}


/* FUNCTION: _list_print
 *
 * ARGUMENTS: FILE*    stream
 *
 *            list_p  list_p
 *
 * prints a list to the stream
 *
 */
int _list_print( FILE* fp, list_p list_p ){

  /* dummies */
  entry_p     dummy_entry_p;
  int         info=0;


  // check list
  if( !list_p ){

    fprintf( PolyCEID_STDERR, "ERROR: no list to print\n" );      
  
    fflush( PolyCEID_STDERR );

    info=1;
  
  }
  else{

    if( !list_p->level ){

      FIELD_NEW_PRINT( fp, &list_p->title );

      fprintf( fp, "{\n" );

    }  

    if( list_p->parent_entry_p ){

      S_ECHO( "parent entry", list_p->parent_entry_p->field.name );

    }   

    if( list_p->first_entry_p ){

      S_ECHO( "first entry", list_p->first_entry_p->field.name );

    }   

    // if there is at least one entry
    if( list_p->first_entry_p ){

      // set dummy_entry_p
      dummy_entry_p = list_p->first_entry_p;

      // main loop
      while( dummy_entry_p ){

        // entry print
        ENTRY_PRINT( fp,  dummy_entry_p ); 

        // update dummy_entry_p;
        dummy_entry_p = dummy_entry_p->next_entry_p;
      
      } /* end main loop */

    } /* end if there is at least one entry */
    
    if( !list_p->level ){

      fprintf( fp, "}\n" );

    }  

  } /* end check list */

  
  return info;

}

/* FUNCTION: _list_cut
 *
 * ARGUMENTS: list_p   list_p
 *
 * cuts list_p (it might be a sublist)
 *
 */
int _list_cut( list_p list_p ){

  /* dummies */
  int info=0;


  // pointer check
  if( !list_p ){

    fprintf( PolyCEID_STDERR, "ERROR: list_p %p points to NULL\n", list_p );

    fflush( PolyCEID_STDERR );

    info=1;

  }        
  else{

    S_ECHO( "cutting list", list_p->title.name );

    // if thre is at least one entry
    if( list_p->first_entry_p ){

      // cut all the entries from the second on
      while( list_p->first_entry_p->next_entry_p ){

        ENTRY_CUT( list_p->first_entry_p->next_entry_p );
        
      }        

      // cut the first entry and
      // prune the parent list, if not the root
      S_ECHO( "cutting first entry", list_p->first_entry_p->field.name );
      
      // check N_entries consistency
      if( list_p->N_entries != 1 ){

        fprintf( PolyCEID_STDERR, "ERROR: N_entries %d is not 1 at the end of list %s cutting\n", list_p->N_entries, list_p->title.name );

	fflush( PolyCEID_STDERR );

	info=1;

      }
      else{

        ENTRY_CUT( list_p->first_entry_p );
        
      }        

    }

  }


  return info;

}

/* FUNCTION: _list_search
 *
 * ARGUMENTS: char*  name
 *
 *            list_p   l_p
 *
 *  searches for the entry containing the field 'name'
 *
 *  and gives a pointer to it
 *
 * NOTE: the function must not be preprocessed
 *
 */
entry_p  _list_search( list_p l_p, char* name ){

  /* dummies */
  entry_p         e_p=NULL;
  entry_p         dummy_entry_p=NULL;


  // pointer check
  if( !l_p ){

    fprintf( PolyCEID_STDERR, "ERROR: no list %s to be searched\n", name );

    fflush( PolyCEID_STDERR );

  }
  else{

    // if there's at least one entry
    if( l_p->first_entry_p ){

      S_ECHO( "searching for", name );

      // pointer update
      dummy_entry_p = l_p->first_entry_p;
      
      // main loop 
      while( dummy_entry_p ){

        S_ECHO( "searching in field", dummy_entry_p->field.name );

        // check field name
	if( IS_FIELD_NEW( name, &dummy_entry_p->field ) ){

           S_ECHO( "field found", dummy_entry_p->field.name );

           // pointer update
           e_p = dummy_entry_p;

           break;

	}	
        else{
		
          // if there's a sublist
	  if( dummy_entry_p->sublist_p ){

             S_ECHO( "searching sublist", dummy_entry_p->sublist_p->title.name );

             // recurrency
	     e_p = LIST_SEARCH( dummy_entry_p->sublist_p, name );

             // check if name has been found
	     if( e_p ){
		     
               S_ECHO( "field found [sublist]", e_p->field.name );

	       break;

	     }  
	     
	  }	

          // update pointer
          dummy_entry_p = dummy_entry_p->next_entry_p;

        } /* end check field name*/

      } /* end main loop */       

    }
    
  }
  
  
  return e_p;

}  

/*****************************
 *********UTLITIES************
 *****************************/

/* FUNCTION: _list_init
 *
 * ARGUMENTS: list_p
 *
 * initialise an existing list
 *
 */
int   _list_init( list_p list_p ){

  /* dummies */	
  int info=0;


  // if the list exists
  if( list_p ){

    FIELD_NEW_ASSIGN( "root", &list_p->title );

    S_ECHO( "Initialize the list", list_p->title.name );

    list_p->level=0;
 
    list_p->N_entries=0;
 
    list_p->first_entry_p=NULL;
 
    list_p->parent_entry_p=NULL;
 
  }
  else{

    info=1;

  } /* end if the list exists */        
 
 
  return info;

}

/* FUNCTION:  _sublist_add
 *
 * ARGUMENTS: entry_p entry_p
 *
 * add a sublist the list
 *
 */
int _sublist_add( entry_p entry_p ){

  /* dummies */
  int      info=0;


  // if the entry exists
  if( entry_p ){
    
    // allocate sublist
    entry_p->sublist_p = (list_p) calloc( ( size_t ) 1, sizeof( list ) );

    // pointer check
    if( !entry_p->sublist_p ){

      fprintf( PolyCEID_STDERR, "ERROR: failed sublist allocation of entry '%s'\n", entry_p->field.name );        
       
      fflush( PolyCEID_STDERR );
       
      info=1;

    }        
    else{
	      
      FIELD_NEW_ASSIGN( entry_p->field.name, &entry_p->sublist_p->title ); // sublist's name is the same of the parent entry    
      entry_p->sublist_p->level         = entry_p->level +1;

      entry_p->sublist_p->N_entries     = 0;

      entry_p->sublist_p->first_entry_p = NULL;

      entry_p->sublist_p->parent_entry_p = entry_p;

    } /* end if pointer check */

  } /* end if the list exist */        

  return info;

}  

/* FUNCTION: _entry_print
 *
 * ARGUMENTS: FILE*    stream
 *
 *            entry_p  entry_p
 *
 * prints an entry to the stream
 *
 */
int _entry_print( FILE* fp, entry_p entry_p ){

  /* dummies */
  int level=0;
  int i;
  int info=0;


  // check list
  if( !entry_p ){

    fprintf( PolyCEID_STDERR, "ERROR: no entry to print\n" );      
  
    fflush( PolyCEID_STDERR );

    info=1;
  
  }
  else{

    // set level
    level = entry_p->level;  

    // indentation
    for( i=0; i< level+1; i++ ){

       fprintf( fp, "  " ); 

    }        

    FIELD_NEW_PRINT( fp, &entry_p->field );
    
    if( entry_p->parent_list_p ){

      S_ECHO( "parent list", entry_p->parent_list_p->title.name );

    }	

    if( entry_p->previous_entry_p ){

      S_ECHO( "previous entry", entry_p->previous_entry_p->field.name );

    }	

    if( entry_p->next_entry_p ){

      S_ECHO( "next entry", entry_p->next_entry_p->field.name );

    }	

    if( entry_p->sublist_p ){

      S_ECHO( "sublist", entry_p->sublist_p->title.name );

    }	

    // check for a sublist
    if( entry_p->sublist_p ){

      // indentation
      for( i=0; i< level+1; i++ ){

        fprintf( fp, "  " ); 

      }        

      fprintf( fp, "{\n" ); 

      // recurrency 
      LIST_PRINT( fp, entry_p->sublist_p );

      // indentation
      for( i=0; i< level+1; i++ ){

        fprintf( fp, "  " ); 

      }        

      fprintf( fp, "}\n" ); 

    } /* end check for a sublist */

  }


  return info;

}

/* FUNCTION:  _entry_add
 *
 * ARGUMENTS: char*   name
 *
 *            list_p  list_p
 *
 * add an entry the list, name --> entry_p->field.name
 *
 */
int _entry_add( char* name, list_p list_p ){

  /* dummies */
  entry_p  dummy_entry_p;
  int      i_list;
  int      info=0;


  // if the list exists
  if( list_p ){
    
    // New entry allocation

    // if first entry
    if( !list_p->N_entries ){

      S_ECHO( "Adding the first entry", list_p->title.name );

      // allocate new entry
      list_p->first_entry_p = (entry_p) calloc( ( size_t ) 1, sizeof( entry ) );
            
      // pointer check
      if( !list_p->first_entry_p ){

        fprintf( PolyCEID_STDERR, "ERROR: failed allocation of the first entry of the list '%s'\n", list_p->title.name );        
        fflush( PolyCEID_STDERR );
       
        info=1;

      } /* pointer check*/        
	     
      // update pointer
      dummy_entry_p = list_p->first_entry_p;

      dummy_entry_p->previous_entry_p = NULL; // is the first entry

    }
    else{
    
      dummy_entry_p=list_p->first_entry_p;
    
      D_ECHO( "adding an entry", list_p->N_entries );
    
      // find the last entry of the list
      for( i_list=0; i_list<list_p->N_entries -1; i_list++ ){

        dummy_entry_p=dummy_entry_p->next_entry_p;

        // pointer check
        if( !dummy_entry_p ){

          fprintf( PolyCEID_STDERR, "ERROR: entry %d of list '%s' points to NULL\n", i_list, list_p->title.name );        
       
          fflush( PolyCEID_STDERR );
       
          info=1;

          break;

        } /* end if pointer check */        

      } /* end for i_list */

      if( !info ){
    
        S_ECHO( "Adding an entry", list_p->title.name );

        // allocate new entry
        dummy_entry_p->next_entry_p = (entry_p) calloc( ( size_t ) 1, sizeof( entry ) );
            
        // pointer check
        if( !dummy_entry_p->next_entry_p ){

          fprintf( PolyCEID_STDERR, "ERROR: failed allocation of the %dth of the list '%s'\n", i_list, list_p->title.name );        
          fflush( PolyCEID_STDERR );
       
          info=1;

        } /* pointer check*/        
	
	// setting the previou entry
        dummy_entry_p->next_entry_p->previous_entry_p = dummy_entry_p;

        // update pointer
        dummy_entry_p = dummy_entry_p->next_entry_p;

      }

    } /* end first entry*/

    if( !info ){
	      
        FIELD_NEW_ASSIGN( name, &dummy_entry_p->field );

        dummy_entry_p->level         = list_p->level; // the level is inherited

        dummy_entry_p->next_entry_p  = NULL;

	dummy_entry_p->sublist_p     = NULL;

	dummy_entry_p->parent_list_p = list_p;

        // update N_entries
        list_p->N_entries++;
    
    }
    
  }
  else{

    info=1;

  } /* end if the list exist */        


  return info;

}  

/* FUNCTION: _entry_cut
 *
 * ARGUMENTS: entry_p   entry_p
 *
 * cuts entry_p from the list
 *
 */
int _entry_cut( entry_p e_p ){

  /* dummies */
  list_p   parent_list_p;
  int      info=0;


  // pointer check
  if( !e_p ){

    fprintf( PolyCEID_STDERR, "ERROR: entry %p points to NULL\n", e_p );        
       
    fflush( PolyCEID_STDERR );
       
    info=1;

  }
  else{

    S_ECHO( "Cutting entry", e_p->field.name );

    // set parent list
    if( e_p->parent_list_p ){

      parent_list_p = e_p->parent_list_p;

      // set previous entry
      if( !e_p->previous_entry_p ){

        S_ECHO( "This is the first entry of the list", e_p->field.name );

        // set the new first entry of the list
        if( !e_p->next_entry_p ){

          S_ECHO( "This is also the last entry of the list", e_p->field.name );

          parent_list_p->first_entry_p = NULL; 

        }  
        else{

          e_p->next_entry_p->previous_entry_p = NULL;

          parent_list_p->first_entry_p = e_p->next_entry_p;

        }

      }
      else{

        // set the new next entry of the list
        if( !e_p->next_entry_p ){

          S_ECHO( "This is the last entry of the list", e_p->field.name );

          e_p->previous_entry_p->next_entry_p = NULL; 

        }  
        else{ 

          e_p->next_entry_p->previous_entry_p = e_p->previous_entry_p;

          e_p->previous_entry_p->next_entry_p = e_p->next_entry_p; 

        }

      }

      // reset N_entries
      parent_list_p->N_entries--;
    
      // if there is a sublist
      if( e_p->sublist_p ){
	      
        // cut all the entries of the sublist
        LIST_CUT( e_p->sublist_p );

      }

      // free allocated memory
      S_ECHO( "deallocating entry", e_p->field.name );

      free( e_p );

      // if the parent list has no remaining entry, 
      // cut the parent list
      if( !parent_list_p->N_entries && parent_list_p->level ){

        parent_list_p->parent_entry_p->sublist_p=NULL;

	free( parent_list_p );

      }

    } 
    else{
    
      fprintf( PolyCEID_STDERR, "ERROR: entry %p points has no parent list\n", e_p );        
       
      fflush( PolyCEID_STDERR );
       
      info=1;

    }
    /* end set parent list */

  }        
  /* end pointer check */


  return info;

}


