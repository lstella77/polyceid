
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

#ifndef _PolyCEID_LIST_
#define _PolyCEID_LIST_

// headers
#include "../input/PolyCEID_field.h"

/*****************************
 **********STRUCTURES*********
 *****************************/

struct  _PolyCEID_list
{

  field    title;

  int      level;

  int      N_entries;

  struct  _PolyCEID_entry*  first_entry_p;

  struct  _PolyCEID_entry*  parent_entry_p; // new

};

typedef  struct  _PolyCEID_list    list;

typedef  struct  _PolyCEID_list*   list_p;

typedef  struct  _PolyCEID_list**  list_pp;

struct  _PolyCEID_entry
{

  field    field;

  int      level;

  struct  _PolyCEID_entry*  next_entry_p;

  struct  _PolyCEID_list*   sublist_p;

  struct  _PolyCEID_entry*  previous_entry_p; // new

  struct  _PolyCEID_list*   parent_list_p;    // new

};

typedef  struct  _PolyCEID_entry    entry;

typedef  struct  _PolyCEID_entry*   entry_p;

typedef  struct  _PolyCEID_entry**  entry_pp;

/*****************************
 **********CONSTANTS**********
 *****************************/

#define OPEN_SEP   "{"

#define CLOSE_SEP  "}"

#define FILE_SEP   "FILE"

/*****************************
 **********FUNCTIONS**********
 *****************************/

/* FUNCTION: _list_read
 *
 * ARGUMENTS: FILE*   stream
 *
 *            list_p  list_p
 *
 * reads a list from the stream
 *
 */
int _list_read( FILE* fp, list_p lits_p );

#define LIST_READ( fp, lp )  FUNCTION_CHECK_NEW( _list_read( (fp), (lp) ) )

/* FUNCTION: _list_print
 *
 * ARGUMENTS: FILE*   stream
 *
 *            list_p  list_p
 *
 * prints a list to the stream
 *
 */
int _list_print( FILE* fp, list_p lits_p );

#define LIST_PRINT( fp, lp )  FUNCTION_CHECK_NEW( _list_print( (fp), (lp) ) )

/* FUNCTION: _list_cut
 *
 * ARGUMENTS: list_p   list_p
 *
 * cuts list_p (it might be a sublist)
 *
 */
int _list_cut( list_p l_p );

#define LIST_CUT( lp )   FUNCTION_CHECK_NEW( _list_cut( (lp) ) )

/* FUNCTION: _list_search
 *
 * ARGUMENTS: list_p   l_p
 *
 *            char*  name
 *
 *  searches for the entry containing the field 'name'
 *
 *  and gives a pointer to it
 *
 *  NOTE: the function must not be preprocessed
 *
 */
entry_p  _list_search( list_p l_p, char* name );

#define LIST_SEARCH( lp, name )    _list_search( (lp), (name) )

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
int   _list_init( list_p list_p );

#define LIST_INIT( lp )   FUNCTION_CHECK_NEW( _list_init( (lp) ) )

/* FUNCTION:  _sublist_add
 *
 * ARGUMENTS: entry_p   entry_p
 *
 * add a sublist the list
 *
 */
int _sublist_add( entry_p entry_p );

#define SUBLIST_ADD( ep )   FUNCTION_CHECK_NEW( _sublist_add( (ep) ) )

/* FUNCTION: _entry_print
 *
 * ARGUMENTS: FILE*    stream
 *
 *            entry_p  entry_p
 *
 * prints an entry to the stream
 *
 */
int _entry_print( FILE* fp, entry_p entry_p );

#define ENTRY_PRINT( fp, ep )  FUNCTION_CHECK_NEW( _entry_print( (fp), (ep) ) )  

/* FUNCTION:  _entry_add
 *
 * ARGUMENTS: char*   name
 *
 *            list_p  list_p
 *
 * add an entry the list, name --> entry_p->field.name
 *
 */
int _entry_add( char* name, list_p list_p );

#define ENTRY_ADD( name, lp )   FUNCTION_CHECK_NEW( _entry_add( (name), (lp) ) )

/* FUNCTION: _entry_cut
 *
 * ARGUMENTS: entry_p   entry_p
 *
 * cuts entry_p from the list
 *
 */
int _entry_cut( entry_p entry_p );
		
#define ENTRY_CUT( ep )    FUNCTION_CHECK_NEW( _entry_cut( (ep) ) )

#endif /* _PolyCEID_LIST */
