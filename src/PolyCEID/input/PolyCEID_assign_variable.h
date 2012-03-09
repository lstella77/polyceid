
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

/* FILE: PolyCEID_assign_variable.h
 *
 * It contains function for parsing
 * The input file is made by string fields
 * followed by variables
 *
 */

#ifndef _PolyCEID_ASSIGN_VARIABLE_
#define _PolyCEID_ASSIGN_VARIABLE_

#include "../input/PolyCEID_list.h"
#include "../structures/PolyCEID_structures_atoms.h"
#include "../utils/PolyCEID_error.h"
#include <float.h>
#include <limits.h>

/*****************************
 **********FUNCTIONS**********
 *****************************/

/* FUNCTION: _assign_int_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 *            int     default value
 *
 * RETURN:    int
 *
 * search for the int variable in the list and assign 
 * its value
 *
 */
int _assign_int_variable( list_p, char*, int, int*, int* );

#define ASSING_INT_VARIABLE( lp, name, n, ip, idp )  FUNCTION_CHECK_NEW( _assign_int_variable( (lp), (name), (n), (ip), (idp) ) )

/* FUNCTION: _assign_double_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 *            double     default value
 *
 * RETURN:    int
 *
 * search for the double variable in the list and assign 
 * its value
 *
 */
int _assign_double_variable( list_p, char*, int, double*, double* );

#define ASSING_DOUBLE_VARIABLE( lp, name, n, dp, ddp )  FUNCTION_CHECK_NEW( _assign_double_variable( (lp), (name), (n), (dp), (ddp) ) )

/* FUNCTION: _assign_string_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 *            char*   default value
 *
 * RETURN:    char*
 *
 * search for the string variable in the list and assign 
 * its value
 *
 */
int _assign_string_variable( list_p, char*, char*, char* );

#define ASSING_STRING_VARIABLE( lp, name, sp, sdp )  FUNCTION_CHECK_NEW( _assign_string_variable( (lp), (name), (sp), (sdp) ) )

/* FUNCTION: _assign_stringcat_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 *            char*   default value
 *
 * RETURN:    char*
 *
 * search for the string variable in the list and assign 
 * its value
 *
 */
int _assign_stringcat_variable( list_p, char*, char*, char* );

#define ASSING_STRINGCAT_VARIABLE( lp, name, sp, sdp )  FUNCTION_CHECK_NEW( _assign_stringcat_variable( (lp), (name), (sp), (sdp) ) )

/* FUNCTION: _assign_flag_variable
 *
 * ARGUMENTS: list_p  list_p
 *
 *            char*   name of the int variable
 *
 * RETURN:    unsigned short int
 *
 * search for the flag variable in the list and assign 
 * its value
 *
 * WARNING: not to be preprocessed!
 *
 */
unsigned short int _assign_flag_variable( list_p, char* );

#define ASSING_FLAG_VARIABLE( lp, name )  _assign_flag_variable( (lp), (name) )

#endif /* _PolyCEID_ASSIGN_VARIABLE_ */
