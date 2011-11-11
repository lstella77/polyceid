
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

#ifndef __PolyCEID_CI_TABLE__
#define __PolyCEID_CI_TABLE__


#include <stdio.h>
#include "../utils/my_error.h"
#include "../algebra/PolyCEID_matrix.h"
#include "../algebra/PolyCEID_combinatorics.h"
#include "../structures/PolyCEID_structures_constants.h"
#include "../initial/PolyCEID_indexing.h"
#include "../hamiltonians/PolyCEID_hamiltonian_single.h"
#include "../utils/PolyCEID_sorting.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int CI_table_no_spin( constants_p );

#define CI_TABLE_NO_SPIN( c )  FUNCTION_CHECK( CI_table_no_spin( &(c) ), CI_table_no_spin )

//------------------------------------------

int CI_table_spin( constants_p );

#define CI_TABLE_SPIN( c )  FUNCTION_CHECK( CI_table_spin( &(c) ), CI_table_spin )

//------------------------------------------

/* utilities */

//------------------------------------------

int CI_table_choose( const int, const int, const int, const int, imatrix_p );

#define CI_TABLE_CHOOSE( i, j, l, n, imat )  FUNCTION_CHECK( CI_table_choose( (i), (j), (l), (n), &(imat) ), CI_table_choose )

//------------------------------------------

int purge_spin_flip( const imatrix, imatrix_p );

#define PURGE_SPIN_FLIP( imat1, imat2 ) FUNCTION_CHECK( purge_spin_flip( (imat1), &(imat2) ), purge_spin_flip )

//------------------------------------------

int energy_pruning( constants_p );

#define ENERGY_PRUNING( c ) FUNCTION_CHECK( energy_pruning( &(c) ), energy_pruning )

//------------------------------------------

#endif /* __PolyCEID_CI_TABLE__ */
