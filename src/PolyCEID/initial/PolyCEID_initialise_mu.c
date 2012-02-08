
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
#include "PolyCEID_initialise_mu.h"


/*********************
  FUNCTIONS & MACROS
*********************/

int initialise_mu( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  unsigned short int  flag_normal_mode_expansion;
  /* dummies */
  int  info=0;


  flag_normal_mode_expansion  =  constants.flag_normal_mode_expansion;


  /* mu_def initialisation */
  if( MU_UPDATE( constants, *state_p, *config_p ) ) info=1;


  /* final transformation/copy */
  if( flag_normal_mode_expansion ){

    if( TRANSFORM_MU( constants, *state_p, *config_p ) ) info=1;

  }


  return info;

}

//------------------------------------------
