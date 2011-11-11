
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
#include "PolyCEID_centre_of_mass.h"



/*********************
  FUNCTIONS & MACROS
*********************/


/* compute_centre_of_mass */

int compute_centre_of_mass( const constants constants, atoms_p atoms_p ){

  /* constants */
  int   sdim;
  int   N_atoms;
  /* state */
  rvector_p  positions_p;
  rvector_p  centre_of_mass_p;
  rvector_p  masses_p;
  double     mass_tot;
  /* dummies */
  int i, comp;
  int info=0;


  N_atoms          =  constants.N_atoms;
  sdim             =  constants.spacial_dimension; 
  positions_p      = &(atoms_p->positions);
  centre_of_mass_p = &(atoms_p->centre_of_mass);
  masses_p         = &(atoms_p->masses);
  mass_tot         = atoms_p->mass_tot;

  /* set to zero */
  if( RVECTOR_ZERO( *centre_of_mass_p ) ) info=1;

  for( i=0; i<N_atoms; i++ ){

    for( comp=0; comp<sdim; comp++ ){

      centre_of_mass_p->rvector[ comp ] += masses_p->rvector[ comp +i *sdim ] *positions_p->rvector[ comp +i *sdim ];

    }
    
  }

  for( comp=0; comp<sdim; comp++ ){

    centre_of_mass_p->rvector[ comp ] /= mass_tot;

    // BUGFIX: is it really needed
    if( constants.flag_periodic_boundary_condition ){

      centre_of_mass_p->rvector[ comp ] = PBC_DIFF( centre_of_mass_p->rvector[ comp ], constants.cell_dim.rvector[ comp ] );

    }

  }


#ifdef __DEBUG_PLUS__

  fprintf( stdout, "centre_of_mass\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *centre_of_mass_p  ) ) info=1;
  fprintf( stdout, "\n" );
  
#endif /* __DEBUG_PLUS__ */

  
  return info;

}

//------------------------------------------

/* centre_of_mass_transform */

//BUGFIX: this was originally TMP
int centre_of_mass_transform( const constants constants, atoms_p atoms_p ){

  /* constants */
  int        N_atoms;
  /* state */
  rvector_p  positions_p;
  /* dummies */
  int        sdim;
  int        i ,j;
  int        index;
  int        info=0;


  N_atoms     =  constants.N_atoms;
  sdim        =  constants.spacial_dimension; 
  positions_p = &(atoms_p->positions);


#ifdef __DEBUG__

  fprintf( stdout, "DOING: centre_of_mass_transform\n" );

#endif /* __DEBUG__ */


  /* compute centre of mass */
  if( COMPUTE_CENTRE_OF_MASS( constants, *atoms_p ) ) info=1;

  
#ifdef __DEBUG__

  fprintf( stdout, "centroid\n" );
  if( RVECTOR_PRINT_PLUS( stdout, atoms_p->centre_of_mass ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


  for( i=0; i<N_atoms; i++ ){
    
    for( j=0; j<sdim; j++ ){
      
      index = j +i *sdim;

      positions_p->rvector[ index ] -= atoms_p->centre_of_mass.rvector[ j ];

      // BUGFIX: is it really needed
      if( constants.flag_periodic_boundary_condition ){

        positions_p->rvector[  index ] = PBC_DIFF( positions_p->rvector[ index ], constants.cell_dim.rvector[ j ] );

      }
      
    }
    
  }

#ifdef __DEBUG__

  fprintf( stdout, "transformed positions\n" );
  if( RVECTOR_PRINT_PLUS( stdout, *positions_p ) ) info=1;
  fprintf( stdout, "\n" );

#endif /* __DEBUG__ */


  return info;

}

//------------------------------------------


