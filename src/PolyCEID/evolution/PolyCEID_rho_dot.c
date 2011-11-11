
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
#include "PolyCEID_rho_dot.h"


#define EPS_LOC   EPS


/*********************
  FUNCTIONS & MACROS
*********************/

/* rho_dot_update */

int rho_dot_update( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  unsigned short int  flag_no_Ehrenfest_frame;
  int             N_coor_red;
  int             max_rho_index;
  ivector         relevant_modes;
  rvector         masses_aux;
  /* state */
  matrix_array_p  rho_p;
  matrix_array_p  rho_dot_p;
  matrix_p        H_matrix_p;
  matrix_array_p  delta_F_matrix_p;
  matrix_array_p  K_matrix_p;
  matrix_p        dummy_matrix1_p;
  matrix_p        dummy_matrix2_p;
  matrix_p        dummy_matrix3_p;
  matrix_p        dummy_matrix4_p;
  rvector_p       ar_p;
  rvector_p       ap_p;
  /* dummies */
  int             i_rho;
  int             i_mode;
  int             j_mode;
  int             i_coor1;
  int             i_coor2;
  int             index;
  int             K_matrix_index;
  int             info=0;
  complex         dummy;
#ifdef __DEBUG_PLUS__
  double          norm;
#endif /* __DEBUG_PLUS__ */


#ifdef __DEBUG_PLUS__
  fprintf( stdout, "DOING: rho_dot_update\n" );
#endif /* __DEBUG_PLUS__ */


  flag_no_Ehrenfest_frame =  constants.flag_no_Ehrenfest_frame;
  N_coor_red            =  constants.N_coor_red;
  max_rho_index         =  constants.max_rho_index;
  relevant_modes        =  constants.relevant_modes;
  masses_aux            =  config_p->atoms.masses_aux;

  rho_p                 = &(config_p->electrons.rho);
  rho_dot_p             = &(state_p->rho_dot);
  H_matrix_p            = &(config_p->electrons.H_matrix);
  delta_F_matrix_p      = &(config_p->electrons.delta_F_matrix);
  K_matrix_p            = &(config_p->electrons.K_matrix);
  dummy_matrix1_p       = &(state_p->dummy_matrix1);
  dummy_matrix2_p       = &(state_p->dummy_matrix2);
  dummy_matrix3_p       = &(state_p->dummy_matrix3);
  dummy_matrix4_p       = &(state_p->dummy_matrix4);
  ar_p               = &(state_p->phonons.ar);
  ap_p               = &(state_p->phonons.ap);


  /* main loop */
  for( i_rho=0; i_rho<max_rho_index; i_rho++ ){


#ifdef __DEBUG_PLUS__
    fprintf( stdout, "->i_rho = %d\n", i_rho );
#endif /* __DEBUG_PLUS__ */


    /* global */
    if( RHO_INDEX_INVERSE( i_rho, dummy_rho_index1, dummy_rho_index2 ) ) info=1;

    /*
      fprintf( stdout, "  dummy_rho_indices\n" );
      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index1 ) ) info=1;
      if( IVECTOR_PRINT_PLUS( stdout, dummy_rho_index2 ) ) info=1;
    */

    /* set the rho_dot component to zero */
    if( MATRIX_ZERO( rho_dot_p->array[ i_rho ] ) ) info=1;



    /*----------------------
      First part
      ----------------------*/


    //    fprintf( stdout, "  ---------- First part -----------\n" );


    /* set dummy matrix to zero */
    if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;

    /* coor loop */
    for( i_mode=0; i_mode<N_coor_red; i_mode++ ){


      //      fprintf( stdout, "-->i_mode = %d\n", i_mode );

      /* set dummy matrix to zero */
      if( MATRIX_ZERO( *dummy_matrix2_p ) ) info=1;

      i_coor1 = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate


      //      fprintf( stdout, "   i_coor1 = %d, i_atom = %d \n", i_coor1, i_atom );


      // term 1
      index = RHO_INDEX_CHANGE1_PLUS2_FIRST( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix3_p  ) ) info=1;
	
	if( MATRIX_SUM( *dummy_matrix2_p, *dummy_matrix3_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      }

      // term 2
      index = i_rho; // for compatibility


      //      fprintf( stdout, "   index = %d \n", index );


      dummy = CMPLX( 2.0e0 *( dummy_rho_index1.ivector[ i_mode ] ) +1.0e0 );

      /*
	fprintf( stdout, "dummy\n" );
	if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	fprintf( stdout, "\n" );
      */

      if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix3_p  ) ) info=1;

      if( MATRIX_DIF( *dummy_matrix2_p, *dummy_matrix3_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      // term 3
      index = RHO_INDEX_CHANGE1_MINUS2_FIRST( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] ) *( dummy_rho_index1.ivector[ i_mode ] -1.0e0 ) ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/
	
	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix3_p  ) ) info=1;
	
	if( MATRIX_SUM( *dummy_matrix2_p, *dummy_matrix3_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      }

      // term 4
      index = RHO_INDEX_CHANGE1_MINUS2_SECOND( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] ) *( dummy_rho_index2.ivector[ i_mode ] -1.0e0 ) ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix3_p  ) ) info=1;
	
	if( MATRIX_DIF( *dummy_matrix2_p, *dummy_matrix3_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      }

      // term 5
      index = i_rho; // for compatibility


      //      fprintf( stdout, "   index = %d \n", index );


      dummy = CMPLX( 2.0e0 *( dummy_rho_index2.ivector[ i_mode ] ) +1.0e0 );

      /*
	fprintf( stdout, "dummy\n" );
	if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	fprintf( stdout, "\n" );
      */

      if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix3_p  ) ) info=1;

      if( MATRIX_SUM( *dummy_matrix2_p, *dummy_matrix3_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!


      // term 6
      index = RHO_INDEX_CHANGE1_PLUS2_SECOND( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ i_mode ] +1.0e0 ) ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/
	
	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix3_p  ) ) info=1;
	
	if( MATRIX_DIF( *dummy_matrix2_p, *dummy_matrix3_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      }

      // final

      dummy = CMPLX( ( ap_p->rvector[ i_coor1 ] ) *( ap_p->rvector[ i_coor1 ] ) /( masses_aux.rvector[ i_coor1 ] ) );

      if( MATRIX_SCALAR_PRODUCT( *dummy_matrix2_p, dummy, *dummy_matrix2_p  ) ) info=1;  //WARNING: overwriting!

      if( MATRIX_SUM( *dummy_matrix1_p, *dummy_matrix2_p, *dummy_matrix1_p  ) ) info=1;  //WARNING: overwriting!

    } /* i_mode loop */

    dummy = CMPLX_DIVISION( CMPLX( 0.25e0 ), IHBAR );

    if( MATRIX_SCALAR_PRODUCT( *dummy_matrix1_p, dummy, *dummy_matrix1_p ) );  //WARNING: overwriting!

    if( MATRIX_DIF( rho_dot_p->array[ i_rho ], *dummy_matrix1_p, rho_dot_p->array[ i_rho ] ) ) info=1;  //WARNING: overwriting!


#ifdef __DEBUG_PLUS__
    fprintf( stdout, "   --------------\n" );
    norm = MATRIX_NORM( *dummy_matrix1_p );
    fprintf( stdout, "   Term's norm [1] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */




    if( flag_no_Ehrenfest_frame ){

      /*----------------------
        Second part [Ehrenfest]
        ----------------------*/
    

      //    fprintf( stdout, "  ------------ Second part [Ehrenfest] ------------\n" );


      index = i_rho; // for compatibility

    
      //    fprintf( stdout, "   index = %d \n", index );

      /*
        fprintf( stdout, "H_matrix\n" );
        if( MATRIX_PRINT_PLUS( stdout, *H_matrix_p ) ) info=1;
        fprintf( stdout, "\n" );


        fprintf( stdout, "rho\n" );
        if( MATRIX_PRINT_PLUS( stdout, rho_p->array[ index ] ) ) info=1;
        fprintf( stdout, "\n" );
      */
    
      if( MATRIX_COMMUTATOR( *H_matrix_p, rho_p->array[ index ], *dummy_matrix1_p ) ) info=1;
    
      /*
        fprintf( stdout, "dummy_matrix1\n" );
        if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p ) ) info=1;
        fprintf( stdout, "\n" );
      */

      dummy = CMPLX_DIVISION( CMPLX( 1.0e0 ), IHBAR );

      /*
        fprintf( stdout, "dummy\n" );
        if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
        fprintf( stdout, "\n" );
      */ 
    
      if( MATRIX_SCALAR_PRODUCT( *dummy_matrix1_p, dummy, *dummy_matrix1_p ) );  //WARNING: overwriting!

      // final    
      if( MATRIX_SUM( rho_dot_p->array[ i_rho ], *dummy_matrix1_p, rho_dot_p->array[ i_rho ] ) ) info=1;  //WARNING: overwriting!


#ifdef __DEBUG_PLUS__
      fprintf( stdout, "   --------------\n" );
      norm = MATRIX_NORM( *dummy_matrix1_p );
      fprintf( stdout, "   Term's norm [2] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */

    }


    /*----------------------
      Third part
      ----------------------*/


    //    fprintf( stdout, "  -------------- Third part -----------\n" );


    /* set dummy matrix to zero */
    if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;

    /* coor loop */
    for( i_mode=0; i_mode<N_coor_red; i_mode++ ){


      //      fprintf( stdout, "-->i_mode = %d\n", i_mode );

      /* set dummy matrix to zero */
      if( MATRIX_ZERO( *dummy_matrix2_p ) ) info=1;

      i_coor1 = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate


      //      fprintf( stdout, "   i_coor1 = %d \n", i_coor1 );


      // term 1

      /* set dummy matrix to zero */
      if( MATRIX_ZERO( *dummy_matrix3_p ) ) info=1;

      // A
      index = RHO_INDEX_CHANGE1_PLUS1_FIRST( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      }

      // B
      index = RHO_INDEX_CHANGE1_MINUS1_FIRST( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( dummy_rho_index1.ivector[ i_mode ] ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      }

      // 
      if( MATRIX_MATRIX_PRODUCT( delta_F_matrix_p->array[ i_coor1 ] , *dummy_matrix3_p,  *dummy_matrix4_p ) ) info=1;

      if( MATRIX_SUM( *dummy_matrix2_p, *dummy_matrix4_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      // term 2

      /* set dummy matrix to zero */
      if( MATRIX_ZERO( *dummy_matrix3_p ) ) info=1;

      // A
      index = RHO_INDEX_CHANGE1_MINUS1_SECOND( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( dummy_rho_index2.ivector[ i_mode ] ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      }

      // B
      index = RHO_INDEX_CHANGE1_PLUS1_SECOND( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( dummy_rho_index2.ivector[ i_mode ] +1.0e0 ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      }

      //
      if( MATRIX_MATRIX_PRODUCT( *dummy_matrix3_p, delta_F_matrix_p->array[ i_coor1 ] , *dummy_matrix4_p ) ) info=1;

      if( MATRIX_DIF( *dummy_matrix2_p, *dummy_matrix4_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      // final

      dummy = CMPLX( ar_p->rvector[ i_coor1 ] );

      if( MATRIX_SCALAR_PRODUCT( *dummy_matrix2_p, dummy, *dummy_matrix2_p  ) ) info=1;  //WARNING: overwriting!
    
      if( MATRIX_SUM( *dummy_matrix1_p, *dummy_matrix2_p, *dummy_matrix1_p  ) ) info=1;  //WARNING: overwriting!

    } /* i_mode loop */

    dummy = CMPLX_DIVISION( CMPLX( OSQRT2 ), IHBAR );

    if( MATRIX_SCALAR_PRODUCT( *dummy_matrix1_p, dummy, *dummy_matrix1_p ) );  //WARNING: overwriting!

    if( MATRIX_DIF( rho_dot_p->array[ i_rho ], *dummy_matrix1_p, rho_dot_p->array[ i_rho ] ) ) info=1;  //WARNING: overwriting!


#ifdef __DEBUG_PLUS__
    fprintf( stdout, "   --------------\n" );
    norm = MATRIX_NORM( *dummy_matrix1_p );
    fprintf( stdout, "   Term's norm [3] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */
 


    /*----------------------
      Fourth part
      ----------------------*/
    

    //    fprintf( stdout, "  ------------ Fourth part -----------\n" );


    /* set dummy matrix to zero */
    if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;

    /* coor loop */
    for( i_mode=0; i_mode<N_coor_red; i_mode++ ){


      //      fprintf( stdout, "-->i_mode = %d\n", i_mode );

      /* set dummy matrix to zero */
      if( MATRIX_ZERO( *dummy_matrix2_p ) ) info=1;

      i_coor1 = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate

      K_matrix_index = COORDINATE_INDEX( i_coor1, i_coor1 );

      //      fprintf( stdout, "   i_coor1 = %d, K_matrix_index = %d \n", i_coor1, K_matrix_index );


      // term 1

      /* set dummy matrix to zero */
      if( MATRIX_ZERO( *dummy_matrix3_p ) ) info=1;

      // A
      index = RHO_INDEX_CHANGE1_PLUS2_FIRST( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) ) ); //WARNING: already computed...

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      }

      // B
      index = i_rho;


      //      fprintf( stdout, "   index = %d \n", index );


      dummy = CMPLX( 2.0e0 *( dummy_rho_index1.ivector[ i_mode ] ) +1.0e0 );

      /*
	fprintf( stdout, "dummy\n" );
	if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	fprintf( stdout, "\n" );
      */

      if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

      if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      // C
      index = RHO_INDEX_CHANGE1_MINUS2_FIRST( i_rho, i_mode );


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] ) *( dummy_rho_index1.ivector[ i_mode ] -1.0e0 ) ) ); //WARNING: already computed...

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      }

      // 
      if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index ] , *dummy_matrix3_p,  *dummy_matrix4_p ) ) info=1;

      if( MATRIX_SUM( *dummy_matrix2_p, *dummy_matrix4_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      // term 2

      /* set dummy matrix to zero */
      if( MATRIX_ZERO( *dummy_matrix3_p ) ) info=1;

      // A
      index = RHO_INDEX_CHANGE1_MINUS2_SECOND( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] ) *( dummy_rho_index2.ivector[ i_mode ] -1.0e0 ) ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      }

      // B
      index = i_rho;


      //      fprintf( stdout, "   index = %d \n", index );


      dummy = CMPLX( 2.0e0 *( dummy_rho_index2.ivector[ i_mode ] ) +1.0e0 );
      
      /*
	fprintf( stdout, "dummy\n" );
	if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	fprintf( stdout, "\n" );
	*/
      
      if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;
      
      if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      // C
      index = RHO_INDEX_CHANGE1_PLUS2_SECOND( i_rho, i_mode ); 


      //      fprintf( stdout, "   index = %d \n", index );


      if( index < max_rho_index ){

	dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] +2.0e0 ) *( dummy_rho_index2.ivector[ i_mode ] +1.0e0 ) ) );

	/*
	  fprintf( stdout, "dummy\n" );
	  if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	  fprintf( stdout, "\n" );
	*/

	if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

      }

      //
      if( MATRIX_MATRIX_PRODUCT( *dummy_matrix3_p, K_matrix_p->array[ K_matrix_index ] , *dummy_matrix4_p ) ) info=1;

      if( MATRIX_DIF( *dummy_matrix2_p, *dummy_matrix4_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

      // final

      dummy = CMPLX( ( ar_p->rvector[ i_coor1 ] ) *( ar_p->rvector[ i_coor1 ] ) );

      if( MATRIX_SCALAR_PRODUCT( *dummy_matrix2_p, dummy, *dummy_matrix2_p  ) ) info=1;  //WARNING: overwriting!
    
      if( MATRIX_SUM( *dummy_matrix1_p, *dummy_matrix2_p, *dummy_matrix1_p  ) ) info=1;  //WARNING: overwriting!

    } /* i_mode loop */

    dummy = CMPLX_DIVISION( CMPLX( 0.25e0 ), IHBAR );

    if( MATRIX_SCALAR_PRODUCT( *dummy_matrix1_p, dummy, *dummy_matrix1_p ) );  //WARNING: overwriting!

    if( MATRIX_SUM( rho_dot_p->array[ i_rho ], *dummy_matrix1_p, rho_dot_p->array[ i_rho ] ) ) info=1;  //WARNING: overwriting!


#ifdef __DEBUG_PLUS__
    fprintf( stdout, "   --------------\n" );
    norm = MATRIX_NORM( *dummy_matrix1_p );
    fprintf( stdout, "   Term's norm [4] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */


#ifndef __ENERGY_MOD__

    /*----------------------
      Fifth part
      ----------------------*/


    //    fprintf( stdout, "  -------------- Fifth part -------------\n" );


    /* set dummy matrix to zero */
    if( MATRIX_ZERO( *dummy_matrix1_p ) ) info=1;

    /* coor loop */
    for( i_mode=0; i_mode<N_coor_red; i_mode++ ){


      //      fprintf( stdout, "-->i_mode = %d\n", i_mode );
      

      for( j_mode=(i_mode+1); j_mode<N_coor_red; j_mode++ ){ // WARNING: j_mode != i_mode


	//	fprintf( stdout, "-->j_mode = %d\n", i_mode );


	/* set dummy matrix to zero */
	if( MATRIX_ZERO( *dummy_matrix2_p ) ) info=1;
	
	i_coor1 = relevant_modes.ivector[ i_mode ]; // i_mode to cartesian coordinate

	i_coor2 = relevant_modes.ivector[ j_mode ]; // j_mode to cartesian coordinate

	K_matrix_index = COORDINATE_INDEX( i_coor1, i_coor2 );


	//	fprintf( stdout, "   i_coor1 = %d, i_coor2 = %d, K_matrix_index = %d \n", i_coor1, i_coor2, K_matrix_index );


	// term 1

	/* set dummy matrix to zero */
	if( MATRIX_ZERO( *dummy_matrix3_p ) ) info=1;

	// A
	index = RHO_INDEX_CHANGE2_PLUS1_PLUS1_FIRST( i_rho, i_mode, j_mode ); 


	//	fprintf( stdout, "   index = %d \n", index );

	
	if( index < max_rho_index ){
	
	  dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ j_mode ] +1.0e0 ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	  if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

	}

	// B
	index = RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( i_rho, i_mode, j_mode ); 


	//	fprintf( stdout, "   index = %d \n", index );

	
	if( index < max_rho_index ){
	
	  dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index1.ivector[ j_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	  if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

	}

	// C
	index = RHO_INDEX_CHANGE2_PLUS1_MINUS1_FIRST( i_rho, j_mode, i_mode ); //WARNING: note the index inversion


	//	fprintf( stdout, "   index = %d \n", index );

	
	if( index < max_rho_index ){
	
	  dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] ) *( dummy_rho_index1.ivector[ j_mode ] +1.0e0 ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	  if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

	}

	// D
	index = RHO_INDEX_CHANGE2_MINUS1_MINUS1_FIRST( i_rho, i_mode, j_mode ); 



	//	fprintf( stdout, "   index = %d \n", index );

	
	if( index < max_rho_index ){
	
	  dummy = CMPLX( sqrt( ( dummy_rho_index1.ivector[ i_mode ] ) *( dummy_rho_index1.ivector[ j_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	  if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

	}

	/*
	  fprintf( stdout, "dummy_matrix3\n" );
	  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix3_p ) ) info=1;

	  fprintf( stdout, "K_matrix\n" );
	  if( MATRIX_PRINT_PLUS( stdout, K_matrix_p->array[ K_matrix_index ] ) ) info=1;
	*/

	// 
	if( MATRIX_MATRIX_PRODUCT( K_matrix_p->array[ K_matrix_index ], *dummy_matrix3_p, *dummy_matrix4_p ) ) info=1;

	/*
	  fprintf( stdout, "dummy_matrix4\n" );
	  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix4_p ) ) info=1;
	*/

	if( MATRIX_SUM( *dummy_matrix2_p, *dummy_matrix4_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

	// term 2

	/* set dummy matrix to zero */
	if( MATRIX_ZERO( *dummy_matrix3_p ) ) info=1;
      
	// A
	index = RHO_INDEX_CHANGE2_MINUS1_MINUS1_SECOND( i_rho, i_mode, j_mode ); 


	//	fprintf( stdout, "   index = %d \n", index );


	if( index < max_rho_index ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] ) *( dummy_rho_index2.ivector[ j_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	  if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

	}

	// B
	index = RHO_INDEX_CHANGE2_PLUS1_MINUS1_SECOND( i_rho, j_mode, i_mode ); //WARNING: note the index inversion


	//	fprintf( stdout, "   index = %d \n", index );


	if( index < max_rho_index ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	  if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

	}

	// C
	index = RHO_INDEX_CHANGE2_PLUS1_MINUS1_SECOND( i_rho, i_mode, j_mode ); 


	//	fprintf( stdout, "   index = %d \n", index );


	if( index < max_rho_index ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	  if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

	}
	
	// D
	index = RHO_INDEX_CHANGE2_PLUS1_PLUS1_SECOND( i_rho, i_mode, j_mode ); 


	//	fprintf( stdout, "   index = %d \n", index );


	if( index < max_rho_index ){

	  dummy = CMPLX( sqrt( ( dummy_rho_index2.ivector[ i_mode ] +1.0e0 ) *( dummy_rho_index2.ivector[ j_mode ] +1.0e0 ) ) );

	  /*
	    fprintf( stdout, "dummy\n" );
	    if( CMPLX_PRINT_PLUS( stdout, dummy ) ) info=1;
	    fprintf( stdout, "\n" );
	  */

	  if( MATRIX_SCALAR_PRODUCT( rho_p->array[ index ], dummy, *dummy_matrix4_p  ) ) info=1;

	  if( MATRIX_SUM( *dummy_matrix3_p, *dummy_matrix4_p, *dummy_matrix3_p  ) ) info=1; //WARNING: overwriting!

	}

	/*
	  fprintf( stdout, "dummy_matrix3\n" );
	  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix3_p ) ) info=1;
	  
	  fprintf( stdout, "K_matrix\n" );
	  if( MATRIX_PRINT_PLUS( stdout, K_matrix_p->array[ K_matrix_index ] ) ) info=1;
	*/

	//
	if( MATRIX_MATRIX_PRODUCT( *dummy_matrix3_p, K_matrix_p->array[ K_matrix_index ] , *dummy_matrix4_p ) ) info=1;

	/*
	  fprintf( stdout, "dummy_matrix4\n" );
	  if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix4_p ) ) info=1;
	*/

	if( MATRIX_DIF( *dummy_matrix2_p, *dummy_matrix4_p, *dummy_matrix2_p  ) ) info=1; //WARNING: overwriting!

	// final

	dummy = CMPLX( ( ar_p->rvector[ i_coor1 ] ) *( ar_p->rvector[ i_coor2 ] ) );

	if( MATRIX_SCALAR_PRODUCT( *dummy_matrix2_p, dummy, *dummy_matrix2_p  ) ) info=1;  //WARNING: overwriting!
    
	if( MATRIX_SUM( *dummy_matrix1_p, *dummy_matrix2_p, *dummy_matrix1_p  ) ) info=1;  //WARNING: overwriting!

      } /* j_mode loop */

    } /* i_mode loop */

    dummy = CMPLX_DIVISION( CMPLX( 0.5e0 ), IHBAR ); // WARNING: this coefficient depends on the way the double loop i_mode j_mode is handled

    if( MATRIX_SCALAR_PRODUCT( *dummy_matrix1_p, dummy, *dummy_matrix1_p ) );  //WARNING: overwriting!

    /*
      fprintf( stdout, "dummy_matrix1\n" );
      if( MATRIX_PRINT_PLUS( stdout, *dummy_matrix1_p ) ) info=1;
    */

    if( MATRIX_SUM( rho_dot_p->array[ i_rho ], *dummy_matrix1_p, rho_dot_p->array[ i_rho ] ) ) info=1;  //WARNING: overwriting!


#ifdef __DEBUG_PLUS__
    fprintf( stdout, "   --------------\n" );
    norm = MATRIX_NORM( *dummy_matrix1_p );
    fprintf( stdout, "   Term's norm [5] = %le\n", norm );
#endif /* __DEBUG_PLUS__ */


#endif /* __ENERGY_MOD__ */


    /* if rho_dot is too small, set it to zero */
#ifdef __RHO_DOT_CHECK__

    if( MATRIX_NORM( rho_dot_p->array[ i_rho ] ) < EPS_LOC ){

      if( MATRIX_ZERO( rho_dot_p->array[ i_rho ] ) ) info=1;

    } 

#endif /* __RHO_DOT_CHECK__ */


#ifdef __DEBUG__

    if( CHECK_PARTICLE_HOLE_MATRIX( state_p->rho_dot.array[ i_rho ] ) ) info=1;

#endif /* __DEBUG__ */


    if( info ) break;

  } /* end main loop */


#ifdef __DEBUG_PLUS__
  fprintf( stdout, "\n" );
  fprintf( stdout, "-----------\n" );
  fprintf( stdout, "rho_dot\n" );
  if( MATRIX_ARRAY_PRINT_PLUS( stdout, *rho_dot_p ) ) info=1;
  fprintf( stdout, "-----------\n" );
  fprintf( stdout, "\n" );
#endif /* __DEBUG_PLUS__ */


  return info;

}

//------------------------------------------

/* utilities */

//------------------------------------------

#ifdef __DEBUG__

/* rho_dot_is_zero */

int rho_dot_is_zero( const constants constants, state_p state_p, config_p config_p ){

  /* constants */
  int             max_rho_index;
  /* state */
  matrix_array_p  rho_dot_p;
  /* dummies */
  int             index;
  int             info=0;
  double          norm;


  max_rho_index =  constants.max_rho_index;

  rho_dot_p     = &(state_p->rho_dot);


  for( index=0; index<max_rho_index; index++ ){

    norm = MATRIX_NORM( rho_dot_p->array[ index ] );

    if( norm > EPS_LOC ){

      fprintf( stderr, "ERROR: the entry %d of rho_dot is non-zero [ %le ]\n", index, norm );
      fflush( stderr );

      info=1;

      break;

    }


  }




  return info;

}

//------------------------------------------

#endif /* __DEBUG__ */
