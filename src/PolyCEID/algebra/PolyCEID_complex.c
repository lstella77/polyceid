
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
#include "PolyCEID_complex.h"


complex __my_complex_cast( const double real ){

  complex z;

  z.z[0] = real;
  z.z[1] = 0.0e0;

  return z;

}

complex __my_complex_initialize( const double real, const double imaginary ){

  complex z;

  z.z[0] = real;
  z.z[1] = imaginary;

  return z;

}

double __my_complex_real_part( const complex z ){

  return z.z[0];

}

double __my_complex_imaginary_part( const complex z ){

  return z.z[1];

}

/*
  double __my_complex_norm( const complex z ){

  int    n=1;
  int    incx=1;
  double norm;


  //  norm = z.z[0] *z.z[0] +z.z[1] *z.z[1];

  //  norm = sqrt( norm );

  return norm;

  }
*/

complex __my_complex_times_I( const complex z ){

  complex zI;


  zI.z[0] = -z.z[1];
  zI.z[1] =  z.z[0];

  return zI;

}

complex __my_complex_conjugate( const complex z ){


  complex zC;

  zC.z[0] =  z.z[0];
  zC.z[1] = -z.z[1];

  return zC;

}

complex __my_complex_inverse( const complex z ){

  complex z_inverse;
  double  norm2;


  norm2 = z.z[0] *z.z[0] +z.z[1] *z.z[1];

  if( fabs( norm2  ) < DBL_EPSILON ) return CMPLX_HUGE_VAL;

  z_inverse.z[0] =  z.z[0] /norm2;
  z_inverse.z[1] = -z.z[1] /norm2;
  
  return z_inverse;

}

complex __my_complex_sum( const complex z1, const complex z2 ){

  complex z3;


  z3.z[0] = z1.z[0] + z2.z[0];
  z3.z[1] = z1.z[1] + z2.z[1];

  return z3;

}

complex __my_complex_dif( const complex z1, const complex z2 ){

  complex z3;


  z3.z[0] = z1.z[0] - z2.z[0];
  z3.z[1] = z1.z[1] - z2.z[1];

  return z3;

}

complex __my_complex_product( const complex z1, const complex z2 ){

  complex z3;


  z3.z[0] = z1.z[0] *z2.z[0] -z1.z[1] *z2.z[1] ;
  z3.z[1] = z1.z[0] *z2.z[1] +z1.z[1] *z2.z[0];


  return z3;

}

complex __my_complex_division( const complex z1 , const complex z2 ){

  /* dummies */
  complex z3, z2_inverse;
  int     info=0;


  /*
    fprintf( stdout, "z2\n" ); 
    if( CMPLX_PRINT_PLUS( stdout, z2 ) ) info=1;
    fprintf( stdout, "\n" ); 
  */


  z2_inverse = CMPLX_INVERSE( z2 );


  /*
    fprintf( stdout, "z2_inverse\n" ); 
    if( CMPLX_PRINT_PLUS( stdout, z2_inverse ) ) info=1;
    fprintf( stdout, "\n" ); 
  */

  z3  = CMPLX_PRODUCT( z1, z2_inverse );


  if( !info ){

    return z3;

  }
  else{

    return CMPLX_HUGE_VAL;

  }

}

complex __my_complex_phase_fix( const complex a ){

  /* dummies */
  double   phase;
  complex  z;

  
  if( fabs( a.z[1] ) > EPS ){
  
    // phase of complex
    phase = atan( ( a.z[1] ) /( a.z[0] ) ); 

    // rotating complex to have its Im(complex)=0
    z.z[0] = ( a.z[0] ) *cos( phase ) +( a.z[1] ) *sin( phase );
    z.z[1] = 0.0e0;

  }
  else{

    z.z[0]=a.z[0];
    z.z[1]=0.0e0;

  }        


  return z;

}	

int     __my_complex_copy( complex_p z1_p, const complex z2 ){

  z1_p->z[0] = z2.z[0];
  z1_p->z[1] = z2.z[1];


  return 0;

}

int     __my_complex_compare( const complex z1, const complex z2 ){

  /* dummies */
  int info=0;


  if( fabs( z1.z[0] -z2.z[0] ) > FLT_EPSILON ||
      fabs( z1.z[1] -z2.z[1] ) > FLT_EPSILON ){

    fprintf( stderr, "ERROR: the two complex numbers are not equal\n" );
    fflush( stderr );

    info=1;

  }


  return info;

}

int __my_complex_read( FILE* file_p, complex_p zp ){

  /* dummies */
  int info=0;


  if( FILE_CHECK( file_p, __my_complex_print  ) ) return 1;

  if( fscanf( file_p, "%le", &zp->z[0] ) < 1 ) info=1;
  if( fscanf( file_p, "%le", &zp->z[1] ) < 1  ) info=1;
  //  fflush( file_p ); // Why?


  return info;

}

int __my_complex_print( FILE* file_p, const complex z ){

  if( FILE_CHECK( file_p, __my_complex_print  ) ) return 1;

  fprintf( file_p, DOUBLE_FORMAT"  "DOUBLE_FORMAT, z.z[0], z.z[1] );
  fflush( file_p );


  return 0;

}

int __my_complex_print_plus( FILE* file_p, const complex z ){

  if( FILE_CHECK( file_p, __my_complex_print  ) ) return 1;

  fprintf( file_p, "# ( "DOUBLE_FORMAT", "DOUBLE_FORMAT" )", z.z[0], z.z[1] );
  fflush( file_p );


  return 0;

}

int __my_complex_print_plus_var( FILE* file_p, const complex z ){

  if( FILE_CHECK( file_p, __my_complex_print  ) ) return 1;

  fprintf( file_p, "( "DOUBLE_FORMAT", "DOUBLE_FORMAT" )", z.z[0], z.z[1] );
  fflush( file_p );


  return 0;

}
