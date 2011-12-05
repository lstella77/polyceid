
dnl  Copyright (C) 2011 by Lorenzo Stella <lorenzo DOT stella77 AT gmail DOT com>

dnl  Permission is hereby granted, free of charge, to any person obtaining a copy
dnl  of this software and associated documentation files (the "Software"), to deal
dnl  in the Software without restriction, including without limitation the rights
dnl  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
dnl  copies of the Software, and to permit persons to whom the Software is
dnl  furnished to do so, subject to the following conditions:

dnl  The above copyright notice and this permission notice shall be included in
dnl  all copies or substantial portions of the Software.

dnl  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
dnl  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
dnl  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
dnl  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
dnl  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
dnl  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
dnl  THE SOFTWARE.

AC_INIT([PolyCEID],[0.1],[lorenzo DOT stella77 AT gmail DOT com])
#AC_CONFIG_SRCDIR([./src/PolyCEID/PolyCEID_main.c])
AC_CONFIG_MACRO_DIR([m4])
AC_CANONICAL_TARGET
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

AC_ARG_ENABLE([check-symmetries],AC_HELP_STRING([--enable-check-symmetries],[check the symmetries]))
if test "${enable_check_symmetries}" == ""; then
	enable_check_symmetries=no
fi
if test "${enable_check_symmetries}" == "yes"; then
	AC_DEFINE([__CHECK_SYMMETRIES__],1,[set to 1 if you want to check symmetries])	
fi

AC_ARG_ENABLE([k-matrix-update],AC_HELP_STRING([--enable-k-matrix-update],[update the K_matrix]))
if test "${enable_k_matrix_update}" == ""; then
	enable_k_matrix_update=no
fi
if test "${enable_k_matrix_update}" == "yes"; then
	AC_DEFINE([__K_MATRIX_UPDATE__],1,[set to 1 if you want to update the K_matrix])	
fi

AC_ARG_ENABLE([time-scale],AC_HELP_STRING([--enable-time-scale],[use scaled time]))
if test "${enable_time_scale}" == ""; then
	enable_time_scale=no
fi
if test "${enable_time_scale}" == "yes"; then
	AC_DEFINE([__TIME_SCALE__],1,[set to 1 if you want to use scaled time])	
fi

AC_ARG_ENABLE([forces_non_cons],AC_HELP_STRING([--enable-forces_non_cons],[use non-conservative forces]))
if test "${enable_forces_non_cons}" == ""; then
	enable_forces_non_cons=no
fi
if test "${enable_forces_non_cons}" == "yes"; then
	AC_DEFINE([__FORCES_NON_CONS__],1,[set to 1 if you want to use non-conservative forces])	
fi

AC_PROG_CC
#AC_PROG_F77
AC_PROG_FC
F77="${FC}"
FFLAGS="${FCFLAGS}"
AC_PROG_RANLIB
ACX_BLAS 
ACX_LAPACK 

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile \
                  examples/Makefile \
                  examples/example01/Makefile \
                  examples/example01/reference/Makefile \
                  src/Makefile \
		  src/PolyCEID/Makefile \
		  src/PolyCEID/algebra/Makefile \
		  src/PolyCEID/evolution/Makefile \
		  src/PolyCEID/hamiltonians/Makefile \
		  src/PolyCEID/initial/Makefile \
		  src/PolyCEID/input/Makefile \
		  src/PolyCEID/output/Makefile \
		  src/PolyCEID/structures/Makefile \
		  src/PolyCEID/utils/Makefile])
AC_OUTPUT

