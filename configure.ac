#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([slater_stoch],[1.0.0],[https://github.com/scemama/slater_stoch/issues],[],[])
AC_CONFIG_SRCDIR([configure.ac])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_AUX_DIR([build])

AM_INIT_AUTOMAKE([subdir-objects color-tests parallel-tests silent-rules 1.11])
AM_SILENT_RULES(yes)

AC_LANG(Fortran)

# Checks for programs.
#AC_PROG_F77
AX_PROG_FC_MPI
F77=$FC

AS_CASE(["$FC"],
[ifort*], [FCFLAGS="$FCFLAGS -mkl=sequential"], 
[mpiifort*],  [FCFLAGS="$FCFLAGS -mkl=sequential"],
[

# Checks for libraries.
## BLAS
AC_LANG_PUSH(C)
AX_BLAS([], [AC_MSG_ERROR([BLAS was not found.])])
AC_LANG_POP()
LIBS="$LIBS $BLAS_LIBS"
])




# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_LANG_PUSH(C)
AC_CHECK_FUNC([vderfcinv],[ AC_DEFINE([HAVE_VDERFINVC],[1],[vdErfInvc found]) ])

AC_ARG_ENABLE([mmap],
            [AS_HELP_STRING([--enable-mmap],
                            [Enable support for mmap])],
            [enable_mmap=$enableval],
            [enable_mmap=no])

AM_CONDITIONAL([HAVE_MMAP], [test "$enable_mmap" = "yes"])

AS_IF([test "$enable_mmap" = "yes"], [AC_DEFINE([HAVE_MMAP], [1], [Use mmap in SVD cleaning])])

AC_ARG_WITH([trexio],
            [AS_HELP_STRING([--without-trexio],
                            [disable support for TREXIO])],
            [with_trexio=$withval],
            [with_trexio=yes])

AS_CASE([$with_trexio],
  [no], [],
  [yes],
  [AC_CHECK_LIB([trexio], [trexio_open], [],
                [PKG_CHECK_MODULES([TREXIO], [trexio]) ])
   CFLAGS="$CFLAGS $TREXIO_CFLAGS"
   LIBS="$LIBS $TREXIO_LIBS"
   AC_DEFINE([HAVE_TREXIO], [1], [Define if TREXIO support is available])
   ARGS="${ARGS} trexio"],
  [AS_IF([test ! -d "$withval"],
         [AC_MSG_ERROR([--with-trexio=$withval: argument is not a directory])])
   TREXIO_LIBS="-L$with_trexio/lib -ltrexio"
   TREXIO_CFLAGS="-I$with_trexio/include"
   CFLAGS="$CFLAGS $TREXIO_CFLAGS"
   LIBS="$LIBS $TREXIO_LIBS"
   AC_DEFINE([HAVE_TREXIO], [1], [Define if TREXIO support is available])
   ARGS="${ARGS} trexio"
  ])
AC_LANG_POP


AC_CONFIG_FILES([Makefile])
AC_OUTPUT
