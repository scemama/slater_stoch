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
AC_PROG_F77
AX_PROG_FC_MPI


# Checks for libraries.
## BLAS
AC_LANG_PUSH(C)
AX_BLAS([], [AC_MSG_ERROR([BLAS was not found.])])
AC_LANG_POP()




# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
