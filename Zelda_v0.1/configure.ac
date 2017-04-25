#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.68])
AC_INIT([Zelda], [0.2], [philipp.kaempfer@h-its.org])
AM_INIT_AUTOMAKE([Zelda],[0.2])
AC_PREFIX_DEFAULT ([$(pwd]))
AC_CONFIG_SRCDIR([include/readDB.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([math.h stdint.h stdlib.h string.h unistd.h time.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_SIZE_T
AC_HEADER_STDBOOL
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime memset mkdir pow strchr strrchr])


AC_OUTPUT(Makefile src/Makefile)
