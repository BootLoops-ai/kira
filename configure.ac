#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([kira],[3.1],[jusovitsch@googlemail.com])
AM_INIT_AUTOMAKE([subdir-objects])
AC_CONFIG_SRCDIR([src/kira/dataBase.cpp])
AC_CONFIG_HEADERS([src/pyred/config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC


AC_ARG_ENABLE(firefly,[AS_HELP_STRING([--enable-firefly], [enable support for algorithms to interpolate rational functions via the finite field meethods])])
AC_MSG_CHECKING([whether to enable firefly])
AS_IF([test "x${enable_firefly}" = "xyes" ], AC_MSG_RESULT([yes]), AC_MSG_RESULT([no]))
AM_CONDITIONAL([FIREFLY],[test "x${enable_firefly}" = "xyes"])

AC_ARG_ENABLE(debug,[AS_HELP_STRING([--enable-debug], [use valgrind modus])])
AC_MSG_CHECKING([whether to enable debug])
AS_IF([test "x${enable_debug}" = "xyes" ], AC_MSG_RESULT([yes]), AC_MSG_RESULT([no]))
AM_CONDITIONAL([DEBUG],[test "x${enable_debug}" = "xyes"])

AC_ARG_ENABLE(weight128,[AS_HELP_STRING([--enable-weight128], [use weight bit 128])])
AC_MSG_CHECKING([whether to enable weight bit 128])
AS_IF([test "x${enable_weight128}" = "xyes" ], AC_MSG_RESULT([yes]), AC_MSG_RESULT([no]))
AM_CONDITIONAL([WEIGHT128],[test "x${enable_weight128}" = "xyes"])

# Checks for libraries.
AC_CHECK_LIB([dl], [main])
PKG_CHECK_MODULES(GINAC, ginac)
PKG_CHECK_MODULES(CLN, cln)
PKG_CHECK_MODULES(YAML_CPP, yaml-cpp)
PKG_CHECK_MODULES(ZLIB, zlib)
AM_COND_IF( [FIREFLY],[
PKG_CHECK_MODULES(FIREFLY, firefly)
    ]
)

AC_CHECK_LIB(gmp, __gmpz_init, , [AC_MSG_ERROR([GNU MP not found, see http://gmplib.org/])])
#AC_CHECK_LIB(flint, n_mulmod2_preinv, , [AC_MSG_ERROR([Flint not found])])
#AC_CHECK_LIB(FIREFLY, firefly_exists, , [AC_MSG_ERROR([FireFly not found])])
# Checks for header files.
AC_CHECK_HEADERS([fcntl.h inttypes.h limits.h malloc.h stddef.h stdint.h stdlib.h string.h sys/file.h sys/ioctl.h sys/mount.h sys/param.h sys/statvfs.h sys/time.h unistd.h utime.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UID_T
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_MODE_T
AC_TYPE_OFF_T
AC_TYPE_PID_T
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_CHECK_MEMBERS([struct stat.st_blksize])
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRERROR_R
AC_FUNC_STRTOD
AC_CHECK_FUNCS([clock_gettime dup2 fdatasync getpagesize gettimeofday localtime_r memmove memset mkdir pow select strchr strerror strstr strtol utime])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

AC_MSG_WARN([The Autotools build system is deprecated and will be removed in a future Kira release. Please use the Meson build system.])
