dnl @synopsis VL_PROG_CC_WARNINGS([ANSI])
dnl
dnl From http://ac-archive.sourceforge.net/ac-archive/vl_prog_cc_warnings.html
dnl
dnl Enables a reasonable set of warnings for the C compiler.
dnl Optionally, if the first argument is nonempty, turns on flags which
dnl enforce and/or enable proper ANSI C if such are known with the
dnl compiler used.
dnl
dnl Currently this macro knows about GCC, Solaris C compiler, Digital
dnl Unix C compiler, C for AIX Compiler, HP-UX C compiler, IRIX C
dnl compiler, NEC SX-5 (Super-UX 10) C compiler, and Cray J90 (Unicos
dnl 10.0.0.8) C compiler.
dnl
dnl @category C
dnl @author Ville Laurikari <vl@iki.fi>
dnl Updated by Rob Davies
dnl @version 2002-04-04
dnl @license AllPermissive
dnl Copying and distribution of this file, with or without modification,
dnl are permitted in any medium without royalty provided the copyright notice
dnl and this notice are preserved. Users of this software should generally
dnl follow the principles of the MIT License including its disclaimer.
dnl Original Copyright (c) Ville Laurikari 2002
dnl Modifications Copyright (c) Genome Research Limited 2015

AC_DEFUN([VL_PROG_CC_WARNINGS], [
  AC_ARG_ENABLE([warnings],
    [AS_HELP_STRING([--disable-warnings], [turn off compiler warnings])],
    [],
    [enable_warnings=yes])

  AS_IF([test "x$enable_warnings" != xno],[
    AC_PROG_GREP

    ansi="$1"
    AS_IF([test "x$ansi" = "x"],
          [msg="for C compiler warning flags"],
          [msg="for C compiler warning and ANSI conformance flags"])

    AC_CACHE_CHECK($msg, vl_cv_prog_cc_warnings, [
      vl_cv_prog_cc_warnings=""
      AS_IF([test "x$CC" != "x"],[
        cat > conftest.c <<EOF
int main(int argc, char **argv) { return 0; }
EOF

        dnl Most compilers print some kind of a version string with some command
        dnl line options (often "-V").  The version string should be checked
        dnl before doing a test compilation run with compiler-specific flags.
        dnl This is because some compilers (like the Cray compiler) only
        dnl produce a warning message for unknown flags instead of returning
        dnl an error, resulting in a false positive.  Also, compilers may do
        dnl erratic things when invoked with flags meant for a different
        dnl compiler.

        dnl GCC
        AS_IF([test "x$GCC" = "xyes"],[
	  AS_IF([test "x$ansi" = "x"],
	        [vl_cv_prog_cc_warnings="-Wall"],
		[vl_cv_prog_cc_warnings="-Wall -ansi -pedantic"])
	],
        dnl Solaris C compiler
        ["$CC" -V 2>&1 | $GREP -i "WorkShop" > /dev/null 2>&1 &&
         "$CC" -c -v -Xc conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[
	  AS_IF([test "x$ansi" = "x"],
	        [vl_cv_prog_cc_warnings="-v"],
		[vl_cv_prog_cc_warnings="-v -Xc"])
	],
        dnl Digital Unix C compiler
        ["$CC" -V 2>&1 | $GREP -i "Digital UNIX Compiler" > /dev/null 2>&1 &&
         "$CC" -c -verbose -w0 -warnprotos -std1 conftest.c > /dev/null 2>&1 &&
         test -f conftest.o], [
	  AS_IF([test "x$ansi" = "x"],
	   	[vl_cv_prog_cc_warnings="-verbose -w0 -warnprotos"],
		[vl_cv_prog_cc_warnings="-verbose -w0 -warnprotos -std1"])
	],
        dnl C for AIX Compiler
        ["$CC" 2>&1 | $GREP -i "C for AIX Compiler" > /dev/null 2>&1 &&
         "$CC" -c -qlanglvl=ansi -qinfo=all conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[
	  AS_IF([test "x$ansi" = "x"],
	        [vl_cv_prog_cc_warnings="-qsrcmsg -qinfo=all:noppt:noppc:noobs:nocnd"],
		[vl_cv_prog_cc_warnings="-qsrcmsg -qinfo=all:noppt:noppc:noobs:nocnd -qlanglvl=ansi"])
	],
        dnl IRIX C compiler
        ["$CC" -version 2>&1 | $GREP -i "MIPSpro Compilers" > /dev/null 2>&1 &&
         "$CC" -c -fullwarn -ansi -ansiE conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[
	  AS_IF([test "x$ansi" = "x"],
	        [vl_cv_prog_cc_warnings="-fullwarn"],
		[vl_cv_prog_cc_warnings="-fullwarn -ansi -ansiE"])
	],
        dnl HP-UX C compiler
        [what "$CC" 2>&1 | $GREP -i "HP C Compiler" > /dev/null 2>&1 &&
         "$CC" -c -Aa +w1 conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[
	  AS_IF([test "x$ansi" = "x"],
	        [vl_cv_prog_cc_warnings="+w1"],
		[vl_cv_prog_cc_warnings="+w1 -Aa"])
	],
        dnl The NEC SX-5 (Super-UX 10) C compiler
        ["$CC" -V 2>&1 | $GREP "/SX" > /dev/null 2>&1 &&
         "$CC" -c -pvctl[,]fullmsg -Xc conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[
	  AS_IF([test "x$ansi" = "x"],
	        [vl_cv_prog_cc_warnings="-pvctl[,]fullmsg"],
		[vl_cv_prog_cc_warnings="-pvctl[,]fullmsg -Xc"])
	],
        dnl The Cray C compiler (Unicos)
        ["$CC" -V 2>&1 | $GREP -i "Cray" > /dev/null 2>&1 &&
         "$CC" -c -h msglevel 2 conftest.c > /dev/null 2>&1 &&
         test -f conftest.o],[
	  AS_IF([test "x$ansi" = "x"],
	        [vl_cv_prog_cc_warnings="-h msglevel 2"],
		[vl_cv_prog_cc_warnings="-h msglevel 2 -h conform"])
        ])
        rm -f conftest.*
      ])
    ])
    AS_IF([test "x$vl_cv_prog_cc_warnings" != "x"],
          [CFLAGS="$vl_cv_prog_cc_warnings $CFLAGS"])
  ])
])dnl
