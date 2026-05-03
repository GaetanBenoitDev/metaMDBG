# The idea is the same used as in AC_SEARCH_LIB, but unlike that the
# no-library case is the last scenario instead of the first one.
# The reason for this is to handle scenarios where the external library
# is preferred over the internal C library implementation.  An example of
# this is FreeBSD's pthread functionality, which has stub functions in the C
# library that do not work.  See htscodecs issue#64

# AX_SEARCH_LIBS_REV(FUNCTION, SEARCH-LIB-PATH
AC_DEFUN([AX_SEARCH_LIBS_REV],
[
dnl Create an input to test linking a file calling $1.
dnl Used by AC_LINK_IFELSE.
AC_LANG_CONFTEST([AC_LANG_CALL([], [$1])])
_found=no
_LIBS=$LIBS
AC_MSG_CHECKING([for function $1])
for xlib in $2 ""
do
    if test "x$xlib" != "x"
    then
       LIBS="-l$xlib $_LIBS"
       _res="-l$xlib"
    else
       LIBS="$_LIBS"
       _res="(no library needed)"
    fi

    AC_LINK_IFELSE([], [_found=yes; break])
done

if test "$_found" = "yes"
then
    AC_MSG_RESULT([$_res])
else
    AC_MSG_RESULT([not found])
    AC_MSG_ERROR([Function $1 not found])
fi

unset _found
unset _res
unset _LIBS
])
