Reference implementation files
==============================

This directory contains javascript implementations of the custom
codecs using in CRAM 3.1, capable of being run under node.js.

These is not written for speed, but for clarity and as an exercise in
checking the pseudocode in the CRAM specification.  It is written as
close to this pseudocode as is possible.


Prerequisites: minimist package for command line parsing and bzip2 for
part of the arith_gen.js code.

    npm install minimist
    npm install bzip2


iostream.js
-----------

Makes a buffer appear to be a stream with ReadByte, ReadITF8, etc
functions.


rans.js
-------

Implements the order-0 and order-1 rans (4x8) decoder as used in CRAM3.0.


main_rans.js
------------

A command line tool to exercise the rans.js code, included for debug
purposes.


rans4x16.js, main_rans4x16.js
-----------------------------

A 16-bit renormalising variant of rANS above.  This also includes
transforms for RLE, bit-packing and 4-way interleaving.


arith_sh.js
-----------

Arithmetic (range) coding with Schindler carry handling.

byte_model.js
-------------

An adaptive model for keeping track of symbol frequencies.

arith_gen.js, main_arith_gen.js
-------------------------------

Wrapper around arith_sh.js to perform order-0/1 encoding with RLE and
bit-packing.  Plus debug command line tool


fqzcomp.js, main_fqzcomp.js
---------------------------

Implements the fqzcomp quality compression codec. Plus debug command
line tool.


tok3.js, main_tok3.js
---------------------

Implements the tokenise_name3 read identifier compression codec.
Plus debug command line tool.


Testing
=======

The various main js files can be used for adhoc testing.  There is
also a Makefile which performs checks against known defined data
streams and does round-trip testing in both Javascript and if compiled
the C variant.  You can set CORPUS make variable to a larger data set
such htscodecs-corpus.

eg.

    make check CORPUS=../tests/htscodecs-corpus/
