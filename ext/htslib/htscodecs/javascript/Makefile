all: test

CBIN=../build/tests
NODE=node
NODE_OPTS=--use-strict

TESTS=test_r4x8 test_r4x16 test_arith test_fqzcomp test_tok3
CORPUS=../tests

modules: node_modules/bzip2 node_modules/minimist

node_modules/bzip2:
	-mkdir node_modules
	npm install bzip2

node_modules/minimist:
	-mkdir node_modules
	npm install minimist

test check: modules ${TESTS}

test_r4x8:
	@echo
	@echo === Checking r4x8
	@if [ ! -e ${CBIN}/rans4x8 ]; then echo "Set CBIN if you wish to validate against C version"; fi
	@for i in ${CORPUS}/dat/r4x8/*; do \
            echo $$i;\
	    base=`echo $$i | sed 's/\.[0-9]*$$//;s#/q#/../q#'`; \
	    level=`echo $$i | sed 's/.*\.//'`;\
	    # Decode predefined data set\
	    a=`${NODE} ${NODE_OPTS} main_rans.js -d -r $$i 2>/dev/null | md5sum`; \
	    b=`cut -f 1 < $$base | tr -d '\012' | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Mismatch for $$i; \
	    # Round trip. \
	    cut -f 1 < $$base | tr -d '\012' > _nonl; \
	    ${NODE} ${NODE_OPTS} main_rans.js -o $$level -r _nonl 2>/dev/null > _; \
	    a=`${NODE} ${NODE_OPTS} main_rans.js -d -r _ 2>/dev/null | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Fail round-trip for $$base with level $$level; \
	    if [ -e ${CBIN}/rans4x8 ]; \
	    then \
	      a=`${CBIN}/rans4x8 -d -r < _ 2>/dev/null | tr '\000' '\012' | md5sum`; \
	      test "$$a" = "$$b" || echo $<: Fail JS to C round-trip for $$base; \
	    fi; \
	done

test_r4x16:
	@echo
	@echo === Checking r4x16
	@if [ ! -e ${CBIN}/rans4x16pr ]; then echo "Set CBIN if you wish to validate against C version"; fi
	@for i in ${CORPUS}/dat/r4x16/*; do \
            echo $$i;\
	    base=`echo $$i | sed 's/\.[0-9]*$$//;s#/q#/../q#'`; \
	    level=`echo $$i | sed 's/.*\.//'`;\
	    # Decode predefined data set\
	    a=`${NODE} ${NODE_OPTS} main_rans4x16.js -d -r $$i 2>/dev/null | md5sum`; \
	    b=`cut -f 1 < $$base | tr -d '\012' | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Mismatch for $$i; \
	    # Round trip. \
	    cut -f 1 < $$base | tr -d '\012' > _nonl; \
	    ${NODE} ${NODE_OPTS} main_rans4x16.js -o $$level -r _nonl 2>/dev/null > _; \
	    a=`${NODE} ${NODE_OPTS} main_rans4x16.js -d -r _ 2>/dev/null | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Fail round-trip for $$base with level $$level; \
	    if [ -e ${CBIN}/rans4x16pr ]; \
	    then \
	      a=`${CBIN}/rans4x16pr -d -r < _ 2>/dev/null | tr '\000' '\012' | md5sum`; \
	      test "$$a" = "$$b" || echo $<: Fail JS to C round-trip for $$base; \
	    fi; \
	done

test_arith:
	@echo
	@echo === Checking arith
	@if [ ! -e ${CBIN}/arith_dynamic ]; then echo "Set CBIN if you wish to validate against C version"; fi
	@for i in ${CORPUS}/dat/arith/*; do \
            echo $$i;\
	    base=`echo $$i | sed 's/\.[0-9]*$$//;s#/q#/../q#'`; \
	    level=`echo $$i | sed 's/.*\.//'`;\
	    # Decode predefined data set\
	    a=`${NODE} ${NODE_OPTS} main_arith_gen.js -d -r $$i 2>/dev/null | md5sum`; \
	    b=`cut -f 1 < $$base | tr -d '\012' | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Mismatch for $$i; \
	    # Round trip. \
	    cut -f 1 < $$base | tr -d '\012' > _nonl; \
	    ${NODE} ${NODE_OPTS} main_arith_gen.js -o $$level -r _nonl 2>/dev/null > _; \
	    a=`${NODE} ${NODE_OPTS} main_arith_gen.js -d -r _ 2>/dev/null | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Fail round-trip for $$base with level $$level; \
	    if [ -e ${CBIN}/arith_dynamic ]; \
	    then \
	      a=`${CBIN}/arith_dynamic -d -r < _ 2>/dev/null | tr '\000' '\012' | md5sum`; \
	      test "$$a" = "$$b" || echo $<: Fail JS to C round-trip for $$base; \
	    fi; \
	done

test_fqzcomp:
	@echo
	@echo === Checking fqzcomp
	@if [ ! -e ${CBIN}/fqzcomp_qual ]; then echo "Set CBIN if you wish to validate against C version"; fi
	@for i in ${CORPUS}/dat/fqzcomp/q*; do \
            echo $$i;\
	    base=`echo $$i | sed 's/\.[0-9]*$$//;s#/q#/../q#'`; \
	    level=`echo $$i | sed 's/.*\.//'`;\
	    # Decode predefined data set\
	    a=`${NODE} ${NODE_OPTS} main_fqzcomp.js -d -r $$i 2>/dev/null | md5sum`; \
	    b=`awk '{print $$1}' $$base | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Mismatch for $$i; \
	    # Round trip. \
	    ${NODE} ${NODE_OPTS} main_fqzcomp.js -s $$level -r $$base 2>/dev/null > _; \
	    a=`${NODE} ${NODE_OPTS} main_fqzcomp -d -r _ 2>/dev/null | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Fail round-trip for $$base with level $$level; \
	    if [ -e ${CBIN}/fqzcomp_qual ]; \
	    then \
	      a=`${CBIN}/fqzcomp_qual -d -r < _ 2>/dev/null | tr '\000' '\012' | md5sum`; \
	      test "$$a" = "$$b" || echo $<: Fail JS to C round-trip for $$base; \
	    fi; \
	done

test_tok3:
	@echo
	@echo === Checking tok3
	@if [ ! -e ${CBIN}/tokenise_name3 ]; then echo "Set CBIN if you wish to validate against C version"; fi
	@for base in ${CORPUS}/names/*.names; do \
	    echo -n "$$base  	";\
	    ${NODE} ${NODE_OPTS} main_tok3.js -a -r $$base 2>/dev/null > _; \
	    a=`${NODE} ${NODE_OPTS} main_tok3 -d -r _ 2>/dev/null | md5sum`; \
	    cat _ | wc -c;\
	    b=`cat $$base | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Fail round-trip for $$base; \
	    if [ -e ${CBIN}/tokenise_name3 ]; \
	    then \
	      a=`${CBIN}/tokenise_name3 -d -r < _ | tr '\000' '\012' | md5sum`; \
	      test "$$a" = "$$b" || echo $<: Fail JS to C round-trip for $$base; \
	    fi; \
	done;
	@for i in ${CORPUS}/names/tok3/*; do \
	   echo $$i;\
	    base=`echo $$i | sed 's/\.[0-9]*$$//;s#/tok3##'`; \
	    level=`echo $$i | sed 's/.*\.//'`;\
	    # Decode predefined data set\
	    a=`${NODE} ${NODE_OPTS} main_tok3 -d -r $$i 2>/dev/null | md5sum`; \
	    b=`cat $$base | md5sum`; \
	    test "$$a" = "$$b" || echo $<: Mismatch for $$i; \
	done;
