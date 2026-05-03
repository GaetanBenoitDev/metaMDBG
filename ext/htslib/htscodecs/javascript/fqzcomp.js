/*
 * Copyright (c) 2019 Genome Research Ltd.
 * Author(s): James Bonfield
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice,
 *       this list of conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 *
 *    3. Neither the names Genome Research Ltd and Wellcome Trust Sanger
 *       Institute nor the names of its contributors may be used to endorse
 *       or promote products derived from this software without specific
 *       prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL GENOME RESEARCH
 * LTD OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

const IOStream = require("./iostream");
const ByteModel = require("./byte_model");
const RangeCoder = require("./arith_sh");


//----------------------------------------------------------------------
// Main arithmetic entry function: decodes a compressed src and
// returns the uncompressed buffer.

function read_array(src, tab, size) {
    var j = 0; // array value
    var z = 0; // array index: tab[j]
    var last = -1;

    // Remove first level of run-length encoding
    var R = new Array(1024) // runs
    while (z < size) {
	var run = src.ReadByte()
	R[j++] = run
	z += run

	if (run == last) {
	    var copy = src.ReadByte()
	    z += run * copy
	    while (copy--)
		R[j++] = run
	}
	last = run
    }

    // Now expand runs in R to tab, noting 255 is max run
    var i = 0
    j = 0
    z = 0
    while (z < size) {
	var run_len = 0
	do {
	    var part = R[j++]
	    run_len += part
	} while (part == 255)
	
	while (run_len--)
	    tab[z++] = i;
	i++
    }
}

const QMAX = 256

const FLAG_DEDUP  = 2
const FLAG_FLEN   = 4
const FLAG_SEL    = 8    // whether selector is used in context
const FLAG_QMAP   = 16
const FLAG_PTAB   = 32
const FLAG_DTAB   = 64
const FLAG_QTAB   = 128

const GFLAG_MULTI_PARAM = 1
const GFLAG_HAVE_STAB   = 2
const GFLAG_DO_REV      = 4

// Compute a new context from our current state and qual q
function fqz_update_ctx(params, state, q) {
    var last = params.context
    state.qctx = ((state.qctx << params.qshift) + params.qtab[q]); // >>> 0
    last += ((state.qctx & ((1<<params.qbits)-1)) << params.qloc); // >>> 0

    if (params.do_pos)
	last += params.ptab[Math.min(state.p, 1023)] << params.ploc

    if (params.do_delta) {
	last += params.dtab[Math.min(state.delta, 255)] << params.dloc
	// Is it better to use q here or qtab[q]?
	// If qtab[q] we can map eg [a-z0-9A-Z]->0 ,->1 and have
	// delta being a token number count into comma separated lists?
	state.delta += (state.prevq != q) ? 1 : 0
	state.prevq = q
    }

    if (params.do_sel)
	last += state.s << params.sloc

    state.p--

    return last & 0xffff
}

function decode_fqz_single_param(src) {
    var p = {} // params
    
    // Load FQZ parameters
    p.context = src.ReadUint16()
    p.pflags  = src.ReadByte()

    p.do_dedup  = p.pflags & FLAG_DEDUP
    p.fixed_len = p.pflags & FLAG_FLEN
    p.do_sel    = p.pflags & FLAG_SEL
    p.do_qmap   = p.pflags & FLAG_QMAP
    p.do_pos    = p.pflags & FLAG_PTAB
    p.do_delta  = p.pflags & FLAG_DTAB
    p.do_qtab   = p.pflags & FLAG_QTAB

    p.max_sym = src.ReadByte()

    var x = src.ReadByte()
    p.qbits  = x>>4
    p.qshift = x&15
    x = src.ReadByte()
    p.qloc = x>>4
    p.sloc = x&15
    x = src.ReadByte()
    p.ploc = x>>4
    p.dloc = x&15

    // Qual map, eg to "unbin" Illumina qualities
    p.qmap = new Array(256);
    if (p.pflags & FLAG_QMAP) {
	for (var i = 0; i < p.max_sym; i++)
	    p.qmap[i] = src.ReadByte()
    } else {
	// Useful optimisation to speed up main loop
	for (var i = 0; i < 256; i++)
	    p.qmap[i] = i;  // NOP
    }

    // Read tables
    p.qtab = new Array(1024);
    if (p.qbits > 0 && (p.pflags & FLAG_QTAB)) {
	read_array(src, p.qtab, 256)
    } else {
	// Useful optimisation to speed up main loop
	for (var i = 0; i < 256; i++)
	    p.qtab[i] = i;  // NOP
    }

    p.ptab = new Array(1024);
    if (p.pflags & FLAG_PTAB)
	read_array(src, p.ptab, 1024);

    p.dtab = new Array(256);
    if (p.pflags & FLAG_DTAB)
	read_array(src, p.dtab, 256);

    return p
}

function decode_fqz_params(src) {
    var gparams = {
	max_sym: 0
    }

    // Check fqz format version
    var vers = src.ReadByte()
    if (vers != 5) {
	console.error("Invalid FQZComp version number");
	return;
    }

    var gflags = src.ReadByte()
    var nparam = (gflags & GFLAG_MULTI_PARAM) ? src.ReadByte() : 1
    var max_sel = gflags.nparam > 1 ? gflags.nparam-1 : 0 // Note max_sel, not num_sel

    var stab = new Array(256);
    if (gflags & GFLAG_HAVE_STAB) {
	max_sel = src.ReadByte()
	read_array(src, stab, 256);
    } else {
	for (var i = 0; i < nparam; i++)
	    stab[i] = i;
	for (; i < 256; i++)
	    stab[i] = nparam-1;
    }
    gparams.do_rev = (gflags & GFLAG_DO_REV)
    gparams.stab = stab
    gparams.max_sel = max_sel

    gparams.params = new Array(gparams.nparam)
    for (var p = 0; p < nparam; p++) {
	gparams.params[p] = decode_fqz_single_param(src)
	if (gparams.max_sym < gparams.params[p].max_sym)
	    gparams.max_sym = gparams.params[p].max_sym
    }

    return gparams
}

function fqz_create_models(gparams) {
    var model = {}

    model.qual = new Array(1<<16)
    for (var i = 0; i < (1<<16); i++)
	model.qual[i] = new ByteModel(gparams.max_sym+1) // +1 as max value not num. values

    model.len = new Array(4)
    for (var i = 0; i < 4; i++)
	model.len[i] = new ByteModel(256)

    model.rev   = new ByteModel(2)
    model.dup   = new ByteModel(2)

    if (gparams.max_sel > 0)
	model.sel = new ByteModel(gparams.max_sel+1) // +1 as max value not num. values

    return model
}

// Initialise a new record, updating state.
// Returns 1 if dup, otherwise 0
function decode_fqz_new_record(src, rc, gparams, model, state, rev) {
    // Parameter selector
    if (gparams.max_sel > 0) {
	state.s = model.sel.ModelDecode(src, rc)
    } else {
	state.s = 0;
    }
    state.x = gparams.stab[state.s]

    var params = gparams.params[state.x]

    // Reset contexts at the start of each new record
    if (params.fixed_len >= 0) {
	// Not fixed or fixed but first record
	var len = model.len[0].ModelDecode(src, rc)
	len |= model.len[1].ModelDecode(src, rc) << 8
	len |= model.len[2].ModelDecode(src, rc) << 16
	len |= model.len[3].ModelDecode(src, rc) << 24
	if (params.fixed_len > 0)
	    params.fixed_len = -len
    } else {
	len = -params.fixed_len
    }
    state.len = len

    if (gparams.do_rev)
	rev[state.rec] = model.rev.ModelDecode(src, rc)

    state.is_dup = 0
    if (params.pflags & FLAG_DEDUP) {
	if (model.dup.ModelDecode(src, rc))
	    state.is_dup = 1
    }

    state.p = len;  // number of remaining bytes in this record
    state.delta = 0
    state.qctx = 0
    state.prevq = 0
    state.rec++
}

function decode_fqz(src, q_lens) {
    // Decode parameter block
    var n_out = src.ReadUint7()
    var gparams = decode_fqz_params(src)
    if (!gparams) return
    var params = gparams.params
    var rev = new Array(q_lens.length)

    // Create initial models
    var model = fqz_create_models(gparams)

    // Create our entropy encoder and output buffers
    var rc = new RangeCoder(src)
    rc.RangeStartDecode(src)
    var output = new Buffer.allocUnsafe(n_out)

    // Internal FQZ state
    var state = {
	qctx:0,   // Qual-only sub-context
	prevq:0,  // Previous quality value
	delta:0,  // Running delta (q vs prevq)
	p:0,      // Number of bases left in current record
	s:0,      // Current parameter selector value (0 if unused)
	x:0,      // "stab" tabulated copy of s
	len:0,    // Length of current string
	is_dup:0, // This string is a duplicate of last
	rec:0     // Record number
    }

    // The main decode loop itself
    var i = 0     // position in output buffer
    while (i < n_out) {
	if (state.p == 0) {
	    decode_fqz_new_record(src, rc, gparams, model, state, rev)
	    if (state.is_dup > 0) {
		if (model.dup.ModelDecode(src, rc)) {
		    // Duplicate of last line
		    for (var x = 0; x < len; x++)
			output[i+x] = output[i+x-state.len]
		    i += state.len
		    state.p = 0
		    continue
		}
	    }
	    q_lens.push(state.len)

	    var params = gparams.params[state.x]
	    var last = params.context
	}

	// Decode the current quality (possibly mapped via qmap)
	var Q = model.qual[last].ModelDecode(src, rc)

	//if (params.do_qmap)
	//    output[i++] = params.qmap[Q];
	//else
	//    output[i++] = Q
	output[i++] = params.qmap[Q]; // optimised version of above
	last = fqz_update_ctx(params, state, Q)
    }

    if (gparams.do_rev)
	reverse_qualities(output, n_out, rev, q_lens)

    return output;
}

function reverse_qualities(qual, qual_len, rev, len) {
    var rec = 0
    var i = 0
    while (i < qual_len) {
	if (rev[rec]) {
	    var j = 0
	    var k = len[rec]-1
	    while (j < k) {
		var tmp   = qual[i+j]
		qual[i+j] = qual[i+k]
		qual[i+k] = tmp
		j++
		k--
	    }
	}

	i += len[rec++];
    }
}

function decode(src, q_lens) {
    var stream = new IOStream(src);

    //var n_out = stream.ReadUint32(); stream.ReadUint32(); // move to main

    return decode_fqz(stream, q_lens);
}
    
//----------------------------------------------------------------------
// FQZComp encoder.

function pick_fqz_params(src, q_lens, q_dirs, qhist) {
    // Find cardinality of q_dirs
    var qd_last = q_dirs[0];
    for (var i = 0; i < q_dirs.length; i++)
	if (q_dirs[i] != qd_last)
	    break;
    var qd_fixed = (i == q_dirs.length) ? 1 : 0

    // Scan input to find number of symbols and max symbol
    var nsym = 0
    var max_sym = 0

    // selector == 0: Assume one single input dataset
    for (var i = 0; i < 256; i++)
	qhist[0][i] = 0;

    var rec = 0;
    var len = 0
    for (var i = 0; i < src.length; i++) {
	if (len == 0) {
	    len = q_lens[rec < q_lens.length-1 ? rec++ : rec]
	}
	qhist[0][src[i]]++;
	len--;
    }
    for (var i = 0; i < 256; i++) {
	if (!qhist[0][i])
	    continue;
	if (max_sym < i)
	    max_sym = i;
	nsym++;
    }

    var qshift = 5
    var do_qmap = 0
    // Reduced symbol frequencies implies lower qshift and
    // a lookup table to go from qual to Q
    if (nsym <= 16) {
	do_qmap = 1 // based on qhist
	if (nsym <= 2)
	    qshift = 1
	else if (nsym <= 4)
	    qshift = 2
	else if (nsym <= 8)
	    qshift = 3
	else
	    qshift = 4
    }

//    // Two params and a 1-bit selector.
//    // This is 1% overhead vs two data sets compressed independently.
//    // It's 6.9% smaller than compressing both together with 1 param.
//    if (0) return [{
//	// q4
//	qbits:     8,
//	qshift:    2,
//	qloc:      7,
//
//	pbits:     7,
//	pshift:    1,
//	ploc:      0,
//
//	dbits:     0,
//	dshift:    0,
//	dloc:      0,
//
//      sbits:     0,
//      sloc:      0,
//
//	//sbits:     2,
//	//do_stab:   1,
//	sbits:     1,
//	do_stab:   0,
//	context:   (0<<15),
//
//	max_sym:   36,
//	nsym:      4,
//
//	do_qmap:   1,
//	do_dedup:  0,
//	fixed_len: 1,
//	do_sel:  0,
//	do_rev:    0,
//	do_pos:    1,
//	do_delta:  0,
//	do_qtab:   0
//    }, {
//	//q40
//	qbits:     9,
//	qshift:    5,
//	qloc:      7,
//
//	pbits:     7,
//	pshift:    0,
//	ploc:      0,
//
//	dbits:     0,
//	dshift:    0,
//	dloc:      0,
//
//      sbits:     0,
//      sloc:      0,
//
//	//sbits:     2,
//	//do_stab:   1,
//	sbits:     1,
//	do_stab:   0,
//	context:   (1<<15),
//
//	max_sym:   44,
//	nsym:      45,
//
//	do_qmap:   0,
//	do_dedup:  0,
//	fixed_len: 1,
//	do_sel:  0,
//	do_rev:    0,
//	do_pos:    1,
//	do_delta:  0,
//	do_qtab:   0
//    }]

    return [{qbits:     8+(qshift>4),
	     qshift:    qshift,
	     qloc:      7,

	     pbits:     7,
	     pshift:    q_lens[0] > 128 ? 1 : 0,
	     ploc:      0,

	     dbits:     qshift>4 ? 0 : 1,
	     dshift:    3,
	     dloc:      15,


	     // NB: Also useful as a way of embedding sel and doing sel
	     // specific contexts. Identical bar context. Eg 0<<15 or 1<<15.
	     sbits:     0,
	     sloc:      15,
	     do_stab:   0,
	     context:   (0<<15),

	     max_sym:   max_sym,
	     nsym:      nsym,

	     do_qmap:   do_qmap,
	     do_dedup:  0,
	     fixed_len: (q_lens.length == 1) ? 1 : 0,
	     do_sel:    0,
	     do_rev:    0,
	     do_pos:    1,
	     do_delta:  (qshift <= 4) ? 1 : 0,
	     do_qtab:   0,

	     // Override above with some attempt at using selectors
	     // when the q_dirs are specific and non-fixed.
	     qbits:     8+(qshift>4)-(qd_fixed==0),
	     sbits:     1,
	     sloc:      15-(qshift<=4), // read1 vs read2
	     do_stab:   1,
	     do_sel:    1,
	     
//	     // q4+dir: 7245769 with, 7353962 without. 1.5% saving
//	     qbits:     6,
//	     dbits:     2,
//	     dshift:    2,
//	     dloc:      13,
//	     sbits:     1,
//	     sloc:      15,
//	     do_stab:   1,
//	     do_sel:    1,

	     // with 20 bits of context, q40 = 31741545
	     // qbits 10, dbits 2, pbits 7, sbits 1
	    }]
}

function store_array(out, tab, size) {
    var i = 0; // index into tab
    var j = 0; // current value in tab[i]

    var tmp1 = new Array(size*2);
    var sz1 = 0;

    // First level of RLE.  Replace all runs of 'j' values
    // with run-lengths, including zeros for missing values.
    // Eg 0 1 2 2 2 3 3 3 4 4 4 5 5 5 5   7 7
    // to 1 1 3     3     3     4       0 2
    while (i < size) {
	// Length of j^{th} element
	var i_start = i
	while (i < size && tab[i] == j)
	    i++;
	var run_len = i - i_start

	// Encode run length to tmp array
	do {
	    var r = Math.min(255, run_len)
	    tmp1[sz1++] = r
	    run_len -= r
	} while (r == 255)
	j++;
    }

    // Second round of RLE on our tmp array, using a different
    // RLE algorithm.
    // Eg 1 1    3 3  3 4 0 2
    // to 1 1 +0 3 3 +1 4 0 2
    var last = -1
    var tmp2 = new Array(size*2)
    var sz2 = 0
    i = 0  // index into tmp1]
    // k is used size of tmp1[]
    while (i < sz1) {
	var curr = tmp1[i++];
	tmp2[sz2++] = curr
	if (curr == last) {
	    var i_start = i;
	    while (i < sz1 && tmp1[i] == last && i - i_start < 255)
		i++;
	    tmp2[sz2++] = i - i_start;
	} else {
	    last = curr
	}
    }

    // Append 2nd RLE, tmp2, to out.
    out.WriteData(tmp2, sz2)
}

				     

// q_lens is an array of quality lengths per record.
// (If they're all the same, just set one value.)
function encode_fqz_params(out, params, qhist, qtab, ptab, dtab, stab) {
    var dsqr = [
        0, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3,
        4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5,
        5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
        6, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7
    ]

    for (var i = 0; i < params.length; i++)
	stab[i] = i; // 1 parameter set per selector value
    for (; i < 256; i++)
	stab[i] = params.length-1;

    // Store global meta-data
    out.WriteByte(5);            // FQZ format number
    var gflags = ((params.length > 1) ? GFLAG_MULTI_PARAM : 0)
	       | ((params[0].do_stab) ? GFLAG_HAVE_STAB   : 0)
    out.WriteByte(gflags)

    if (gflags & GFLAG_MULTI_PARAM)
	out.WriteByte(params.length) // Number of parameter blocks.

    if (gflags & GFLAG_HAVE_STAB) {
	var max_sel = 1<<params[0].sbits;
	if (max_sel > 0) max_sel--;
	out.WriteByte(max_sel)
	store_array(out, stab, 256)
    }

    // Store per-param meta-data
    for (var p = 0; p < params.length; p++) {
	out.WriteUint16(params[p].context)
	out.WriteByte((params[p].do_qtab  ? FLAG_QTAB  : 0) |  // FLAG
		      (params[p].do_delta ? FLAG_DTAB  : 0) |
		      (params[p].do_pos   ? FLAG_PTAB  : 0) |
		      (params[p].do_qmap  ? FLAG_QMAP  : 0) |
		      (params[p].do_sel   ? FLAG_SEL   : 0) |
		      (params[p].fixed_len? FLAG_FLEN  : 0) |
		      (params[p].do_dedup ? FLAG_DEDUP : 0))
	if (params[p].do_qmap)
	    out.WriteByte(params[p].nsym)
	else
	    out.WriteByte(params[p].max_sym)
	out.WriteByte((params[p].qbits << 4) | (params[p].qshift))
	out.WriteByte((params[p].qloc  << 4) | (params[p].sloc))
	out.WriteByte((params[p].ploc  << 4) | (params[p].dloc))

	if (params[p].do_qmap) {
	    params[p].max_sym = params[p].nsym
	    var n = 0;
	    for (var i = 0; i < 256; i++) {
		if (qhist[p][i]) {
		    out.WriteByte(i)
		    qhist[p][i] = n++;
		}
	    }
	    // Ensure we have all matched input params
	    for (; n < params[p].nsym; n++)
		out.WriteByte(0)
	} else {
	    //params[p].nsym = 255;
	    for (var i = 0; i < 256; i++)
		qhist[p][i] = i; // NOP
	}

	if (params[p].qbits > 0) {
	    //	// Eg map 0-44 to a smaller range, to improve context usage.
	    //	// Makes q40 test set go from 33596471 to 33450075 (-0.4%)
	    //	params[p].do_qtab = 1;
	    //	for (var j = i = 0; i < params[p].max_sym; i++) {
	    //	    qtab[i]=j;
	    //	    if ((i%3)!=0 | i >= 28) j++
	    //	    console.log("qtab[",i,"]=",qtab[i]);
	    //	}
	    //	for (; i < 256; i++)
	    //	    qtab[i] = qtab[params[p].max_sym-1]

	    for (var i = 0; i < 256; i++)
		qtab[p][i] = i; // NOP for now

	    if (params[p].do_qtab)
		store_array(out, qtab[p], 256)
	}

	if (params[p].pbits > 0) {
	    for (var i = 0; i < 1024; i++)
		ptab[p][i] = Math.min((1<<params[p].pbits)-1, i >> params[p].pshift)

	    store_array(out, ptab[p], 1024)
	}

	if (params[p].dbits > 0) {
	    for (var i = 0; i < 256; i++)
		if (dsqr[i] > (1<<params[p].dbits) - 1)
		    dsqr[i] = (1<<params[p].dbits) - 1
	    for (var i = 0; i < 256; i++)
		dtab[p][i] = dsqr[Math.min(dsqr.length-1, i >> params[p].dshift)]

	    store_array(out, dtab[p], 256)
	}
    }

    return out
}

function encode_fqz(out, src, q_lens, q_dirs, params, qhist, qtab, ptab, dtab, stab) {
    //console.error("0:",params[0])
    //console.error("1:",params[1])

    var max_sel = 1<<params[0].sbits
    if (max_sel > 0) max_sel--
    var n_in = src.length

    // Create the models
    var max_sym = 0;
    for (var p = 0; p < params.length; p++)
	if (max_sym < params[p].max_sym)
	    max_sym = params[p].max_sym;

    var model_qual = new Array(1<<16)
    for (var i = 0; i < (1<<16); i++)
	model_qual[i] = new ByteModel(max_sym+1)

    var model_len = new Array(4)
    for (var i = 0; i < 4; i++)
	model_len[i] = new ByteModel(256)

    var model_rev    = new ByteModel(2)
    var model_dup    = new ByteModel(2)
    var model_sel    = new ByteModel(max_sel+1)

    // Note: our JavaScript encoder doesn't have a way for reversing
    // some quality strings, so we ignore do_rev for now.
    var rc = new RangeCoder(src)

    // The main encoding loop
    var p = 0; // remaining position along current record
    var i = 0; // index in src data
    var rec = 0;

    while (i < n_in) {
	if (p == 0) {
	    //var s = 0 // single non-mixed sample
	    var s = q_dirs[rec]
	    if (params[0].sbits > 0) {// FIXME: check All params[].do_stab / sbits must be identical
		//console.log("Ssel", s)
	        model_sel.ModelEncode(out, rc, s)
	    }
	    var x = stab[s]

	    // Reset contexts at the statr of each new record
	    var len = q_lens[Math.min(q_lens.length-1, rec++)]
	    if (params[x].fixed_len) {
		if (params[x].fixed_len > 0) { // First length
		    //console.log("Len", len)
		    model_len[0].ModelEncode(out, rc, len       & 0xff)
		    model_len[1].ModelEncode(out, rc, (len>>8)  & 0xff)
		    model_len[2].ModelEncode(out, rc, (len>>16) & 0xff)
		    model_len[3].ModelEncode(out, rc, (len>>24) & 0xff)
		    params[x].fixed_len = -1; // indicate we've stored it once
		}
	    } else {
		//console.log("len", len)
		model_len[0].ModelEncode(out, rc, len       & 0xff)
		model_len[1].ModelEncode(out, rc, (len>>8)  & 0xff)
		model_len[2].ModelEncode(out, rc, (len>>16) & 0xff)
		model_len[3].ModelEncode(out, rc, (len>>24) & 0xff)
	    }

	    if (params[x].do_dedup)
		process.exit(1) // FIXME

	    p = len
	    var delta = 0
	    //var last  = 0
	    var last  = params[x].context
	    var qlast = 0
	    var q1    = 0
	}

	// Encode current quality
	var q = src[i++]
	var Q = qhist[x][q]
	model_qual[last].ModelEncode(out, rc, Q)
	//console.log("Ctx",last,qhist[x][q])

	// Update contexts for next quality
	qlast = (qlast << params[x].qshift) + qtab[x][Q]
	last  = params[x].context
	last += (qlast & ((1<<params[x].qbits)-1)) << params[x].qloc

	// 46.6-48.6 billion cycles with ifs + "<< params[x].?loc" shifts
	// 47.3-47.3 billion cycles with ifs
	// 47.1-47.9 billion cycles without ifs
	if (params[x].pbits > 0)
	    last += ptab[x][Math.min(p, 1023)] << params[x].ploc

	if (params[x].dbits > 0) {
	    last += dtab[x][Math.min(delta, 255)] << params[x].dloc
	    delta += (q1 != Q) ? 1 : 0
	    q1 = Q
	}

	if (params[x].do_sel)
	    last += s << params[x].sloc

	last = (last & 0xffff)
	p--
    }

    rc.RangeFinishEncode(out)
    return out.buf.slice(0, out.pos)
}

function encode(src, q_lens, q_dirs) {
    var qhist = new Array(2)
    var qtab  = new Array(2)
    var ptab  = new Array(2)
    var dtab  = new Array(2)
    var stab  = new Array(256)

    for (var s = 0; s < 2; s++) {
        qhist[s] = new Array(256)
        qtab[s]  = new Array(256)
        ptab[s]  = new Array(1024) 
        dtab[s]  = new Array(256)
    }

    var out = new IOStream("", 0, src.length*1.1 + 100); // FIXME: guestimate worst case

    out.WriteUint7(src.length);
    var params = pick_fqz_params(src, q_lens, q_dirs, qhist)
    var out = encode_fqz_params(out, params, qhist, qtab, ptab, dtab, stab)
    return encode_fqz(out, src, q_lens, q_dirs, params, qhist, qtab, ptab, dtab, stab)
}

module.exports = { decode, encode }
