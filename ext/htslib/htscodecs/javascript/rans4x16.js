/*
 * Copyright (c) 2019,2020 Genome Research Ltd.
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

//----------------------------------------------------------------------
// rANS primitives itself
//
// RansGet* is decoder side

function RansGetCumulativeFreq(R, bits) {
    return R & ((1<<bits)-1)
}

function RansGetSymbolFromFreq(C, f) {
    // NOTE: Inefficient.
    // In practice we would implement this via a precomputed
    // lookup table C2S[f]; see RansBuildC2S below.
    var s = 0;
    while (f >= C[s+1])
	s++;

    //console.error(f, C, s)

    return s;
}

function RansBuildC2S(C, bits) {
    var max = 1<<bits
    var C2S = new Array(max);
    var s = 0;
    for (var f = 0; f < max; f++) {
	while (f >= C[s+1])
	    s++;
	C2S[f] = s;
    }
    return C2S;
}

function RansAdvanceStep(R, c, f, bits) {
    return f * (R >> bits) + (R & ((1<<bits)-1)) - c;
}

function RansRenorm(src, R) {
    if (R < (1<<15))
	R = (R << 16) + src.ReadUint16();

    return R;
}


// RanEnc* is for encoder
function RansEncInit() {
    return 1<<15;
}

function RansEncFlush(R, dst) {
    dst.WriteByteNeg((R >> 24) & 0xff);
    dst.WriteByteNeg((R >> 16) & 0xff);
    dst.WriteByteNeg((R >>  8) & 0xff);
    dst.WriteByteNeg((R >>  0) & 0xff);
}

function RansEncRenorm(R, dst, freq, scale_bits) {
    //var R_max = (((1 << 15) >> scale_bits) << 16) * freq;
    var R_max = (1 << (31-scale_bits)) * freq

    while (R >= R_max) {
	dst.WriteByteNeg((R>>8) & 0xff);
	dst.WriteByteNeg(R & 0xff);
	R >>= 16;
    }
    return R;
}

// Puts a symbol with frequency freq, cumulative freq start
// and total freq 1<<scale_bits.
//
// Note with static probabilities, /freq and %freq could be
// precomputed via multiplies and shifts.
function RansEncPut(R, dst, start, freq, scale_bits) {
    var scale = 1<<scale_bits;
    R = RansEncRenorm(R, dst, freq, scale_bits);
    R = (Math.floor(R / freq) << scale_bits) + (R % freq) + start;
    return R;
}


//----------------------------------------------------------------------
// Run length encoding
function EncodeRLE(src, N) {
    // Step 1: find which symbols benefit from RLE
    var L = new Array(256)
    for (var i = 0; i < 256; i++)
	L[i] = 0

    var last = -1
    for (var i = 0; i < src.length; i++) {
	L[src[i]] += src[i] == last ? 1 : -1
	last = src[i]
    }

    var nrle = 0;
    for (var i = 0; i < 256; i++)
	if (L[i] > 0)
	    nrle++

    if (!nrle) {
	// Format cannot cope with zero RLE symbols, so pick one!
	nrle = 1;
	L[0] = 1;
    }

    // Start meta-data as list of symbols to RLE
    var meta = new IOStream("", 0, nrle+1 + src.length)
    meta.WriteByte(nrle)
    for (var i = 0; i < 256; i++)
	if (L[i] > 0)
	    meta.WriteByte(i)

    // Step 2: Now apply RLE itself
    var data = new Buffer.allocUnsafe(src.length)
    var dpos = 0
    for (var i = 0; i < src.length; i++) {
	data[dpos++] = src[i]
	if (L[src[i]] > 0) {
	    last = src[i];
	    var run = 0;
	    while (i+run+1 < src.length && src[i+run+1] == last)
		run++;
	    meta.WriteUint7(run)
	    i += run
	}
    }

    // Compress the meta-data
    var cmeta = RansEncode0(meta.buf.slice(0, meta.pos), N)
    var hdr = new IOStream("", 0, 16)
    hdr.WriteUint7(meta.pos*2)   // Uncompressed meta-data length + compressed-bit-flag(0)
    hdr.WriteUint7(dpos)         // Length of RLE encoded data
    hdr.WriteUint7(cmeta.length) // Compressed meta-data length
    var meta = Buffer.concat([hdr.buf.slice(0,hdr.pos), cmeta])

    return [meta, data.slice(0, dpos)]
}

function DecodeRLEMeta(src, N) {
    var u_meta_len = src.ReadUint7()
    var rle_len = src.ReadUint7()

    // Decode RLE lengths
    if (u_meta_len & 1) {
	var rle_meta = src.ReadData((u_meta_len-1)/2)
    } else {
	var comp_meta_len = src.ReadUint7()
	var rle_meta = src.ReadData(comp_meta_len)
	rle_meta = RansDecode0(new IOStream(rle_meta), u_meta_len/2, N)
    }

    // Decode list of symbols for which RLE lengths are applied
    var rle_meta = new IOStream(rle_meta)
    var L = new Array(256)
    var n = rle_meta.ReadByte()
    if (n == 0)
	n = 256;
    for (var i = 0; i < n; i++)
	L[rle_meta.ReadByte()] = 1

    return [L, rle_meta, rle_len]
}

function DecodeRLE(buf, L, rle_meta, len) {
    var src = new IOStream(buf);

    var out = new Buffer.allocUnsafe(len)

    // Expand up buf+meta to out; i = buf index, j = out index
    var j = 0;
    for (var i = 0; j < len; i++) {
	var sym = buf[i];
	if (L[sym]) {
	    var run = rle_meta.ReadUint7()
	    for (var r = 0; r <= run; r++)
		out[j++] = sym
	} else {
	    out[j++] = sym
	}
    }

    return out
}

//----------------------------------------------------------------------
// Bit packing

function EncodePack(src) {
    // Step 1: identify number of distinct symbols
    var F = new Array(256)
    for (var i = 0; i < 256; i++)
	F[i] = 0

    for (var i = 0; i < src.length; i++)
	F[src[i]]++

    var P = new Array(256)
    var nsym = 0;
    for (var i = 0; i < 256; i++)
	if (F[i] > 0)
	    P[i] = nsym++

    if (nsym > 16) {
	//console.error("Too many symbols to pack:",nsym)
	return
    }


    // Pack data
    if (nsym <= 1) {
	// Constant
	var data = new Buffer.allocUnsafe(0)
    }

    else if (nsym <= 2) {
	// 1 bit per value
	var data = new Buffer.allocUnsafe(Math.ceil(src.length/8))
	var j = -1
	for (i = 0; i < src.length; i++) {
	    if (i % 8 == 0)
		data[++j] = 0
	    data[j] += P[src[i]] << (i % 8)
	}
    }

    else if (nsym <= 4) {
	// 2 bits per value
	var data = new Buffer.allocUnsafe(Math.ceil(src.length/4))
	var j = -1
	for (i = 0; i < src.length; i++) {
	    if (i % 4 == 0)
		data[++j] = 0
	    data[j] += P[src[i]] << ((i % 4) * 2)
	}
    }

    else {
	// 4 bits per value
	var data = new Buffer.allocUnsafe(Math.ceil(src.length/2))
	var j = -1
	for (i = 0; i < src.length; i++) {
	    if (i % 2 == 0)
		data[++j] = 0
	    data[j] += P[src[i]] << ((i % 2) * 4)
	}
    }


    // Produce pack meta-data
    var meta = new IOStream("", 0, nsym+5)
    meta.WriteByte(nsym)
    var j = 0
    for (var i = 0; i < 256; i++) {
	if (F[i] > 0) {
	    F[i] = j++;
	    meta.WriteByte(i)
	}
    }
    meta.WriteUint7(data.length)

    return [meta.buf.slice(0, meta.pos), data]
}


// Pack meta data is the number and value of distinct symbols plus
// the length of the packed byte stream.
function DecodePackMeta(src) {
    var nsym = src.ReadByte()
    var P = new Array(nsym)

    for (var i = 0; i < nsym; i++)
	P[i] = src.ReadByte()

    var len = src.ReadUint7()

    return [P, nsym, len]
}

// Extract bits from src producing output of length len.
// Nsym is number of distinct symbols used.
function DecodePack(data, P, nsym, len) {
    var out = new Buffer.allocUnsafe(len)
    var j = 0;

    // Constant value
    if (nsym <= 1) {
	for (var i = 0; i < len; i++)
	    out[i] = P[0]
    }

    // 1 bit per value
    else if (nsym <= 2) {
	for (i = 0; i < len; i++) {
	    if (i % 8 == 0)
		var v = data[j++];

	    out[i] = P[v & 1]
	    v >>= 1
	}
    }

    // 2 bits per value
    else if (nsym <= 4) {
	for (i = 0; i < len; i++) {
	    if (i % 4 == 0)
		var v = data[j++];

	    out[i] = P[v & 3]
	    v >>= 2
	}
    }

    // 4 bits per value
    else if (nsym <= 16) {
	for (i = 0; i < len; i++) {
	    if (i % 2 == 0)
		var v = data[j++];

	    out[i] = P[v & 15]
	    v >>= 4
	}
    }

    return out
}


//----------------------------------------------------------------------
// 4 way interleaving.
// This is simply 4 rANS streams interleaved to form bytes 0,4,8...,
// 1,5,9..., 2,6,10... and 3,7,11...
//
// It works well when the distributions differ for each of the 4 bytes,
// for example when compressing a series of 32-bit integers.
//
// Maybe make this more general purpose of X* where we specify the stripe
// size instead of fixing it at 4?
function RansEncodeStripe(hdr, src, N) {
    if (N == 0)
	N = 4; // old default

    // Split into multiple streams
    var part = new Array(N)
    var ulen = new Array(N)
    for (var s = 0; s < N; s++) {
	ulen[s] = Math.floor(src.length / N) + ((src.length % N) > s);
	part[s] = new Array(ulen[s])
    }

    for (var x = 0, i = 0; i < src.length; i+=N, x++) {
	for (var j = 0; j < N; j++)
	    if (x < part[j].length)
		part[j][x] = src[i+j]
    }

    // Compress each part
    var comp = new Array(N)
    var total = 0
    for (var s = 0; s < N; s++) {
	// Example: try O0 and O1 and choose best
	var comp0 = encode(part[s], 0)
	var comp1 = encode(part[s], 1)
	comp[s] = (comp1.length < comp0.length) ? comp1 : comp0
	total += comp[s].length
    }

    // Serialise
    var out = new IOStream("", 0, total+5*N+1)
    out.WriteByte(N)
    for (var s = 0; s < N; s++)
	out.WriteUint7(comp[s].length)

    for (var s = 0; s < N; s++)
	out.WriteData(comp[s], comp[s].length)

    return out.buf.slice(0, out.buf.pos)
}

function RansDecodeStripe(src, len) {
    var N = src.ReadByte()

    // Retrieve lengths
    var clen = new Array(N)
    var ulen = new Array(N)
    for (var j = 0; j < N; j++)
	clen[j] = src.ReadUint7()

    // Decode streams
    var T = new Array(N);
    for (var j = 0; j < N; j++) {
	ulen[j] = Math.floor(len / N) + ((len % N) > j)
	T[j] = RansDecodeStream(src, ulen[j])
    }

    // Transpose
    var out = new Buffer.allocUnsafe(len)
    for (var j = 0; j < N; j++) {
	for (var i = 0; i < ulen[j]; i++) {
	    out[i*N + j] = T[j][i];
	}
    }

    return out;
}


//----------------------------------------------------------------------
// Main rANS entry function: decodes a compressed src and
// returns the uncompressed buffer.
function decode(src) {
    var stream = new IOStream(src)
    return RansDecodeStream(stream, 0)
}

function RansDecodeStream(stream, n_out) {
    var format = stream.ReadByte();
    var order  = format & 1
    var x32    = format & 4
    var stripe = format & 8
    var nosz   = format & 16
    var cat    = format & 32
    var rle    = format & 64
    var pack   = format & 128

    var Nway = x32 ?32 : 4

    if (!nosz)
	n_out = stream.ReadUint7();

    // N-way interleaving
    if (stripe)
	return RansDecodeStripe(stream, n_out)

    // Bit packing
    if (pack) {
	var pack_len = n_out
	var [P, nsym, n_out] = DecodePackMeta(stream)
    }

    // Run length encoding
    if (rle) {
	var rle_len = n_out
	var [L, rle_meta, n_out] = DecodeRLEMeta(stream, Nway)
    }

    // Uncompress data (all, packed or run literals)
    if (cat)
	var buf = stream.ReadData(n_out)
    else if (order == 0)
	var buf = RansDecode0(stream, n_out, Nway)
    else
	var buf = RansDecode1(stream, n_out, Nway)

    // Apply expansion transforms
    if (rle)
	buf = DecodeRLE(buf, L, rle_meta, rle_len)

    if (pack)
	buf = DecodePack(buf, P, nsym, pack_len)

    return buf
}

function encode(src, format) {
    var hdr = new IOStream("", 0, 10);
    hdr.WriteByte(format);

    var order = format & 1
    var x32   = format & 4
    var stripe= format & 8
    var nosz  = format & 16
    var cat   = format & 32
    var rle   = format & 64
    var pack  = format & 128

    var Nway  = x32 ? 32 : 4; // interleaving amount
    var N     = format>>8;    // stripe size

    if (!nosz)
	hdr.WriteUint7(src.length);

    if (stripe)
	return Buffer.concat([hdr.buf.slice(0, hdr.pos), RansEncodeStripe(hdr, src, N)])

    var pack_meta = new Buffer.alloc(0)
    if (pack)
	[pack_meta, src] = EncodePack(src)

    var rle_meta = new Buffer.alloc(0)
    if (rle)
	[rle_meta, src] = EncodeRLE(src, Nway)

    if (src.length < 4 && order == 1) {
	// Protect against short order-1 data due to RLE/Pack
	order = 0
	hdr.buf[0] &= ~1
    }

    if (cat)
	var comp = src
    else if (order == 0)
	var comp = RansEncode0(src, Nway)
    else
	var comp = RansEncode1(src, Nway)

    return Buffer.concat([hdr.buf.slice(0,hdr.pos), pack_meta, rle_meta, comp])
}

//----------------------------------------------------------------------
// Order-0 decoder

function ReadAlphabet(src) {
    var A = new Array(256)
    for (var i = 0; i < 256; i++)
	A[i] = 0;

    var rle = 0
    var sym = src.ReadByte()
    var last_sym = sym

    do {
	A[sym] = 1;
	if (rle > 0) {
	    rle--
	    sym++
	} else {
	    sym = src.ReadByte()
	    if (sym == last_sym+1)
		rle = src.ReadByte()
	}
	last_sym = sym
    } while (sym != 0)

    return A
}

// Decode a single table of order-0 frequences,
// filling out the F and C arrays.
function ReadFrequencies0(src, F, C) {
    // Initialise; not in the specification - implicit?
    for (var i = 0; i < 256; i++)
	F[i] = 0;

    // Fetch alphabet
    var A = ReadAlphabet(src);

    // Fetch frequencies for the symbols listed in our alphabet
    for (var i = 0; i < 256; i++) {
	if (A[i] > 0)
	    F[i] = src.ReadUint7()
    }

    NormaliseFrequencies0_Shift(F, 12)

    // Compute C[] from F[]
    C[0] = 0;
    for (var i = 0; i <= 255; i++)
	C[i+1] = C[i] + F[i];
}

function RansDecode0(src, nbytes, N) {
    // Decode frequencies
    var F = new Array(256);
    var C = new Array(256);
    ReadFrequencies0(src, F, C);

    // Fast lookup to avoid slow RansGetSymbolFromFreq
    var C2S = RansBuildC2S(C, 12);

    // Initialise rANS state
    var R = new Array(N);
    for (var i = 0; i < N; i++)
	R[i] = src.ReadUint32();

    // Main decode loop
    var output = new Buffer.allocUnsafe(nbytes);
    for (var i = 0; i < nbytes; i++) {
	var ix = i&(N-1); // equiv to i%N as N is power of 2
	var f = RansGetCumulativeFreq(R[ix], 12);
	var s = C2S[f]; // Equiv to RansGetSymbolFromFreq(C, f);

	output[i] = s;
	R[ix] = RansAdvanceStep(R[ix], C[s], F[s], 12);
	R[ix] = RansRenorm(src, R[ix]);
    }

    // Main decode loop
    return output;
}

//----------------------------------------------------------------------
// Order-0 encoder

function BuildFrequencies0(src, F) {
    for (var i = 0; i < 256; i++)
	F[i] = 0;

    for (var i = 0; i < src.length; i++)
	F[src[i]]++;
}

function NormaliseFrequencies0(F, bits) {
    // Compute total
    var tot = 0;
    for (var i = 0; i < 256; i++)
	tot += F[i];

    // Scale total of frequencies to max
    const max = (1<<bits);
    var scale = max / tot;
    do {
	var max_val = 0;
	var max_idx = 0;
	var renorm = 0;
	tot = 0;
	for (var i = 0; i < 256; i++) {
	    if (F[i] == 0)
		continue

	    if (max_val < F[i]) {
		max_val = F[i]
		max_idx = i
	    }

	    F[i] = Math.floor(F[i] * scale);
	    if (F[i] == 0)
		F[i] = 1;

	    tot += F[i];
	}

	// Adjust new tot to ensure it matches.
	if (tot < max) {
	    // Too low, boost the most common symbol
	    F[max_idx] += max-tot;
	} else if (tot-max < F[max_idx]/2 && F[max_idx] > 2) {
	    // Too high, reduce the common symbol
	    F[max_idx] -= tot-max;
	} else if (tot != max) {
	    // Much too high, fudge scale and try again.
	    scale = max / tot;
	    renorm = 1;
	}
    } while (renorm)
}

function NormaliseFrequencies0_Shift(F, bits) {
    // Compute total and number of bits to shift by
    var tot = 0;
    for (var i = 0; i < 256; i++)
	tot += F[i];

    if (tot == 0 || tot == (1<<bits))
	return

    var shift = 0;
    while (tot < (1<<bits)) {
	tot *= 2;
	shift++;
    }

    // Scale total of frequencies to (1<<bits)
    for (var i = 0; i < 256; i++)
	F[i] <<= shift;
}

function WriteAlphabet(out, F) {
    var rle = 0;
    for (var i = 0; i < 256; i++) {
	if (!F[i])
	    continue

	if (rle > 0)
	    rle--
	else {
	    out.WriteByte(i)

	    if (i > 0 && F[i-1] > 0) {
		// We've encoded two symbol frequencies in a row.
		// How many more are there?  Store that count so
		// we can avoid writing consecutive symbols.
		for (rle = i+1; rle<256 && F[rle]; rle++)
		    ;
		rle -= i+1;

		out.WriteByte(rle);
	    }
	}
    }
    out.WriteByte(0)
}

function WriteFrequencies0(out, F) {
    WriteAlphabet(out, F)

    for (var i = 0; i < 256; i++) {
	if (F[i])
	    out.WriteUint7(F[i])
    }
}

function RansEncode0(src, N) {
    const nbytes = src.length;
    var output = new IOStream("", 0, 257*3+9);

    // Compute frequencies
    var F = new Array(256)
    BuildFrequencies0(src, F)
    var bit_size = Math.ceil(Math.log2(nbytes));
    if (bit_size > 12)
	bit_size = 12;
    NormaliseFrequencies0(F, bit_size);
    WriteFrequencies0(output, F);
    NormaliseFrequencies0(F, 12);

    // Compute cumulative frequencies
    var C = new Array(256)
    C[0] = 0;
    for (var i = 1; i < 256; i++)
	C[i] = C[i-1] + F[i-1];

    // Initialise rANS state
    var R = new Array(N);
    for (var i = 0; i < N; i++)
	R[i] = RansEncInit();

    // Allow expansion room if trying to compress random data.
    var rans_out = new IOStream("", (nbytes*1.05+100)>>0, (nbytes*1.05+100)>>0);

    // Main encode loop
    for (var i = nbytes-1; i >= 0; i--)
	R[i%N] = RansEncPut(R[i%N], rans_out, C[src[i]], F[src[i]], 12);

    for (var i = N-1; i >= 0; i--)
	RansEncFlush(R[i], rans_out);

    // Stitch blocks together into final output buffer
    //console.error("pos=",rans_out.pos, " len=",rans_out.length)
    //console.error(rans_out.buf.slice(rans_out.pos, rans_out.length))
    return Buffer.concat([output.buf.slice(0, output.pos),
			  rans_out.buf.slice(rans_out.pos, rans_out.length)],
			 output.pos + rans_out.length - rans_out.pos);
}

//----------------------------------------------------------------------
// Order-1 decoder

// Decode a table of order-1 frequences,
// filling out the F and C arrays.
function ReadFrequencies1(src, F, C, shift) {
    // Initialise; not in the specification - implicit?
    for (var i = 0; i < 256; i++) {
	F[i] = new Array(256);
	C[i] = new Array(256);
	for (var j = 0; j < 256; j++)
	    F[i][j] = 0;
    }

    // Fetch alphabet
    var A = ReadAlphabet(src);

    // Read F[]
    for (var i = 0; i < 256; i++) {
	if (!A[i])
	    continue

	var run = 0;
	for (var j = 0; j < 256; j++) {
	    if (!A[j])
		continue

	    if (run > 0) {
		run--
	    } else {
		F[i][j] = src.ReadUint7();
		if (F[i][j] == 0)
		    run = src.ReadByte();
	    }
	}

	NormaliseFrequencies0_Shift(F[i], shift)

	// Compute C[] from F[]
	C[i][0] = 0;
	for (var j = 0; j < 256; j++)
	    C[i][j+1] = C[i][j] + F[i][j];
    }
}

function RansDecode1(src, nbytes, N) {
    // FIXME: this bit is missing from the RansDecode0 pseudocode.

    var comp = src.ReadByte();
    var shift = comp >> 4;

    var freq_src = src
    if (comp & 1) {
	var ulen = src.ReadUint7()
	var clen = src.ReadUint7()
	var comp = new IOStream(src.ReadData(clen))
	var freq_src = new IOStream(RansDecode0(comp, ulen, 4));
    }

    // Decode frequencies
    var F = new Array(256);
    var C = new Array(256);
    ReadFrequencies1(freq_src, F, C, shift);

    // Fast lookup to avoid slow RansGetSymbolFromFreq
    var C2S = new Array(256);
    for (var i = 0; i < 256; i++)
	// Could do only for symbols in alphabet?
	C2S[i] = RansBuildC2S(C[i], shift);

    // Initialise rANS state
    var R = new Array(N);
    var L = new Array(N);
    for (var j = 0; j < N; j++) {
	R[j] = src.ReadUint32();
	L[j] = 0;
    }

    // Main decode loop
    var output = new Buffer.allocUnsafe(nbytes);
    var nbytesx = Math.floor(nbytes/N);
    for (var i = 0; i < nbytesx; i++) {
	for (var j = 0; j < N; j++) {
	    var f = RansGetCumulativeFreq(R[j], shift);

	    //var s = RansGetSymbolFromFreq(C[L[j]], f);
	    var s = C2S[L[j]][f]; // Precomputed version of above

	    output[i+j*nbytesx] = s;
	    R[j] = RansAdvanceStep(R[j], C[L[j]][s], F[L[j]][s], shift);
	    R[j] = RansRenorm(src, R[j]);
	    L[j] = s;
	}
    }

    // Now deal with the remainder if buffer size is not a multiple of N,
    // using the last rANS state exclusively.  (It'd have been nice to have
    // designed this to just act as if we kept going with a bail out.)
    i = N*i;
    while (i < nbytes) {
	var f = RansGetCumulativeFreq(R[N-1], shift);
	var s = RansGetSymbolFromFreq(C[L[N-1]], f);
	output[i++] = s;
	R[N-1] = RansAdvanceStep(R[N-1], C[L[N-1]][s], F[L[N-1]][s], shift);
	R[N-1] = RansRenorm(src, R[N-1]);
	L[N-1] = s;
    }

    return output;
}

//----------------------------------------------------------------------
// Order-1 encoder

function BuildFrequencies1(src, F, F0, N) {
    for (var i = 0; i < 256; i++) {
	F0[i] = 0;
	for (var j = 0; j < 256; j++)
	    F[i][j] = 0;
    }

    var last = 0;
    for (var i = 0; i < src.length; i++) {
	F0[last]++;
	F[last][src[i]]++;
	last = src[i];
    }
    F0[last]++;

    // Also accept we'll be starting at N points, not just byte 0
    for (var i = 1; i < N; i++)
	F[0][src[i * Math.floor(src.length/N)]]++;
    F0[0] += N-1;
}

function NormaliseFrequencies1(F, F0, shift) {

    for (var i = 0; i < 256; i++) {
	if (!F0[i])
	    continue;

	var bit_size = Math.ceil(Math.log2(F0[i]));
	if (bit_size > shift)
	    bit_size = shift;

	NormaliseFrequencies0(F[i], bit_size)
    }
}

function NormaliseFrequencies1_Shift(F, F0, shift) {
    for (var i = 0; i < 256; i++)
	if (F0[i])
	    NormaliseFrequencies0_Shift(F[i], shift)
}

function WriteFrequencies1(out, F, F0) {
    WriteAlphabet(out, F0)

    for (var i = 0; i < 256; i++) {
	if (!F0[i])
	    continue

	var run = 0
	for (var j = 0; j < 256; j++) {
	    if (!F0[j])
		continue

	    if (run) {
		run--
	    } else {
		out.WriteUint7(F[i][j])

		if (!F[i][j]) {
		    // Count how many more zero-freqs we have
		    for (var k = j+1; k < 256; k++) {
			if (!F0[k])
			    continue

			if (F[i][k] == 0)
			    run++
			else
			    break
		    }
		    out.WriteByte(run)
		}
	    }
	}
    }
}

function RansEncode1(src, N) {
    const nbytes = src.length;
    var output = new IOStream("", 0, 257*257*3+9);

    // Compute frequencies
    var F0 = new Array(256)
    var F = new Array(256)
    var C = new Array(256)
    for (var i = 0; i < 256; i++) {
	F[i] = new Array(256);
	C[i] = new Array(256);
    }

    // Frequency precision
    var shift = 12;

    BuildFrequencies1(src, F, F0, N)
    NormaliseFrequencies1(F, F0, shift);

    // Store frequencies, possibly compressed
    var freq = new IOStream("", 0, 257*257*3+9);

    WriteFrequencies1(freq, F, F0);

    var cfreq = RansEncode0(freq.buf.slice(0, freq.pos), 4)
    if (cfreq.length < freq.pos) {
	output.WriteByte(1 | (shift<<4));
	output.WriteUint7(freq.pos)
	output.WriteUint7(cfreq.length)
	output.WriteData(cfreq, cfreq.length);
    } else {
	output.WriteByte(0 | (shift<<4));
	output.WriteData(freq.buf, freq.pos);
    }

    // Normalise and compute cumulative frequencies
    NormaliseFrequencies1_Shift(F, F0, shift);
    for (var i = 0; i < 256; i++) {
	if (!F0[i])
	    continue;

	C[i][0] = 0;
	for (var j = 1; j < 256; j++)
	    C[i][j] = C[i][j-1] + F[i][j-1];
    }

    // Initialise rANS state
    var R = new Array(N);
    var L = new Array(N);
    for (var j = 0; j < N; j++) {
	R[j] = RansEncInit();
	L[j] = 0;
    }
    var rans_out = new IOStream("", (nbytes*1.05+100)>>0, (nbytes*1.05+100)>>0);

    // We have N rans codecs running in parallel on its own 1/Nth of buffer
    var nbytesx = Math.floor(nbytes/N);
    var idx = new Array(N);
    var last = new Array(N)
    for (var j = 0; j < N; j++) {
	idx[j] = (j+1)*nbytesx - 2;
	last[j] = src[idx[j]+1]
    }

    // Deal with the remainder if not a multiple of N
    last[N-1] = src[nbytes-1];
    for (var i = nbytes-2; i > N*nbytesx-2; i--) {
	R[N-1] = RansEncPut(R[N-1], rans_out, C[src[i]][last[N-1]], F[src[i]][last[N-1]], shift);
	last[N-1] = src[i];
    }

    // Main encode loop
    var s = new Array(N);
    while (idx[0] >= 0) {
	// Separate loop as a speed optimisation
	for (var j = N-1; j >= 0; j--)
	    s[j] = src[idx[j]]
	for (var j = N-1; j >= 0; j--) {
	    R[j] = RansEncPut(R[j], rans_out, C[s[j]][last[j]], F[s[j]][last[j]], shift);
	    last[j] = s[j];
	    idx[j]--;
	}
    }

    for (var j = N-1; j >= 0; j--) {
        R[j] = RansEncPut(R[j], rans_out, C[0][last[j]], F[0][last[j]], shift)
    }

    for (var i = N-1; i >= 0; i--)
	RansEncFlush(R[i], rans_out);

    // Stitch blocks together into final output buffer
    return Buffer.concat([output.buf.slice(0, output.pos),
			  rans_out.buf.slice(rans_out.pos, rans_out.length)],
			 output.pos + rans_out.length - rans_out.pos);
}

module.exports = { decode, encode }
