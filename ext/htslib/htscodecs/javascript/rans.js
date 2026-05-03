/*
 * Copyright (c) 2019-2020 Genome Research Ltd.
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

function RansGetCumulativeFreq(R) {
    return R & 0xfff;
}

function RansGetSymbolFromFreq(C, f) {
    // NOTE: Inefficient.
    // In practice we would implement this via a precomputed
    // lookup table C2S[f]; see RansBuildC2S below.
    var s = 0;
    while (f >= C[s+1])
	s++;

    return s;
}

function RansBuildC2S(C) {
    var C2S = new Array(0x1000);
    var s = 0;
    for (var f = 0; f < 0x1000; f++) {
	while (f >= C[s+1])
	    s++;
	C2S[f] = s;
    }
    return C2S;
}

function RansAdvanceStep(R, c, f) {
    return f * (R >> 12) + (R & 0xfff) - c;
}

function RansRenorm(src, R) {
    while (R < (1<<23))
	R = (R << 8) + src.ReadByte();

    return R;
}


// RanEnc* is for encoder
function RansEncInit() {
    return 1<<23;
}

function RansEncFlush(R, dst) {
    dst.WriteByteNeg((R >> 24) & 0xff);
    dst.WriteByteNeg((R >> 16) & 0xff);
    dst.WriteByteNeg((R >>  8) & 0xff);
    dst.WriteByteNeg((R >>  0) & 0xff);
}

function RansEncRenorm(R, dst, freq, scale_bits) {
    var R_max = (((1 << 23) >> scale_bits) << 8) * freq;

    while (R >= R_max) {
	dst.WriteByteNeg(R & 0xff);
	R >>= 8;
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
// Main rANS entry function: decodes a compressed src and
// returns the uncompressed buffer.
function decode(src) {
    var stream = new IOStream(src);
    var order = stream.ReadByte();
    var n_in  = stream.ReadUint32();
    var n_out = stream.ReadUint32();

    if (order == 0) {
	return RansDecode0(stream, n_out)
    } else {
	return RansDecode1(stream, n_out)
    }
}

function encode(src, order) {
    //var stream = new IOStream(src);
    //var n_in  = stream.ReadUint32();
    //var n_out = stream.ReadUint32();

    if (order == 0) {
	return RansEncode0(src)
    } else {
	return RansEncode1(src)
    }
}

//----------------------------------------------------------------------
// Order-0 decoder

// Decode a single table of order-0 frequences,
// filling out the F and C arrays.
function ReadFrequencies0(src, F, C) {
    // Initialise; not in the specification - implicit?
    for (var i = 0; i < 256; i++)
	F[i] = 0;

    var sym = src.ReadByte();
    var last_sym = sym;
    var rle = 0;

    // Read F[]
    do {
	var f = src.ReadITF8();
	F[sym] = f;
	if (rle > 0) {
	    rle--;
	    sym++;
	} else {
	    sym = src.ReadByte();
	    if (sym == last_sym+1)
		rle = src.ReadByte();
	}
	last_sym = sym;
    } while (sym != 0);

    // Compute C[] from F[]
    C[0] = 0;
    for (var i = 0; i <= 255; i++)
	C[i+1] = C[i] + F[i];
}

function RansDecode0(src, nbytes) {
    // Decode frequencies
    var F = new Array(256);
    var C = new Array(256);
    ReadFrequencies0(src, F, C);

    // Fast lookup to avoid slow RansGetSymbolFromFreq
    var C2S = RansBuildC2S(C);

    // Initialise rANS state
    var R = new Array(4);
    for (var i = 0; i < 4; i++)
	R[i] = src.ReadUint32();

    // Main decode loop
    var output = new Buffer.allocUnsafe(nbytes);
    for (var i = 0; i < nbytes; i++) {
	var i4 = i%4;
	var f = RansGetCumulativeFreq(R[i4]);
	var s = C2S[f]; // Equiv to RansGetSymbolFromFreq(C, f);

	output[i] = s;
	R[i4] = RansAdvanceStep(R[i4], C[s], F[s]);
	R[i4] = RansRenorm(src, R[i4]);
    }

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

function NormaliseFrequencies0(F) {
    // Compute total
    var tot = 0;
    for (var i = 0; i < 256; i++)
	tot += F[i];

    // Scale total of frequencies to max
    const max = (1<<12);
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
	    scale = scale * 0.99
	    renorm = 1;
	}
    } while (renorm)
}

function WriteFrequencies0(out, F) {
    var rle = 0;
    for (var i = 0; i < 256; i++) {
	if (!F[i])
	    continue

	// Output Symbol if needed and Frequency
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

	out.WriteITF8(F[i])
    }
    out.WriteByte(0);
}

function RansEncode0(src) {
    const nbytes = src.length
    var output = new IOStream("", 0, 257*3+9);

    output.WriteByte(0);   // Order 0
    output.WriteUint32(0); // compressed size: correct later
    output.WriteUint32(0); // uncompressed size: correct later

    // Compute frequencies
    var F = new Array(256)
    BuildFrequencies0(src, F)
    NormaliseFrequencies0(F);
    WriteFrequencies0(output, F);

    // Compute cumulative frequencies
    var C = new Array(256)
    C[0] = 0;
    for (var i = 1; i < 256; i++)
	C[i] = C[i-1] + F[i-1];

    // Initialise rANS state
    var R = new Array(4);
    for (var i = 0; i < 4; i++)
	R[i] = RansEncInit();

    var alloc = Math.floor(nbytes*1.05+100)
    var rans_out = new IOStream("", alloc, alloc)

    // Main encode loop
    for (var i = nbytes-1; i >= 0; i--)
	R[i%4] = RansEncPut(R[i%4], rans_out, C[src[i]], F[src[i]], 12);

    for (var i = 3; i >= 0; i--)
	RansEncFlush(R[i], rans_out);

    // Stitch blocks together into final output buffer
    var freq_tab = output.pos
    output.buf.writeInt32LE(freq_tab-9 + (rans_out.length - rans_out.pos), 1);
    output.buf.writeInt32LE(nbytes, 5);

    return Buffer.concat([output.buf.slice(0, output.pos),
			  rans_out.buf.slice(rans_out.pos, rans_out.length)],
			 output.pos + rans_out.length - rans_out.pos);
}

//----------------------------------------------------------------------
// Order-1 decoder

// Decode a table of order-1 frequences,
// filling out the F and C arrays.
function ReadFrequencies1(src, F, C) {
    // Initialise; not in the specification - implicit?
    for (var i = 0; i < 256; i++) {
	F[i] = new Array(256);
	C[i] = new Array(256);
	for (var j = 0; j < 256; j++)
	    F[i][j] = 0;
    }

    var sym = src.ReadByte();
    var last_sym = sym;
    var rle = 0;

    // Read F[]
    do {
	ReadFrequencies0(src, F[sym], C[sym]);

	if (rle > 0) {
	    rle--;
	    sym++;
	} else {
	    sym = src.ReadByte();
	    if (sym == last_sym+1)
		rle = src.ReadByte();
	}
	last_sym = sym;
    } while (sym != 0);
}

function RansDecode1(src, nbytes) {
    // Decode frequencies
    var F = new Array(256);
    var C = new Array(256);
    ReadFrequencies1(src, F, C);

    // Fast lookup to avoid slow RansGetSymbolFromFreq
    var C2S = new Array(256);
    for (var i = 0; i < 256; i++)
	C2S[i] = RansBuildC2S(C[i]);

    // Initialise rANS state
    var R = new Array(4);
    var L = new Array(4);
    for (var j = 0; j < 4; j++) {
	R[j] = src.ReadUint32();
	L[j] = 0;
    }

    // Main decode loop
    var output = new Buffer.allocUnsafe(nbytes);
    var nbytes4 = Math.floor(nbytes/4);
    for (var i = 0; i < nbytes4; i++) {
	for (var j = 0; j < 4; j++) {
	    var f = RansGetCumulativeFreq(R[j]);

	    //var s = RansGetSymbolFromFreq(C[L[j]], f);
	    var s = C2S[L[j]][f]; // Precomputed version of above

	    output[i+j*nbytes4] = s;
	    R[j] = RansAdvanceStep(R[j], C[L[j]][s], F[L[j]][s]);
	    R[j] = RansRenorm(src, R[j]);
	    L[j] = s;
	}
    }

    // Now deal with the remainder if buffer size is not a multiple of 4,
    // using rANS state 3 exclusively.  (It'd have been nice to have
    // designed this to just act as if we kept going with a bail out.)
    i = 4*i;
    while (i < nbytes) {
	var f = RansGetCumulativeFreq(R[3]);
	var s = RansGetSymbolFromFreq(C[L[3]], f);
	output[i++] = s;
	R[3] = RansAdvanceStep(R[3], C[L[3]][s], F[L[3]][s]);
	R[3] = RansRenorm(src, R[3]);
	L[3] = s;
    }

    return output;
}

//----------------------------------------------------------------------
// Order-1 encoder

function BuildFrequencies1(src, F, F0) {
    for (var i = 0; i < 256; i++) {
	F0[i] = 0;
	for (var j = 0; j < 256; j++)
	    F[i][j] = 0;
    }

    var last = 0;
    for (var i = 0; i < src.length; i++) {
	F0[src[i]]++;
	F[last][src[i]]++;
	//F[last][src[i]]++;
	last = src[i];
    }

    // Also accept we'll be starting at 4 points, not just byte 0
    F[0][src[1*(src.length >> 2)]]++;
    F[0][src[2*(src.length >> 2)]]++;
    F[0][src[3*(src.length >> 2)]]++;
    F0[0] += 3;
}

function NormaliseFrequencies1(F, F0) {
    for (var i = 0; i < 256; i++)
	if (F0[i])
	    NormaliseFrequencies0(F[i])
}

function WriteFrequencies1(out, F, F0) {
    var rle = 0;
    var last_sym = 0;

    for (var i = 0; i < 256; i++) {
	if (!F0[i])
	    continue

	// Output Symbol if needed and Frequency
	if (rle > 0)
	    rle--
	else {
	    out.WriteByte(i)

	    if (i > 0 && F0[i-1] > 0) {
		for (rle = i+1; rle<256 && F0[rle]; rle++)
		    ;
		rle -= i+1;
		out.WriteByte(rle);
	    }
	}

	WriteFrequencies0(out, F[i]);
    }
    out.WriteByte(0);
}

function RansEncode1(src) {
    const nbytes = src.length;
    var output = new IOStream("", 0, 257*257*3+9);

    output.WriteByte(1);   // Order 0
    output.WriteUint32(0); // compressed size: correct later
    output.WriteUint32(0); // uncompressed size: correct later

    // Compute frequencies
    var F0 = new Array(256)
    var F = new Array(256)
    var C = new Array(256)
    for (var i = 0; i < 256; i++) {
	F[i] = new Array(256);
	C[i] = new Array(256);
    }

    BuildFrequencies1(src, F, F0)
    NormaliseFrequencies1(F, F0);
    WriteFrequencies1(output, F, F0);

    // Compute cumulative frequencies
    for (var i = 0; i < 256; i++) {
	if (!F0[i])
	    continue;

	C[i][0] = 0;
	for (var j = 1; j < 256; j++)
	    C[i][j] = C[i][j-1] + F[i][j-1];
    }

    // Initialise rANS state
    var R = new Array(4);
    var L = new Array(4);
    for (var j = 0; j < 4; j++) {
	R[j] = RansEncInit();
	L[j] = 0;
    }
    var rans_out = new IOStream("", nbytes, nbytes);

    // We have 4 rans codecs running in parallel on its own 1/4tr of buffer
    var nbytes4 = Math.floor(nbytes/4);
    var idx = new Array(4);
    var last = new Array(4)
    for (var j = 0; j < 4; j++) {
	idx[j] = (j+1)*nbytes4 - 2;
	last[j] = src[idx[j]+1]
    }

    // Deal with the remainder if not a multiple of 4
    last[3] = src[nbytes-1];
    for (var i = nbytes-2; i > 4*nbytes4-2; i--) {
	R[3] = RansEncPut(R[3], rans_out, C[src[i]][last[3]], F[src[i]][last[3]], 12);
	last[3] = src[i];
    }

    // Main encode loop
    while (idx[0] >= 0) {
	for (var j = 3; j >= 0; j--) {
	    var s = src[idx[j]]
	    R[j] = RansEncPut(R[j], rans_out, C[s][last[j]], F[s][last[j]], 12);
	    last[j] = s;
	    idx[j]--;
	}
    }

    for (var j = 3; j >= 0; j--) {
        R[j] = RansEncPut(R[j], rans_out, C[0][last[j]], F[0][last[j]], 12)
    }

    for (var i = 3; i >= 0; i--)
	RansEncFlush(R[i], rans_out);

    // Stitch blocks together into final output buffer
    var freq_tab = output.pos;
    output.buf.writeInt32LE(freq_tab-9 + (rans_out.length - rans_out.pos), 1);
    output.buf.writeInt32LE(nbytes, 5);

    return Buffer.concat([output.buf.slice(0, output.pos),
			  rans_out.buf.slice(rans_out.pos, rans_out.length)],
			 output.pos + rans_out.length - rans_out.pos);
}

module.exports = { decode, encode }
