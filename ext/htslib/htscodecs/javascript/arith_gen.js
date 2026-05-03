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

const RangeCoder = require("./arith_sh");
const IOStream = require("./iostream");
const ByteModel = require("./byte_model");
const bzip2 = require("bzip2");

const ARITH_ORDER  = 1
const ARITH_EXT    = 4
const ARITH_STRIPE = 8
const ARITH_NOSIZE = 16
const ARITH_CAT    = 32
const ARITH_RLE    = 64
const ARITH_PACK   = 128

module.exports = class RangeCoderGen {
    decode(src) {
	this.stream = new IOStream(src);
	return this.decodeStream(this.stream)
    }

    decodeStream(stream, n_out=0) {
	var flags = this.stream.ReadByte();
	if (!(flags & ARITH_NOSIZE))
	    n_out = this.stream.ReadUint7();
	var e_len = n_out;

	var order = flags & ARITH_ORDER;

	// 4-way recursion
	if (flags & ARITH_STRIPE)
	    return this.decodeStripe(this.stream, n_out)

	// Meta data
	if (flags & ARITH_PACK) {
	    var P
	    [P, e_len] = this.decodePackMeta(this.stream)
	}

	// NOP, useful for tiny blocks
	if (flags & ARITH_CAT)
	    var data = this.decodeCat(this.stream, e_len)

	// Entropy decode
	else if (flags & ARITH_EXT) {
	    var data = this.decodeExt(this.stream, e_len)
	} else if (flags & ARITH_RLE) {
	    var data = order
		? this.decodeRLE1(this.stream, e_len)
		: this.decodeRLE0(this.stream, e_len)
	} else {
	    var data = order
		? this.decode1(this.stream, e_len)
		: this.decode0(this.stream, e_len)
	}

	// Transforms
	if (flags & ARITH_PACK)
	    data = this.decodePack(data, P, n_out)

	return data
    }

    encode(src, flags) {
	this.stream = new IOStream("", 0, src.length*1.1 + 100); // guestimate worst case!

	this.stream.WriteByte(flags);
	if (!(flags & ARITH_NOSIZE))
	    this.stream.WriteUint7(src.length);

	if (flags & ARITH_STRIPE)
	    return Buffer.concat([this.stream.buf.slice(0, this.stream.pos),
				  this.encodeStripe(this.stream, src, flags>>8)])

	var order = flags & ARITH_ORDER;
	var e_len = src.length;

	// step 1: Encode meta-data
	var pack_meta
	if (flags & ARITH_PACK)
	    [pack_meta, src, e_len] = this.encodePack(src)

	// step 2: Write any meta data
	if (flags & ARITH_PACK)
	    this.stream.WriteStream(pack_meta)

	// step 3: arith encoding below
	if (flags & ARITH_RLE) {
	    return order
		? this.encodeRLE1(src, e_len, this.stream)
		: this.encodeRLE0(src, e_len, this.stream);
	} else {
	    return order
		? this.encode1(src, e_len, this.stream)
		: this.encode0(src, e_len, this.stream);
	}
    }

    //----------------------------------------------------------------------
    // Order-0 codec
    decode0(stream, n_out) {
	var output = new Buffer.allocUnsafe(n_out);

	var max_sym = stream.ReadByte()
	if (max_sym == 0)
	    max_sym = 256

	var byte_model = new ByteModel(max_sym);

	var rc = new RangeCoder(stream);
	rc.RangeStartDecode(stream);

	for (var i = 0; i < n_out; i++)
	    output[i] = byte_model.ModelDecode(stream, rc);

	return output;
    }

    encode0(src, n_in, out) {
	// Count the maximum symbol present
	var max_sym = 0;
	for (var i = 0; i < n_in; i++)
	    if (max_sym < src[i])
		max_sym = src[i]
	max_sym++;  // FIXME not what spec states!?

	var byte_model = new ByteModel(max_sym);
	out.WriteByte(max_sym);
	var rc = new RangeCoder(out);

	for (var i = 0; i < n_in; i++)
	    byte_model.ModelEncode(out, rc, src[i])
	rc.RangeFinishEncode(out)

	return out.buf.slice(0, out.pos);
    }

    //----------------------------------------------------------------------
    // Order-1 codec

    decode1(stream, n_out) {
	var output = new Buffer.allocUnsafe(n_out);

	var max_sym = stream.ReadByte()
	if (max_sym == 0)
	    max_sym = 256

	var byte_model = new Array(max_sym);
	for (var i = 0; i < max_sym; i++)
	    byte_model[i] = new ByteModel(max_sym);

	var rc = new RangeCoder(stream);
	rc.RangeStartDecode(stream);

	var last = 0;
	for (var i = 0; i < n_out; i++) {
	    output[i] = byte_model[last].ModelDecode(stream, rc);
	    last = output[i];
	}

	return output;
    }

    encode1(src, n_in, out) {
	// Count the maximum symbol present
	var max_sym = 0;
	for (var i = 0; i < n_in; i++)
	    if (max_sym < src[i])
		max_sym = src[i]
	max_sym++;  // FIXME not what spec states!

	var byte_model = new Array(max_sym);
	for (var i = 0; i < max_sym; i++)
	    byte_model[i] = new ByteModel(max_sym);
	out.WriteByte(max_sym);
	var rc = new RangeCoder(out);

	var last = 0;
	for (var i = 0; i < n_in; i++) {
	    byte_model[last].ModelEncode(out, rc, src[i])
	    last = src[i]
	}
	rc.RangeFinishEncode(out)

	return out.buf.slice(0, out.pos);
    }

    //----------------------------------------------------------------------
    // External codec
    decodeExt(stream, n_out) {
	// Bzip2 only for now
	var output = new Buffer.allocUnsafe(n_out)
	var bits = bzip2.array(stream.buf.slice(stream.pos))
	var size = bzip2.header(bits)
	var j = 0
	do {
	    var chunk = bzip2.decompress(bits, size);
	    if (chunk != -1) {
	        Buffer.from(chunk).copy(output, j)
	        j += chunk.length
		size -= chunk.length
	    }
	} while(chunk != -1);

	return output
    }

    encodeExt(stream, n_out) {
	// We cannot compress using Bzip2 now as it's
	// absent from bzip2.js, but consider using
	// https://github.com/cscott/compressjs
    }

    //----------------------------------------------------------------------
    // Order-0 RLE codec
    decodeRLE0(stream, n_out) {
	var output = new Buffer.allocUnsafe(n_out);

	var max_sym = stream.ReadByte()
	if (max_sym == 0)
	    max_sym = 256

	var model_lit = new ByteModel(max_sym);
	var model_run = new Array(258);
	for (var i = 0; i <= 257; i++)
	    model_run[i] = new ByteModel(4)

	var rc = new RangeCoder(stream);
	rc.RangeStartDecode(stream);

	var i = 0;
	while (i < n_out) {
	    output[i] = model_lit.ModelDecode(stream, rc)
	    var part = model_run[output[i]].ModelDecode(stream, rc)
	    var run = part
	    var rctx = 256
	    while (part == 3) {
		part = model_run[rctx].ModelDecode(stream, rc)
		rctx = 257
		run += part
	    }
	    for (var j = 1; j <= run; j++)
		output[i+j] = output[i]
	    i += run+1
	}

	return output;
    }

    encodeRLE0(src, n_in, out) {
	// Count the maximum symbol present
	var max_sym = 0;
	for (var i = 0; i < n_in; i++)
	    if (max_sym < src[i])
		max_sym = src[i]
	max_sym++;  // FIXME not what spec states!

	var model_lit = new ByteModel(max_sym);
	var model_run = new Array(258);
	for (var i = 0; i <= 257; i++)
	    model_run[i] = new ByteModel(4)

	out.WriteByte(max_sym);
	var rc = new RangeCoder(out);

	var i = 0
	while (i < n_in) {
	    model_lit.ModelEncode(out, rc, src[i])
	    var run = 1
	    while (i+run < n_in && src[i+run] == src[i])
		run++
	    run--

	    var rctx = src[i]
	    var last = src[i]
	    i += run+1

	    var part = run >= 3 ? 3 : run
	    model_run[rctx].ModelEncode(out, rc, part)
	    run -= part
	    rctx = 256
	    while (part == 3) {
		part = run >= 3 ? 3 : run
		model_run[rctx].ModelEncode(out, rc, part)
		rctx = 257
		run -= part
	    }
	}
	rc.RangeFinishEncode(out)

	return out.buf.slice(0, out.pos);
    }

    //----------------------------------------------------------------------
    // Order-1 RLE codec

    decodeRLE1(stream, n_out) {
	var output = new Buffer.allocUnsafe(n_out);

	var max_sym = stream.ReadByte()
	if (max_sym == 0)
	    max_sym = 256

	var model_lit = new Array(max_sym);
	for (var i = 0; i < max_sym; i++)
	    model_lit[i] = new ByteModel(max_sym);

	var model_run = new Array(258);
	for (var i = 0; i <= 257; i++)
	    model_run[i] = new ByteModel(4)

	var rc = new RangeCoder(stream);
	rc.RangeStartDecode(stream);

	var last = 0;
	var i = 0;
	while (i < n_out) {
	    output[i] = model_lit[last].ModelDecode(stream, rc)
	    last = output[i]
	    var part = model_run[output[i]].ModelDecode(stream, rc)
	    var run = part
	    var rctx = 256
	    while (part == 3) {
		part = model_run[rctx].ModelDecode(stream, rc)
		rctx = 257
		run += part
	    }
	    for (var j = 1; j <= run; j++)
		output[i+j] = output[i]
	    i += run+1
	}

	return output;
    }

    encodeRLE1(src, n_in, out) {
	// Count the maximum symbol present
	var max_sym = 0;
	for (var i = 0; i < n_in; i++)
	    if (max_sym < src[i])
		max_sym = src[i]
	max_sym++;  // FIXME not what spec states!

	var model_lit = new Array(max_sym)
	for (var i = 0; i < max_sym; i++)
	    model_lit[i] = new ByteModel(max_sym);
	var model_run = new Array(258);
	for (var i = 0; i <= 257; i++)
	    model_run[i] = new ByteModel(4)

	out.WriteByte(max_sym);
	var rc = new RangeCoder(out);

	var i = 0
	var last = 0
	while (i < n_in) {
	    model_lit[last].ModelEncode(out, rc, src[i])
	    var run = 1
	    while (i+run < n_in && src[i+run] == src[i])
		run++
	    run--

	    var rctx = src[i]
	    last = src[i]
	    i += run+1

	    var part = run >= 3 ? 3 : run
	    model_run[rctx].ModelEncode(out, rc, part)
	    run -= part
	    rctx = 256
	    while (part == 3) {
		part = run >= 3 ? 3 : run
		model_run[rctx].ModelEncode(out, rc, part)
		rctx = 257
		run -= part
	    }
	}
	rc.RangeFinishEncode(out)

	return out.buf.slice(0, out.pos);
    }

    //----------------------------------------------------------------------
    // Pack method
    decodePackMeta(stream) {
	this.nsym  = stream.ReadByte()

	var M = new Array(this.nsym);
	for (var i = 0; i < this.nsym; i++)
	    M[i] = stream.ReadByte()

	var e_len = stream.ReadUint7(); // Could be derived data from nsym and n_out

	return [M, e_len]
    }

    decodePack(data, M, len) {
	var out = new Buffer.allocUnsafe(len);

	if (this.nsym <= 1) {
	    // Constant value
	    for (var i = 0; i < len; i++)
		out[i] = M[0]

	} else if (this.nsym <= 2) {
	    // 1 bit per value
	    for (var i = 0, j = 0; i < len; i++) {
		if (i % 8 == 0)
		    var v = data[j++]
		out[i] = M[v & 1]
		v >>= 1
	    }

	} else if (this.nsym <= 4) {
	    // 2 bits per value
	    for (var i = 0, j = 0; i < len; i++) {
		if (i % 4 == 0)
		    var v = data[j++]
		out[i] = M[v & 3]
		v >>= 2
	    }

	} else if (this.nsym <= 16) {
	    // 4 bits per value
	    for (var i = 0, j = 0; i < len; i++) {
		if (i % 2 == 0)
		    var v = data[j++]
		out[i] = M[v & 15]
		v >>= 4
	    }

	} else {
	    // 8 bits per value: NOP
	    return data
	}

	return out
    }

    // Compute M array and return meta-data stream
    packMeta(src) {
	var stream = new IOStream("", 0, 1024)

	// Count symbols
	var M = new Array(256)
	for (var i = 0; i < src.length; i++)
	    M[src[i]] = 1

	// Write Map
	for (var nsym = 0, i = 0; i < 256; i++)
	    if (M[i])
		M[i] = ++nsym; // map to 1..N
	stream.WriteByte(nsym);

	// FIXME: add check for nsym > 16?
	// Or just accept it as an inefficient waste of time.
	for (var i = 0; i < 256; i++) {
	    if (M[i]) {
		stream.WriteByte(i) // adjust to 0..N-1
		M[i]--;
	    }
	}

	return [stream, M, nsym]
    }

    encodePack(data) {
	var meta, M, nsym
	[meta, M, nsym] = this.packMeta(data)

	var len = data.length
	var i = 0;
	if (nsym <= 1) {
	    // Constant values
	    meta.WriteUint7(0)
	    return [meta, new Buffer.allocUnsafe(0), 0];
	}

	if (nsym <= 2) {
	    // 1 bit per value
	    var out = new Buffer.allocUnsafe(Math.floor((len+7)/8));
	    for (var i = 0, j = 0; i < (len & ~7); i+=8, j++)
		out[j] = (M[data[i+0]]<<0)
		       + (M[data[i+1]]<<1)
		       + (M[data[i+2]]<<2)
		       + (M[data[i+3]]<<3)
		       + (M[data[i+4]]<<4)
		       + (M[data[i+5]]<<5)
		       + (M[data[i+6]]<<6)
		       + (M[data[i+7]]<<7)
	    if (i < len) {
		out[j] = 0;
		var v = 0;
		while (i < len) {
		    out[j] |= M[data[i++]]<<v;
		    v++;
		}
		j++;
	    }

	    meta.WriteUint7(j)
	    return [meta, out, out.length]
	}

	if (nsym <= 4) {
	    // 2 bits per value
	    var out = new Buffer.allocUnsafe(Math.floor((len+3)/4));
	    for (var i = 0, j = 0; i < (len & ~3); i+=4, j++)
		out[j] = (M[data[i+0]]<<0)
		       + (M[data[i+1]]<<2)
		       + (M[data[i+2]]<<4)
		       + (M[data[i+3]]<<6)

	    if (i < len) {
		out[j] = 0;
		var v = 0;
		while (i < len) {
		    out[j] |= M[data[i++]]<<v;
		    v+=2;
		}
		j++;
	    }

	    meta.WriteUint7(j)
	    return [meta, out, out.length]
	}

	if (nsym <= 16) {
	    // 4 bits per value
	    var out = new Buffer.allocUnsafe(Math.floor((len+1)/2));
	    for (var i = 0, j = 0; i < (len & ~1); i+=2, j++)
		out[j] = (M[data[i+0]]<<0)
		       + (M[data[i+1]]<<4)
	    if (i < len)
		out[j++] = M[data[i++]];

	    meta.WriteUint7(j)
	    return [meta, out, out.length]
	}

	// Otherwise an expensive NOP
	meta.WriteUint7(data.length)
	return [meta, data, data.length]
    }

    //----------------------------------------------------------------------
    // STRIPE method
    encodeStripe(hdr, src, N) {
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
	    var comp0 = this.encode(part[s], 0)
	    var comp1 = this.encode(part[s], 1)
	    comp[s] = (comp1.length < comp0.length) ? comp1 : comp0
	    total += comp[s].length
	}

	// Serialise
	var out = new IOStream("", 0, total+5*N + 1)
	out.WriteByte(N)
	for (var s = 0; s < N; s++)
	    out.WriteUint7(comp[s].length)

	for (var s = 0; s < N; s++)
	    out.WriteData(comp[s], comp[s].length)

	return out.buf.slice(0, out.buf.pos)
    }

    decodeStripe(stream, len) {
	var N = stream.ReadByte()
	
	// Retrieve lengths
	var clen = new Array(N)
	var ulen = new Array(N)
	for (var j = 0; j < N; j++)
	    clen[j] = stream.ReadUint7()

	// Decode streams
	var T = new Array(N);
	for (var j = 0; j < N; j++) {
	    ulen[j] = Math.floor(len / N) + ((len % N) > j)
	    T[j] = this.decodeStream(stream, ulen[j])
	}

	// Transpose
	var out = new Buffer.allocUnsafe(len)
	for (var j = 0; j < N; j++) {
	    for (var i = 0; i < ulen[j]; i++) {
		out[i*N + j] = T[j][i];
	    }
	}

	return out
    }

    //----------------------------------------------------------------------
    // Cat method
    decodeCat(stream, len) {
	var out = new Buffer.allocUnsafe(len);
	for (var i = 0; i < len; i++)
	    out[i] = stream.ReadByte()

	return out
    }
}
