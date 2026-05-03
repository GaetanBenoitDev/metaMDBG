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

// Name tokeniser
//
// This is a reference implementation designed to match the
// written specification as closely as possible.  It is *NOT*
// an efficient implementation, but see comments below.

const IOStream  = require("./iostream");
const rans      = require("./rans4x16");
const arith_gen = require("./arith_gen");

var arith = new arith_gen()

const TOK_TYPE    = 0
const TOK_STRING  = 1
const TOK_CHAR    = 2
const TOK_DIGITS0 = 3
const TOK_DZLEN   = 4
const TOK_DUP     = 5
const TOK_DIFF    = 6
const TOK_DIGITS  = 7
const TOK_DELTA   = 8
const TOK_DELTA0  = 9
const TOK_MATCH   = 10
const TOK_NOP     = 11
const TOK_END     = 12

//----------------------------------------------------------------------
// Token byte streams
function DecodeTokenByteStreams(src, in_size, use_arith, nnames) {
    var t = -1

    var B = new Array(256)

    while (!src.EOF()) {
	var ttype = src.ReadByte()
	var tok_new = ttype & 128
	var tok_dup = ttype & 64
	var type    = ttype & 63

	if (tok_new) {
	    t++
	    B[t] = new Array(13)
	}

	if (type != TOK_TYPE && tok_new) {
	    var M = new Array(nnames-1).fill(TOK_MATCH)
	    B[t][TOK_TYPE] = new IOStream(Buffer.from([type].concat(M)))
        }

	if (tok_dup) {
	    var dup_pos  = src.ReadByte()
	    var dup_type = src.ReadByte()
	    B[t][type] = new IOStream(B[dup_pos][dup_type].buf)
	} else {
	    var clen = src.ReadUint7()
	    var data = src.ReadData(clen)

	    if (use_arith)
		B[t][type] = arith.decode(data)
	    else
		B[t][type] = rans.decode(data)
	    B[t][type] = new IOStream(B[t][type])
	}
    }

    return B
}

//----------------------------------------------------------------------
// Token decode
function LeftPadNumber(val, len) {
    var str = val+""
    while (str.length < len)
	str = "0" + str

    return str
}

function DecodeSingleName(B, N, T, n) {
    var type = B[0][TOK_TYPE].ReadByte()
    var dist = B[0][type].ReadUint32()
    var m = n - dist

    if (type == TOK_DUP) {
	N[n] = N[m]
	T[n] = T[m]
	return N[n]
    }
    
    var t = 1
    N[n] = ""
    T[n] = new Array(256)
    do {
	type = B[t][TOK_TYPE].ReadByte()

	switch(type) {
	case TOK_CHAR:
	    T[n][t] = B[t][TOK_CHAR].ReadChar()
	    break

	case TOK_STRING:
	    T[n][t] = B[t][TOK_STRING].ReadString()
	    break
	
	case TOK_DIGITS:
	    T[n][t] = B[t][TOK_DIGITS].ReadUint32()
	    break

	case TOK_DIGITS0:
	    var d = B[t][TOK_DIGITS0].ReadUint32()
	    var l = B[t][TOK_DZLEN].ReadByte()
	    T[n][t] = LeftPadNumber(d, l)
	    break

	case TOK_DELTA:
	    T[n][t] = (T[m][t]>>0) + B[t][TOK_DELTA].ReadByte()
	    break

	case TOK_DELTA0:
	    var d = (T[m][t]>>0) + B[t][TOK_DELTA0].ReadByte()
	    var l = T[m][t].length
	    T[n][t] = LeftPadNumber(d, l)
	    break

	case TOK_MATCH:
	    T[n][t] = T[m][t]
	    break

	default:
	    T[n][t] = ""
	    break
	}

	N[n] += T[n][t++]
    } while (type != TOK_END)

    return N[n]
}

//----------------------------------------------------------------------
// Main tokeniser decode entry function: decodes a compressed src and
// returns the uncompressed buffer.
function decode(src, len, separator) {
    var src = new IOStream(src)
    var ulen = src.ReadUint32()
    var nnames = src.ReadUint32()
    var use_arith = src.ReadByte()

    var B = DecodeTokenByteStreams(src, len, use_arith, nnames)
    var N = new Array(nnames)
    var T = new Array(nnames)

    var str = ""
    if (typeof separator === 'undefined')
	separator = '\n'
    for (var i = 0; i < nnames; i++)
	str += DecodeSingleName(B, N, T, i) + separator

    return str
}

//----------------------------------------------------------------------
// Main tokeniser encode function

// Encoder is trickier than decode as we have a lot of decisions to make.
// However here we just make a simple guess without anything complex,
// to demonstrate the basic idea.  See the C implementation for further
// expansion on this.
function encode(src, use_arith) {
    // Convert buffer to array of names
    var str = src.toString()
    if (str[str.length-1] == '\n')
	str = str.substring(0,str.length-1)
    var names = str.split("\n")

    var out = new IOStream("", 0, str.length*2 + 10000) // guess max size
    out.WriteUint32(str.length)
    out.WriteUint32(names.length)
    out.WriteByte(use_arith)

    // Tokenise names
    var T = new Array(names.length)
    var H = {}
    var F = new Array(256).fill(0) // DELTA vs DIGIT frequency
    var max_tok = 0
    var max_len = 0
    for (var i = 0; i < names.length; i++) {
	var [ntok,len] = TokeniseName(T, H, F, names[i], i)
	if (max_tok < ntok)
	    max_tok = ntok
	if (max_len < len)
	    max_len = len
    }

    // Convert tokens to byte streams and serialise
    for (var tnum = 0; tnum < max_tok; tnum++) {
	var B = new Array(TOK_END+1)
	for (var type = 0; type <= TOK_END; type++)
	    B[type] = new IOStream("", 0, names.length * max_len)

	FillByteStreams(B, T, tnum, names, max_tok, max_len)
	SerialiseByteStreams(B, tnum, use_arith, out)
    }

    return out.buf.slice(0, out.pos)
}

function FillByteStreams(B, T, tnum, names, max_tok, max_len) {
    // Create byte streams B[]
    for (var n = 0; n < names.length; n++) {
	if (tnum > 0 && T[n][0].type == TOK_DUP)
	    continue

	if (!T[n][tnum])
	    continue

	B[TOK_TYPE].WriteByte(T[n][tnum].type)

	switch (T[n][tnum].type) {
	case TOK_DIFF:
	    B[TOK_DIFF].WriteUint32(T[n][tnum].val)
	    break

	case TOK_DUP:
	    B[TOK_DUP].WriteUint32(T[n][tnum].val)
	    break

	case TOK_STRING:
	    B[TOK_STRING].WriteString(T[n][tnum].val)
	    break

	case TOK_CHAR:
	    B[TOK_CHAR].WriteChar(T[n][tnum].val)
	    break

	case TOK_DIGITS:
	    B[TOK_DIGITS].WriteUint32(T[n][tnum].val)
	    break

	case TOK_DIGITS0:
	    B[TOK_DIGITS0].WriteUint32(T[n][tnum].val)
	    B[TOK_DZLEN].WriteByte(T[n][tnum].val.length)
	    break

	case TOK_DELTA:
	    B[T[n][tnum].type].WriteByte(T[n][tnum].val)
	    break

	case TOK_DELTA0:
	    B[T[n][tnum].type].WriteByte(T[n][tnum].val)
	    break
	}
    }
}

function SerialiseByteStreams(B, tnum, use_arith, out) {
    // Compress and serialise byte streams B[]
    for (var type = 0; type <= TOK_END; type++) {
	if (B[type].pos <= 0)
	    continue

	out.WriteByte(type + ((type == 0) ? 128 : 0))

	// IOStream to sized buffer
	B[type] = B[type].buf.slice(0, B[type].pos)
	var comp = try_compress(B[type], use_arith)

	out.WriteUint7(comp.length)
	out.WriteData(comp, comp.length)
    }
}

function try_compress(src, use_arith) {
    var best = 1<<30
    var comp

    var methods = [0, 1, 64, 65, 128, 129, 193+8]
    for (var i in methods) {
	var lvl = methods[i]
	if ((lvl & 1) && src.length < 100)
	    continue

	if ((lvl & 8) && (src.length % 4) != 0)
	    continue

	try {
	    var tmp = use_arith
		? arith.encode(src, lvl)
		: rans.encode(src, lvl)
	} catch (e) {
	    var tmp = 0
	}
	if (tmp && best > tmp.length) {
	    best = tmp.length
	    comp = tmp
	}
    }

    return comp
}

function TokeniseName(T, H, F, name, n) {
    var max_len = 0

    // Always compare against last name only
    var p = n-1
    T[n] = new Array(256)

    if (H[name]) {
	//console.error(name,H[name],n)
	T[n][0] = {
	    type: TOK_DUP,
	    val:  n - H[name]
	}
    } else {
	T[n][0] = {
	    type: TOK_DIFF,
	    val:  n == 0 ? 0 : 1
	}
    }

    H[name] = n

    // Splits on alphanumerics, punctuation
    var tok = name.match(/([a-zA-Z0-9]{1,9})|([^a-zA-Z0-9]+)/g)
    for (var i = 0; i < tok.length; i++) {
	var t = i+1 // token 0 = DIFF vs DUP
	var type = TOK_STRING
	var val = tok[i]
	if (tok[i].match(/^0+[0-9]*$/g))
	    type = TOK_DIGITS0
	else if (tok[i].match(/^[0-9]+$/g))
	    type = TOK_DIGITS
	else if (tok[i].length == 1)
	    type = TOK_CHAR

	if (p >= 0 && T[p][t]) {
	    if (T[p][t].str == tok[i]) {
		type = TOK_MATCH
		val = ""
	    } else if (T[p][t].type == TOK_DIGITS || T[p][t].type == TOK_DELTA) {
		var d = val - T[p][t].str;
		F[t]++
		if (d >= 0 && d < 256 && F[t] > n/2) {
		    type = TOK_DELTA
		    val = d
		}
	    } else if ((T[p][t].type == TOK_DIGITS0 || T[p][t].type == TOK_DELTA0)
		       && T[p][t].str.length == val.length) {
		var d = val - T[p][t].str;
		F[t]++
		if (d >= 0 && d < 256 && F[t] > n/2) {
		    type = TOK_DELTA0
		    val = d
		}
	    }
	}

	T[n][t] = {
	    str:  tok[i],
	    val:  val,
	    type: type
	}

	if (max_len < T[n][t].val.length+3)  // +3 for integers; 5 -> (Uint32)5
	    max_len = T[n][t].val.length+3

	//console.error(t,T[n][t])
    }
    T[n][++t] = {
	type: TOK_END
    }

    return [t+1, max_len]
}

module.exports = { encode, decode }
