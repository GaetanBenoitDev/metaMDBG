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

// An arithmetic coder, based on Eugene Shelwien's reimplementation of
// Michael Schindler range coder.
//
// Order-0 byte stream of ~/scratch/data/q40b
// C:              3.1s decode  (approx same vs 32-bit and 64-bit)
// Arith_sh.js     6.7s decode  (32-bit with carries)
// Arith.js      317.0s decode  (64-bit no carries); int64 crippling it.

//----------------------------------------------------------------------
// Arithmetic (range) coder
module.exports = class RangeCoder {
    constructor(src) {
	this.low   = 0;
	this.range = 0xffffffff;
	this.code  = 0;
	this.FFnum = 0;
	this.carry = 0;
	this.cache = 0;
    }

    RangeStartDecode(src) {
	for (var i = 0; i < 5; i++)
	    this.code = (this.code << 8) + src.ReadByte();
	this.code &= 0xffffffff;
	this.code >>>= 0; // force to be +ve int
    }

    RangeGetFrequency(tot_freq) {
	this.range = Math.floor(this.range / tot_freq);
	//return this.code / this.range;
	return Math.floor(this.code / this.range);

	// Conceptual scenario; return freq only and don't modify range yet
	//return Math.floor(this.code / (Math.floor(this.range / tot_freq)));
    }

    RangeDecode(src, sym_low, sym_freq, tot_freq) {
	// Conceptually we divide range here, but in practice we cached it earlier
	//this.range = Math.floor(this.range / tot_freq);

	this.code  -= sym_low * this.range;
	this.range *= sym_freq;

	while (this.range < (1<<24)) {
	    this.range *= 256;
	    this.code = (this.code*256 + src.ReadByte());
	}
    }

    RangeShiftLow(dst) {
	// We know range is < (1<<24) as we got here.  We already have a
	// cached copy of 8 bits from low.  Is this correct, or does it need
	// fixing?  Possible scenarios.
	// 1. Low < 0xff000000 thus low+range < 0xffffffff and cache
	//    cannot possibly change.  Output cache and as many ffs as needed.
	// 2. We already detected an overflow in RangeEncode, setting carry.
	//    In this case output cached byte + 1 and any 00s needed.
	// 3. Neither case - range is low but we haven't yet detected if we're
	//    XXffffff or XY000000 scenario.  Increase counter for ff/00s.

	if (this.low < 0xff000000 | this.carry) {
	    // cached byte if no overflow, byte+1 otherwise
	    dst.WriteByte(this.cache + this.carry);

	    // Flush any tracked FFs (no carry) or 00s (carry).
	    while (this.FFnum) {
		dst.WriteByte(this.carry-1);
		this.FFnum--;
	    }

	    // Take a copy of top byte ready for next flush
	    this.cache = this.low >>> 24;
	    this.carry = 0;
	} else {
	    this.FFnum++; // keep track of number of trailing ff/00 bytes to write
	}
	this.low <<= 8;
	this.low >>>= 0; // force to be +ve int
    }

    RangeEncode(dst, sym_low, sym_freq, tot_freq) {
	var old_low = this.low
	this.range  = Math.floor(this.range / tot_freq)
	this.low   += sym_low * this.range;
	this.low >>>= 0; // Truncate to +ve int so we can spot overflow
	this.range *= sym_freq;

	// "low + sym*range < old_low" means we overflow; set carry.
	// NB: can this.low < old_low occur twice before range < (1<<24)?
	// We claim not, but prove it!
	if (this.low < old_low) {
	    if (this.carry != 0) console.log("ERROR: Multiple carry")
	    this.carry = 1
	}

	// Renormalise if range gets too small
	while (this.range < (1<<24)) {
	    this.range *= 256;
	    this.RangeShiftLow(dst);
	}
    }

    RangeFinishEncode(dst) {
	for (var i = 0; i < 5; i++)
	    this.RangeShiftLow(dst)
    }
};
