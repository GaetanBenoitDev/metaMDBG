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

var fs = require("fs");
var rans = require("./rans4x16");
var argv = require('minimist')(process.argv.slice(2), { boolean: ["d", "r"] });

if (argv._.length != 1) {
    processf.stderr.write("Usage: node main_rans4x16.js [-d] [-o order] input-file > output-file\n");
    process.exit(1);
}

var filein  = argv._[0]

var buf = fs.readFileSync(filein);
var blk_size = 1024*1024;
var raw = argv.r

if (!argv.d) {
    var order = argv.o != undefined ? argv.o : 0;
    // -o8.4 => 8+(256*4)
    order += Math.round((order - (order>>0))*10)*256
    var pos = 0;
    var out_len = 0;
    if (raw)
	blk_size = buf.length
    while (pos < buf.length) {
	var buf2 = rans.encode(buf.slice(pos, pos+blk_size), order);
	var header = new Buffer.allocUnsafe(4);
	if (!raw) {
	    header.writeInt32LE(buf2.length, 0);
	    process.stdout.write(header)
	}
	process.stdout.write(buf2)
	pos += blk_size;
	out_len += buf2.length+4;
    }
    process.stderr.write("Compress order "+order+", "+buf.length+" => " + out_len + "\n");

} else {
    var pos = 0;
    var out_len = 0;
    var len = buf.length
    while (pos < buf.length) {
	if (!raw) {
	    len = buf.readInt32LE(pos);
	    pos += 4;
	}
	var buf2 = rans.decode(buf.slice(pos, pos+len));
	process.stdout.write(buf2)
	out_len += buf2.length;
	pos += len;
    }
    process.stderr.write("Decompress " + buf.length + " => " + out_len + "\n");
}
