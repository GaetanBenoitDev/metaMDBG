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

var fs = require("fs");
var fqz = require("./fqzcomp");
var argv = require('minimist')(process.argv.slice(2), { boolean: ["d", "r"] });

if (argv._.length != 1) {
    process.stderr.write("Usage: node main_fqzcomp.js [-d] input-file > output-file\n");
    process.exit(1);
}

var filein  = argv._[0]

var buf = fs.readFileSync(filein);
var raw = argv.r

if (!argv.d) {
    // Line breaks to get sequence length, but then stitch together into
    // a single non-breaking buffer.
    var len = 0;
    var j = 0;
    var q_lens = new Array
    var q_dirs = new Array
    var q_len = 0
    for (var i = 0; i < buf.length; i++) {
	if (buf[i] == "\n".charCodeAt(0) || buf[i] == "\t".charCodeAt(0)) {
	    q_lens.push(len)
	    if (q_len == 0)
		q_len = len
	    else if (q_len != len)
		q_len = -1 // marker for multiple lengths
	    len = 0;

	    if (buf[i] == "\t".charCodeAt(0)) {
		// parse 2nd token for read1/read2 status
		var dir = ""
		for (i++; i < buf.length && buf[i] != "\n".charCodeAt(0); i++)
		    dir += String.fromCharCode(buf[i])
		q_dirs.push(dir)
	    } else {
		q_dirs.push(0)
	    }
	} else {
	    buf[j++] = buf[i] - 33; // ASCII phred to raw
	    len++;
	}
    }
    buf = buf.slice(0, j)
    if (q_len > 0)
	q_lens = [q_lens[0]]

    var buf2 = fqz.encode(buf, q_lens, q_dirs);
    process.stderr.write("Compress " +buf.length + " => " + buf2.length + "\n");
    if (!raw) {
	var hdr = new Buffer.allocUnsafe(8);
	hdr.writeInt32LE(buf.length, 0);
	hdr.writeInt32LE(buf2.length, 4);
	process.stdout.write(hdr);
    }
    process.stdout.write(buf2);

} else {
    var q_lens = new Array
    // Consume ulen and clen from outer test harness (pointless as non-blocking atm)
    var buf2
    if (raw)
	buf2 = fqz.decode(buf, q_lens);
    else
	buf2 = fqz.decode(buf.slice(8), q_lens);

    // Split into newlines so we can do easy data comparison
    var buf3 = new Buffer.allocUnsafe(buf2.length + q_lens.length)
    var rec = 0;
    var len = q_lens[rec++]
    var j = 0;
    for (var i = 0; i < buf2.length; i++) {
	buf3[j++] = buf2[i] + 33;
	if (--len == 0) {
	    buf3[j++] = "\n".charCodeAt(0)
	    len = q_lens[rec++]
	}
    }

    process.stderr.write("Decompress " + buf.length + " => " + buf3.length + "\n");
    process.stdout.write(buf3);
}
