/*
 * Copyright (c) 2019,2021 Genome Research Ltd.
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

/*
 * Data transpose by N.  Common to rANS4x16 and arith_dynamic decoders.
 *
 * Tuned for specific common cases of N.
 */
static inline void unstripe(unsigned char *out, unsigned char *outN,
			    unsigned int ulen, unsigned int N,
			    unsigned int idxN[256]) {
    int j = 0, k;

    if (ulen >= N) {
	switch (N) {
	case 4:
	    while (j < ulen-4) {
		for (k = 0; k < 4; k++)
		    out[j++] = outN[idxN[k]++];
	    }
	    break;

	case 2:
	    while (j < ulen-2) {
		for (k = 0; k < 2; k++)
		    out[j++] = outN[idxN[k]++];
	    }
	    break;

	default:
	    // General case, around 25% slower overall decode
	    while (j < ulen-N) {
		for (k = 0; k < N; k++)
		    out[j++] = outN[idxN[k]++];
	    }
	    break;
	}
    }
    for (k = 0; j < ulen; k++)
	out[j++] = outN[idxN[k]++];
}
