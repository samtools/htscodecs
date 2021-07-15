#ifndef RANS_INTERNAL_H
#define RANS_INTERNAL_H

#include "config.h"
#include "varint.h"

/*
 * Copyright (c) 2017-2021 Genome Research Ltd.
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

// Internal: common parts to all the rANSNx16pr implementations.

// As per standard rANS_static but using optional RLE or bit-packing
// techniques prior to entropy encoding.  This is a significant
// reduction in some data sets.

// top bits in order byte
#define X_PACK   0x80    // Pack 2,4,8 or infinite symbols into a byte.
#define X_RLE    0x40    // Run length encoding with runs & lits encoded separately
#define X_CAT    0x20    // Nop; for tiny segments where rANS overhead is too big
#define X_NOSZ   0x10    // Don't store the original size; used by STRIPE mode
#define X_STRIPE 0x08    // For N-byte integer data; rotate & encode N streams.
#define X_32     0x04    // 32-way unrolling instead of 4-way

// Not part of the file format, but used to direct the encoder
#define X_SIMD_AUTO 0x100 // automatically enable X_32 if we deem it worthy
#define X_SW32_ENC  0x200 // forcibly use the software version of X_32
#define X_SW32_DEC  0x400 // forcibly use the software version of X_32
#define X_NO_AVX512 0x800 // turn off avx512, but permits AVX2

#define TF_SHIFT 12
#define TOTFREQ (1<<TF_SHIFT)


// 9-11 is considerably faster in the O1 variant due to reduced table size.
// We auto-tune between 10 and 12 though.  Anywhere from 9 to 14 are viable.
#ifndef TF_SHIFT_O1
#define TF_SHIFT_O1 12
#endif
#ifndef TF_SHIFT_O1_FAST
#define TF_SHIFT_O1_FAST 10
#endif
#define TOTFREQ_O1 (1<<TF_SHIFT_O1)
#define TOTFREQ_O1_FAST (1<<TF_SHIFT_O1_FAST)

unsigned char *rans_compress_O0_4x16(unsigned char *in, unsigned int in_size,
				     unsigned char *out, unsigned int *out_size);
unsigned char *rans_uncompress_O0_4x16(unsigned char *in, unsigned int in_size,
				       unsigned char *out, unsigned int out_sz);

int compute_shift(uint32_t *F0, uint32_t (*F)[256], uint32_t *T, int *S);

// Rounds to next power of 2.
// credit to http://graphics.stanford.edu/~seander/bithacks.html
static inline uint32_t round2(uint32_t v) {
    v--;
    v |= v >> 1;
    v |= v >> 2;
    v |= v >> 4;
    v |= v >> 8;
    v |= v >> 16;
    v++;
    return v;
}

static inline int normalise_freq(uint32_t *F, int size, uint32_t tot) {
    int m, M, j, loop = 0;
    uint64_t tr;
    if (!size)
	return 0;

 again:
    tr = ((uint64_t)tot<<31)/size + (1<<30)/size;

    for (size = m = M = j = 0; j < 256; j++) {
	if (!F[j])
	    continue;

	if (m < F[j])
	    m = F[j], M = j;

	if ((F[j] = (F[j]*tr)>>31) == 0)
	    F[j] = 1;
	size += F[j];
//	if (F[j] == tot)
//	    F[j]--;
    }

    int adjust = tot - size;
    if (adjust > 0) {
	F[M] += adjust;
    } else if (adjust < 0) {
	if (F[M] > -adjust && (loop == 1 || F[M]/2 >= -adjust)) {
	    F[M] += adjust;
	} else {
	    if (loop < 1) {
		loop++;
		goto again;
	    }
	    adjust += F[M]-1;
	    F[M] = 1;
	    for (j = 0; adjust && j < 256; j++) {
		if (F[j] < 2) continue;

		int d = F[j] > -adjust;
		int m = d ? adjust : 1-F[j];
		F[j]   += m;
		adjust -= m;
	    }
	}
    }

    //printf("F[%d]=%d\n", M, F[M]);
    return F[M]>0 ? 0 : -1;
}

// A specialised version of normalise_freq_shift where the input size
// is already normalised to a power of 2, meaning we can just perform
// shifts instead of hard to define multiplications and adjustments.
static inline void normalise_freq_shift(uint32_t *F, uint32_t size,
					uint32_t max_tot) {
    if (size == 0 || size == max_tot)
	return;

    int shift = 0, i;
    while (size < max_tot)
	size*=2, shift++;

    for (i = 0; i < 256; i++)
	F[i] <<= shift;
}

// symbols only
static inline int encode_alphabet(uint8_t *cp, uint32_t *F) {
    uint8_t *op = cp;
    int rle, j;

    for (rle = j = 0; j < 256; j++) {
	if (F[j]) {
	    // j
	    if (rle) {
		rle--;
	    } else {
		*cp++ = j;
		if (!rle && j && F[j-1])  {
		    for(rle=j+1; rle<256 && F[rle]; rle++)
			;
		    rle -= j+1;
		    *cp++ = rle;
		}
		//fprintf(stderr, "%d: %d %d\n", j, rle, N[j]);
	    }
	}
    }
    *cp++ = 0;
    
    return cp - op;
}

static inline int decode_alphabet(uint8_t *cp, uint8_t *cp_end, uint32_t *F) {
    if (cp == cp_end)
	return 0;

    uint8_t *op = cp;
    int rle = 0;
    int j = *cp++;
    if (cp+2 >= cp_end)
	goto carefully;

    do {
	F[j] = 1;
	if (!rle && j+1 == *cp) {
	    j = *cp++;
	    rle = *cp++;
	} else if (rle) {
	    rle--;
	    j++;
	    if (j > 255)
		return 0;
	} else {
	    j = *cp++;
	}
    } while(j && cp+2 < cp_end);

 carefully:
    if (j) {
	do {
	    F[j] = 1;
	    if(cp >= cp_end) return 0;
	    if (!rle && j+1 == *cp) {
		if (cp+1 >= cp_end) return 0;
		j = *cp++;
		rle = *cp++;
	    } else if (rle) {
		rle--;
		j++;
		if (j > 255)
		    return 0;
	    } else {
		if (cp >= cp_end) return 0;
		j = *cp++;
	    }
	} while(j && cp < cp_end);
    }

    return cp - op;
}

static inline int encode_freq(uint8_t *cp, uint32_t *F) {
    uint8_t *op = cp;
    int j;

    cp += encode_alphabet(cp, F);

    for (j = 0; j < 256; j++) {
	if (F[j])
	    cp += var_put_u32(cp, NULL, F[j]);
    }

    return cp - op;
}

static inline int decode_freq(uint8_t *cp, uint8_t *cp_end, uint32_t *F,
			      uint32_t *fsum) {
    if (cp == cp_end)
	return 0;

    uint8_t *op = cp;
    cp += decode_alphabet(cp, cp_end, F);

    int j, tot = 0;
    for (j = 0; j < 256; j++) {
	if (F[j]) {
	    cp += var_get_u32(cp, cp_end, (unsigned int *)&F[j]);
	    tot += F[j];
	}
    }

    *fsum = tot;
    return cp - op;
}


// Use the order-0 freqs in F0 to encode the order-1 stats in F.
// All symbols present in F are present in F0, but some in F0 will
// be empty in F.  Thus we run-length encode the 0 frequencies.
static inline int encode_freq_d(uint8_t *cp, uint32_t *F0, uint32_t *F) {
    uint8_t *op = cp;
    int j, dz;

    for (dz = j = 0; j < 256; j++) {
	if (F0[j]) {
	    if (F[j] != 0) {
		if (dz) {
		    // Replace dz zeros with zero + dz-1 run length
		    cp -= dz-1;
		    *cp++ = dz-1;
		}
		dz = 0;
		cp += var_put_u32(cp, NULL, F[j]);
	    } else {
		//fprintf(stderr, "2: j=%d F0[j]=%d, F[j]=%d, dz=%d\n", j, F0[j], F[j], dz);
		dz++;
		*cp++ = 0;
	    }
	} else {
	    assert(F[j] == 0);
	}
    }
    
    if (dz) {
	cp -= dz-1;
	*cp++ = dz-1;
    }

    return cp - op;
}

static inline int decode_freq_d(uint8_t *cp, uint8_t *cp_end, uint32_t *F0,
				uint32_t *F, uint32_t *total) {
    if (cp == cp_end)
	return 0;

    uint8_t *op = cp;
    int j, dz, T = 0;

    for (j = dz = 0; j < 256 && cp < cp_end; j++) {
	//if (F0[j]) fprintf(stderr, "F0[%d]=%d\n", j, F0[j]);
	if (!F0[j])
	    continue;

	uint32_t f;
	if (dz) {
	    f = 0;
	    dz--;
	} else {
	    if (cp >= cp_end) return 0;
	    cp += var_get_u32(cp, cp_end, &f);
	    if (f == 0) {
		if (cp >= cp_end) return 0;
		dz = *cp++;
	    }
	}
	F[j] = f;
	T += f;
    }

    if (total) *total = T;
    return cp - op;
}

#endif // RANS_INTERNAL_H
