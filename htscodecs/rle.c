#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "varint.h"
#include "rle.h"

#define MAGIC 8

//-----------------------------------------------------------------------------
// Auto compute rle_syms / rle_nsyms
static void rle_find_syms(uint8_t *data, uint64_t data_len,
			  int64_t *saved, // dim >= 256 
			  uint8_t *rle_syms, int *rle_nsyms) {
    int last = -1, n;
    uint64_t i;

    if (data_len > 256) {
	// 186/450
	// Interleaved buffers to avoid cache collisions
	int64_t saved2[256+MAGIC] = {0};
	int64_t saved3[256+MAGIC] = {0};
	int64_t saved4[256+MAGIC] = {0};
	int64_t len4 = data_len&~3;
	for (i = 0; i < len4; i+=4) {
	    int d1 = (data[i+0] == last)     <<1;
	    int d2 = (data[i+1] == data[i+0])<<1;
	    int d3 = (data[i+2] == data[i+1])<<1;
	    int d4 = (data[i+3] == data[i+2])<<1;
	    last = data[i+3];
	    saved [data[i+0]] += d1-1;
	    saved2[data[i+1]] += d2-1;
	    saved3[data[i+2]] += d3-1;
	    saved4[data[i+3]] += d4-1;
	}
	while (i < data_len) {
	    int d = (data[i] == last)<<1;
	    saved[data[i]] += d - 1;
	    last = data[i];
	    i++;
	}
	for (i = 0; i < 256; i++)
	    saved[i] += saved2[i] + saved3[i] + saved4[i];
    } else {
	// 163/391
	for (i = 0; i < data_len; i++) {
	    if (data[i] == last) {
		saved[data[i]]++;
	    } else {
		saved[data[i]]--;
		last = data[i];
	    }
	}
    }

    // Map back to a list
    for (i = n = 0; i < 256; i++) {
	if (saved[i] > 0)
	    rle_syms[n++] = i;
    }
    *rle_nsyms = n;
}

uint8_t *rle_encode(uint8_t *data, uint64_t data_len,
		    uint8_t *run,  uint64_t *run_len,
		    uint8_t *rle_syms, int *rle_nsyms,
		    uint8_t *out, uint64_t *out_len) {
    uint64_t i, j, k;
    if (!out)
	if (!(out = malloc(data_len*2)))
	    return NULL;

    // Two pass:  Firstly compute which symbols are worth using RLE on.
    int64_t saved[256+MAGIC] = {0};

    if (*rle_nsyms) {
	for (i = 0; i < *rle_nsyms; i++)
	    saved[rle_syms[i]] = 1;
    } else {
	// Writes back to rle_syms and rle_nsyms
	rle_find_syms(data, data_len, saved, rle_syms, rle_nsyms);
    }

    // 2nd pass: perform RLE itself to out[] and run[] arrays.
    for (i = j = k = 0; i < data_len; i++) {
	out[k++] = data[i];
	if (saved[data[i]] > 0) {
	    int rlen = i;
	    int last = data[i];
	    while (i < data_len && data[i] == last)
		i++;
	    i--;
	    rlen = i-rlen;

	    j += var_put_u32(&run[j], NULL, rlen);
	}
    }
    
    *run_len = j;
    *out_len = k;
    return out;
}

// On input *out_len holds the allocated size of out[].
// On output it holds the used size of out[].
uint8_t *rle_decode(uint8_t *lit, uint64_t lit_len,
		    uint8_t *run, uint64_t run_len,
		    uint8_t *rle_syms, int rle_nsyms,
		    uint8_t *out, uint64_t *out_len) {
    uint64_t j;
    uint8_t *run_end = run + run_len;

#ifdef FUZZING_BUILD_MODE_UNSAFE_FOR_PRODUCTION
    if (*out_len > 100000)
	return NULL;
#endif

    int saved[256] = {0};
    for (j = 0; j < rle_nsyms; j++)
	saved[rle_syms[j]] = 1;

    j = 0;
    uint8_t *lit_end = lit + lit_len;
    uint8_t *out_end = out + *out_len;
    uint8_t *outp = out;
    while (lit < lit_end) {
	uint32_t rlen;
	uint8_t b = *lit++;
	if (saved[b]) {
	    uint32_t v;
	    if ((v = var_get_u32(run, run_end, &rlen)) == 0)
		return NULL;
	    run += v;
	    rlen++;
	    if (outp + rlen > out_end)
		goto err;
	    memset(outp, b, rlen);
	    outp += rlen;
	} else {
	    if (outp >= out_end)
		goto err;
	    *outp++ = b;
	}
    }

    *out_len = outp-out;
    return out;

 err:
    free(out);
    return NULL;
}

