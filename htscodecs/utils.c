/*
 * Copyright (c) 2023 Genome Research Ltd.
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

#include "config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <string.h>

#include "utils.h"
#include "rANS_static4x16.h"

#ifndef NO_THREADS
#include <pthread.h>
#endif

//#define TLS_DEBUG

#ifndef NO_THREADS
/*
 * Thread local storage per thread in the pool.
 *
 * We have some large memory blocks for rANS which we cannot store on the
 * stack due to various system limitations.  Allocaitng them can be
 * expensive as some OSes use mmap and will pass the pages back to the OS
 * on each free.  This unfortunately then means zeroing the pages out again
 * on each new malloc, plus additional switching into the kernel.
 *
 * Instead where available, we use pthread_once to allocate a small arena
 * of memory buffers and we continually reuse these same buffers.  We don't
 * need to memset it (calloc equivalent) either as we're sure that any
 * leakage of data is simply an earlier set of precomputed frequency
 * lookups, and not something more sinister such as an encryption key.
 *
 * If we don't have pthreads, then we have to fall back to the slow
 * traditional calloc instead.
 */

#define MAX_TLS_BUFS 10
typedef struct {
    void   *bufs[MAX_TLS_BUFS];
    size_t sizes[MAX_TLS_BUFS];
    int     used[MAX_TLS_BUFS];
} tls_pool;

static pthread_once_t rans_once = PTHREAD_ONCE_INIT;
static pthread_key_t rans_key;

/*
 * Frees all local storage for this thread.
 * Note: this isn't a function to free a specific allocated item.
 */
static void htscodecs_tls_free_all(void *ptr) {
    tls_pool *tls = (tls_pool *)ptr;
    if (!tls)
        return;

    int i;
    for (i = 0; i < MAX_TLS_BUFS; i++) {
#ifdef TLS_DEBUG
        if (tls->bufs[i])
            fprintf(stderr, "Free %ld = %p\n", tls->sizes[i], tls->bufs[i]);
#endif
        if (tls->used[i]) {
            fprintf(stderr, "Closing thread while TLS data is in use\n");
        }
        free(tls->bufs[i]);
    }

    free(tls);
}

static void htscodecs_tls_init(void) {
    pthread_key_create(&rans_key, htscodecs_tls_free_all);
}

/*
 * Allocates size bytes from the global Thread Local Storage pool.
 * This is shared by all subsequent calls within this thread.
 *
 * An simpler alternative could be possible where we have a fixed number
 * of types of alloc, say 5, and specify the correct label when allocating.
 * Eg histogram, name_context, fqzcomp, rans.  We can have multiple types
 * in use in different stack frames (such name_context + hist + rans), but
 * the number is very limited.  That then paves the way to simply check and
 * realloc without needing to keep track of use status or overflowing
 * the maximum number permitted.
 */
void *htscodecs_tls_alloc(size_t size) {
    int i;

    int err = pthread_once(&rans_once, htscodecs_tls_init);
    if (err != 0) {
        fprintf(stderr, "Initialising TLS data failed: pthread_once: %s\n",
                strerror(err));
        return NULL;
    }

    // Initialise tls_pool on first usage
    tls_pool *tls = pthread_getspecific(rans_key);
    if (!tls) {
        if (!(tls = calloc(1, sizeof(*tls))))
            return NULL;
        pthread_setspecific(rans_key, tls);
    }

    // Query pool for size
    int avail = -1;
    for (i = 0; i < MAX_TLS_BUFS; i++) {
        if (!tls->used[i]) {
            if (size <= tls->sizes[i]) {
                tls->used[i] = 1;
#ifdef TLS_DEBUG
                fprintf(stderr, "Reuse %d: %ld/%ld = %p\n",
                        i, size, tls->sizes[i], tls->bufs[i]);
#endif
                return tls->bufs[i];
            } else if (avail == -1) {
                avail = i;
            }
        }
    }

    if (i == MAX_TLS_BUFS && avail == -1) {
        // Shouldn't happen given our very limited use of this function
        fprintf(stderr, "Error: out of rans_tls_alloc slots\n");
        return NULL;
    }

    if (tls->bufs[avail])
        free(tls->bufs[avail]);
    if (!(tls->bufs[avail] = calloc(1, size)))
        return NULL;
#ifdef TLS_DEBUG
    fprintf(stderr, "Alloc %d: %ld = %p\n", avail, size, tls->bufs[avail]);
#endif
    tls->sizes[avail] = size;
    tls->used[avail] = 1;

    return tls->bufs[avail];
}

void *htscodecs_tls_calloc(size_t nmemb, size_t size) {
#ifdef TLS_DEBUG
    fprintf(stderr, "htscodecs_tls_calloc(%ld)\n", nmemb*size);
#endif
    void *ptr = htscodecs_tls_alloc(nmemb * size);
    if (ptr)
        memset(ptr, 0, nmemb * size);
    return ptr;
}

void htscodecs_tls_free(void *ptr) {
    if (!ptr)
        return;

    tls_pool *tls = pthread_getspecific(rans_key);

    int i;
    for (i = 0; i < MAX_TLS_BUFS; i++) {
        if (tls->bufs[i] == ptr)
            break;
    }
#ifdef TLS_DEBUG
    fprintf(stderr, "Fake free %d size %ld ptr %p\n",
            i, tls->sizes[i], tls->bufs[i]);
#endif
    if (i == MAX_TLS_BUFS) {
        fprintf(stderr, "Attempt to htscodecs_tls_free a buffer not allocated"
                " with htscodecs_tls_alloc\n");
        return;
    }
    if (!tls->used[i]) {
        fprintf(stderr, "Attempt to htscodecs_tls_free a buffer twice\n");
        return;
    }
    tls->used[i] = 0;
}

#else
/*
 * Calloc/free equivalents instead.
 *
 * We use calloc instead of malloc as a sufficiently malformed set of input
 * frequencies may not sum to the expected total frequency size, leaving
 * some elements uninitialised.  It's unlikely, but potentially a crafty
 * attacker could somehow exploit this to pull out parts of this allocated
 * buffer and leak them into the decompressed data stream, potentially
 * compromising previous buffers such as encryption keys.  (Although
 * frankly any well-written crypto library should be zeroing such memory
 * before freeing it to ensure it's never visible to a subsequent malloc.)
 */
void *htscodecs_tls_alloc(size_t size) {
    return calloc(1, size);
}

void *htscodecs_tls_calloc(size_t nmemb, size_t size) {
    return calloc(nmemb, size);
}

void htscodecs_tls_free(void *ptr) {
    free(ptr);
}
#endif

/*
 * Given a compressed block of data in a specified compression method,
 * fill out the 'cm' field with meta-data gleaned from the compressed
 * block.
 *
 * If comp is HTS_COMP_UNKNOWN, we attempt to auto-detect the compression
 * format, but this doesn't work for all methods.
 *
 * Retuns the detected or specified comp method, and fills out *cm
 * if non-NULL.
 */
enum hts_comp_method hts_expand_method(uint8_t *data, int32_t size,
                                       enum hts_comp_method comp,
                                       hts_comp_method_t *cm) {
    if (cm) {
        memset(cm, 0, sizeof(*cm));
        cm->method = comp;
    }

    const char *xz_header = "\xFD""7zXZ"; // including nul

    if (comp == HTS_COMP_UNKNOWN) {
        // Auto-detect
        if (size >= 2 && data[0] == 0x1f && data[1] == 0x8b)
            comp = HTS_COMP_GZIP;
        else if (size >= 3 && data[1] == 'B' && data[2] == 'Z'
                 && data[3] == 'h')
            comp = HTS_COMP_BZIP2;
        else if (size >= 6 && memcmp(xz_header, data, 6) == 0)
            comp = HTS_COMP_LZMA;
        else
            return HTS_COMP_UNKNOWN;
    }

    if (!cm)
        return comp;

    // Interrogate the compressed data stream to fill out additional fields.
    switch (comp) {
    case HTS_COMP_GZIP:
        if (size > 8) {
            if (data[8] == 4)
                cm->level = 1;
            else if (data[8] == 2)
                cm->level = 9;
            else
                cm->level = 5;
        }
        break;

    case HTS_COMP_BZIP2:
        if (size > 3 && data[3] >= '1' && data[3] <= '9')
            cm->level = data[3]-'0';
        break;

    case HTS_COMP_RANS4x8:
        cm->Nway = 4;
        if (size > 0 && data[0] == 1)
            cm->order = 1;
        else
            cm->order = 0;
        break;

    case HTS_COMP_RANSNx16:
        if (size > 0) {
            cm->order  = data[0] & 1;
            cm->Nway   = data[0] & RANS_ORDER_X32    ? 32 : 4;
            cm->rle    = data[0] & RANS_ORDER_RLE    ?  1 : 0;
            cm->pack   = data[0] & RANS_ORDER_PACK   ?  1 : 0;
            cm->cat    = data[0] & RANS_ORDER_CAT    ?  1 : 0;
            cm->stripe = data[0] & RANS_ORDER_STRIPE ?  1 : 0;
            cm->nosz   = data[0] & RANS_ORDER_NOSZ   ?  1 : 0;
        }
        break;

    case HTS_COMP_ARITH:
        if (size > 0) {
            // Not in a public header, but the same transforms as rANSNx16
            cm->order  = data[0] & 3;
            cm->rle    = data[0] & RANS_ORDER_RLE    ?  1 : 0;
            cm->pack   = data[0] & RANS_ORDER_PACK   ?  1 : 0;
            cm->cat    = data[0] & RANS_ORDER_CAT    ?  1 : 0;
            cm->stripe = data[0] & RANS_ORDER_STRIPE ?  1 : 0;
            cm->nosz   = data[0] & RANS_ORDER_NOSZ   ?  1 : 0;
            cm->ext    = data[0] & 4 /*external*/    ?  1 : 0;
        }
        break;

    case HTS_COMP_TOK3:
        if (size > 8) {
            if (data[8] == 1)
                cm->level = 11;
            else if (data[8] == 0)
                cm->level = 1;
        }
        break;

    default:
        break;
    }

    return comp;
}

/*
 * A short single-letter code associated with the expanded compression
 * method cm.
 */
char hts_comp_method_short(hts_comp_method_t *cm) {
    switch (cm->method) {
    case HTS_COMP_RAW:
        return '.';

    case HTS_COMP_GZIP:
        switch (cm->level) {
        case 1:  return '_';
        case 9:  return 'G';
        }
        return 'g';

    case HTS_COMP_BZIP2:
        return cm->level >= 9 ? 'B' : 'b';

    case HTS_COMP_LZMA:
        return 'l';

    case HTS_COMP_RANS4x8:
        return cm->order ? 'R' : 'r';

    case HTS_COMP_RANSNx16: {
        char c = cm->order ? '1' : '0';
        c += (cm->Nway == 32)*4;
        if (cm->stripe) c = '8';
        if (cm->cat) c = '2';
        return c;
    }

    case HTS_COMP_ARITH:
        return cm->order ? 'A' : 'a';

    case HTS_COMP_FQZ:
        return 'f';

    case HTS_COMP_TOK3:
        return cm->level >= 10 ? 'N' : 'n';

    default:
        break;
    }

    return '?';
}

/*
 * Fills out a string associated with the expanded compression
 * method cm.  The string should be at least HTS_COMP_METHOD_STR_SIZE
 * bytes long.
 */
void hts_comp_method_long(hts_comp_method_t *cm, char *str) {
    switch (cm->method) {
    case HTS_COMP_RAW:
        strcpy(str, "raw");
        break;

    case HTS_COMP_GZIP:
        strcpy(str, "gzip");
        if (cm->level == 1)
            strcat(str, "-min");
        else if (cm->level >= 9)
            strcat(str, "-max");
        break;

    case HTS_COMP_BZIP2:
        sprintf(str, "bzip2-%d", cm->level);
        break;

    case HTS_COMP_LZMA:
        strcpy(str, "lzma");
        break;

    case HTS_COMP_RANS4x8:
        sprintf(str, "rans4x8-o%d", cm->order);
        break;

    case HTS_COMP_RANSNx16: {
        if (cm->cat) {
            strcpy(str, "ransNx16-cat");
            break;
        }
        if (cm->stripe) {
            strcpy(str, "ransNx16-stripe");
            break;
        }

        if (cm->Nway == 32)
            sprintf(str, "rans32x16-o%d", cm->order);
        else
            sprintf(str, "rans4x16-o%d", cm->order);
        break;
    }

    case HTS_COMP_ARITH:
        if (cm->cat) {
            strcpy(str, "arith-cat");
            break;
        }
        if (cm->stripe) {
            strcpy(str, "arith-stripe");
            break;
        }
        if (cm->ext) {
            strcpy(str, "arith-ext");
            break;
        }

        sprintf(str, "arith-o%d", cm->order);
        break;

    case HTS_COMP_FQZ:
        strcpy(str, "fqzcomp");
        break;

    case HTS_COMP_TOK3:
        if (cm->level >= 10)
            strcpy(str, "tok3-arith");
        else
            strcpy(str, "tok3-rans");
        break;

    default:
        strcpy(str, "?");
        break;
    }
}
