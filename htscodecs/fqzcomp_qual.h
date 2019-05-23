/*
 * Copyright (c) 2011-2013, 2018-2019 Genome Research Ltd.
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

#ifndef FQZ_COMP_QUAL_H
#define FQZ_COMP_QUAL_H

#ifdef __cplusplus
extern "C" {
#endif

/* Bit flags, deliberately mirroring BAM ones */
#define FQZ_FREVERSE 16
#define FQZ_FREAD2 128

/* Maximum length of an individual quality string.
 * If longer than this, simply break it up into
 * smaller portions.
 */
#ifndef MAX_SEQ
#  define MAX_SEQ 100000
#endif

/* A single record.
 * To compress we need to know the junction from one quality string to
 * the next (len, qual), whether it is first/second read and whether it is
 * reverse complemented (flags).
 */
typedef struct {
    int len;
    int qual; // FIXME: merge len and qual.  Artificial
    int flags;
} fqz_rec;

typedef struct {
    int num_records;
    fqz_rec *crecs;
} fqz_slice;

char *fqz_compress(int vers, fqz_slice *s, char *in, size_t uncomp_size,
                   size_t *comp_size, int level);
char *fqz_decompress(char *in, size_t comp_size, size_t *uncomp_size, int *lengths);

#ifdef __cplusplus
}
#endif

#endif
