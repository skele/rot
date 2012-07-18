/*
 * Copyright (c) 2008, David J. Vanecek.  All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * 3. The name of the author may not be used to endorse or promote
 * products derived from this software without specific prior written
 * permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 * GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 * IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
 * OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include <sys/types.h>
#include <sys/mman.h>		       /* must follow sys/types */
#include <sys/wait.h>
#include <sys/resource.h>
#include <sys/time.h>

#include <assert.h>
#include <err.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>


#include "rot.h"


#define INCR 20000

static void *
exparray(struct data * a, int *nmax, size_t * totalloc)
{
    void *r;
    size_t newalloc;

    if (a == NULL) {
	/*
	 * This is the initial allocation
	 */
	*totalloc = 0;
    }
    newalloc = *totalloc + INCR * sizeof(struct data);
    r = realloc(a, newalloc);
    if (r == NULL)
	err(1, "could not build dynamic array %lu",
	    newalloc / sizeof(double) / 3);
    *totalloc = newalloc;
    *nmax += INCR;
    fprintf(stderr, "return from exparray, r is %p, %ld\n", r, *totalloc);
    fprintf(stderr, "records allocated %ld\n",
	    (*totalloc) / sizeof(struct data));
    return r;
}

int
readdblm(FILE * in, struct data * d, int nd, int nskip, int nread)
{

/*
 * Read in at most nd records (the count includes comments).
 *                             ^^^^  WHAT DOES THIS MEAN??
 *
 * Each record is either a comment (# in column 1) or a data
 * record.
 *
 * A data record consists of floating point fields separated by white space.
 * (See the constructed format string for the gritty details.
 *
 * Return the actual number of data records read.
 */

    int i;
    int c;
    int nrec;
    size_t formatlen;
    size_t st;
    char *format;

    /*
     * Construct the format.  It consists of nskip fields "%*lf " and nread
     * fields "%lf ", and a terminating null.  Its size is, therefore,
     * nskip*5 + nread*4 + 1;
     */

    formatlen = (size_t) (nskip * 5 + nread * 4 + 1) * sizeof(char);
    if ((format = malloc(formatlen)) == NULL)
	err(1, "readdbl");

    *format = '\0';
    st = 0;
    for (i = 0; i < nskip; i++) {
	st = strlcat(format, "%*lf ", formatlen);
	if (st >= formatlen)
	    err(1, "internal error in readdbl %lu %lu", st, formatlen);
    }
    for (i = 0; i < nread; i++) {
	st = strlcat(format, "%lf ", formatlen);
	if (st >= formatlen)
	    err(1, "internal error in readdbl %lu %lu", st, formatlen);
    }

    /* Fix Aug 27, 2011:  Format has spurious trailing blank.  Remove it. */

    format[strlen(format)-1]='\0';

    /* End fix of Aug 27, 2011 */


    nrec = 0;
    for (i = 0; i < nd; i++) {	       /* Store input data in dynamic array */
	c = getc(in);
	if (c != '#') {
	    if (c == EOF)
		break;
	    ungetc(c, in);
	    if ((c = fscanf(in, format,
			    &d[nrec].m, &d[nrec].x, &d[nrec].y, &d[nrec].z,
			    &d[nrec].vx, &d[nrec].vy, &d[nrec].vz)) == EOF)
		break;
	    else
		nrec++;
	}
	for (;;) {		       /* drain record */
	    c = getc(in);
	    if (c == EOF || c == '\n')
		break;
	}
	if (c == EOF)
	    break;
	if (ferror(in))
	    err(1, "IO error in readdblm");
	d[i].i = i;		       /* this is useful for sorting, after a
				        * fashion */
    }
    free(format);
    return nrec;
}


struct data *
OLDDEADreaddbl(FILE * in, int *np, int nskip, int nread)
{
/*
 * OBSOLETE ROUTINE.  NOT USED.
 * POSSIBLY BROKEN.
 */

/*
 * This routine reads in a file of data, consiting of nskip+nread+unknown
 * floating point fields.  It reads the data into a dynamically allocated
 * array, which is initialized in this routine.
 * When all the data in the file is read, it returns the address of
 * the array.
 */

    int i;
    int c;
    int nmax;
    size_t totalloc;
    size_t formatlen;
    size_t st;

    struct data *p;
    char *format;

    /*
     * Construct the format.  It consists of nskip fields "%*lf " and nread
     * fields "%lf ", and a terminating null.  Its size is, therefore,
     * nskip*5 + nread*4 + 1;
     */

    formatlen = (size_t) (nskip * 5 + nread * 4 + 1) * sizeof(char);
    if ((format = malloc(formatlen)) == NULL)
	err(1, "readdbl");

    *format = '\0';
    st = 0;
    for (i = 0; i < nskip; i++) {
	st = strlcat(format, "%*lf ", formatlen);
	if (st >= formatlen)
	    err(1, "internal error in readdbl %lu %lu", st, formatlen);
    }
    for (i = 0; i < nread; i++) {
	st = strlcat(format, "%lf ", formatlen);
	if (st >= formatlen)
	    err(1, "internal error in readdbl %lu %lu", st, formatlen);
    }

    /* Fix Aug 27, 2011:  Format has spurious trailing blank.  Remove it. */

    format[strlen(format)-1]='\0';

    /* End fix of Aug 27, 2011 */

    nmax = 0;
    p = exparray(NULL, &nmax, &totalloc);
    *np = 0;
    for (i = 0;; i++) {		       /* Store input data in dynamic array */

	/*
	 * See if we're going to EOF, if so get out now, rather than possibly
	 * allocating an unneeded block.
	 */

	c = getc(in);
	if (c == EOF)
	    break;
	ungetc(c, in);

	if (i >= nmax)
	    p = exparray(p, &nmax, &totalloc);
	if ((c = fscanf(in, format,
			&p[i].m, &p[i].x, &p[i].y, &p[i].z,
			&p[i].vx, &p[i].vy, &p[i].vz)) == EOF)
	    break;
	for (;;) {		       /* drain record */
	    c = getc(in);
	    if (c == EOF || c == '\n')
		break;
	}
	if (c == EOF)
	    break;
	if (ferror(in))
	    err(1, "IO error in readdbl, total allocation %lu",
		(unsigned long) totalloc);
	p[i].i = i;		       /* this is useful for sorting, after a
				        * fashion */
	if (*np != i)
	    errx(1, "internal error in readdbl, %d != %d", *np, i);
	(*np)++;
    }
    free(format);
    fprintf(logfile, "Rough allocation, %d records %lu bytes\n", *np, totalloc);
    p = realloc(p, (size_t) (*np) * sizeof(struct data));
    fprintf(logfile, "Trimmed to %lu bytes, %d records of %lu bytes each.\n",
	    (*np) * sizeof(struct data), *np, sizeof(struct data));
    if (p == NULL)
	err(1, "trimming p failed in readdbl");
    return p;
}
