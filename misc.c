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

/*
 * Notes
 *
 * Unit vectors
 *   Spherical
 *	e_r:	( cos(phi)*sin(theta), sin(phi)*sin(theta),  cos(theta))
 *	e_theta:( cos(phi)*cos(theta), sin(phi)*cos(theta), -sin(theta))
 *	e_phi:	(           -sin(phi),            cos(phi),          0 )
 *
 *   Cylindrical
 *      e_rho:	(            cos(phi),            sin(phi),          0 )
 *	e_phi:	(           -sin(phi),            cos(phi),          0 )
 *      e_z:	(                   0,                   0,          1 )
 *
 */

void
xlate(struct data * d, int nd, double x0, double y0, double z0)
{
/*
 * Make (x0,y0,z0) the center of the coordinate system.
 */
    int i;

    for (i = 0; i < nd; i++) {
	d[i].x -= x0;
	d[i].y -= y0;
	d[i].z -= z0;
    }
    return;
}

void
derive(struct data * d, int nd)
{
/*
 * given (x,y,z) and (vx,vy,vz), calculate various transformed
 * quantities.
 */

    int i;

    for (i = 0; i < nd; i++) {
	d[i].rho = hypot(d[i].x, d[i].y);	/* cylindrical distance */
	d[i].r = hypot(d[i].rho, d[i].z);	/* radial distance */
	d[i].theta = atan2(d[i].rho, d[i].z);	/* zenith angle */
	d[i].phi = atan2(d[i].y, d[i].x);	/* azimuth angle */

	d[i].csphi = cos(d[i].phi);    /* why not x/rho? */
	d[i].snphi = sin(d[i].phi);
	d[i].sntheta = sin(d[i].theta);
	d[i].cstheta = cos(d[i].theta);

	d[i].vrho = dot(d[i].vx, d[i].vy, d[i].vz,
			d[i].csphi, d[i].snphi, 0.);

	d[i].vr = dot(d[i].vx, d[i].vy, d[i].vz,
		      d[i].csphi * d[i].sntheta, d[i].snphi * d[i].sntheta,
		      d[i].cstheta);

	d[i].vtheta = dot(d[i].vx, d[i].vy, d[i].vz,
		       d[i].csphi * d[i].cstheta, d[i].snphi * d[i].cstheta,
			  -d[i].sntheta);

	d[i].vphi = dot(d[i].vx, d[i].vy, d[i].vz,
			-d[i].snphi, d[i].csphi, 0.);
    }
}

void
cross(double *c1, double *c2, double *c3, double a1, double a2, double a3,
      double b1, double b2, double b3, double f)
{
/*  C = A cross B */

    *c1 = f * (a2 * b3 - a3 * b2);
    *c2 = f * (a3 * b1 - a1 * b3);
    *c3 = f * (a1 * b2 - a2 * b1);
    return;
}

double
dot(double a1, double a2, double a3, double b1, double b2, double b3)
{
/*  return A dot B  */
    return a1 * b1 + a2 * b2 + a3 * b3;
}


double
powi(double x, int i)
{
/* Linux style :-(
 *
 * This may be marginally better than calling pow(3m)
 */

    double xx;

    if (i > 10)
	return pow(x, (double) i);

    if (i == 0 && x == 0.)
	return x / x;
    if (i == 0)
	return 1.;
    if (x == 0.)
	return 0.;

    if (i < 0) {
	x = 1. / x;
	i = -i;
    }
    xx = x;
    for (--i; i; i--)
	x *= xx;
    return x;
}


static struct rusage rusage, rusagec;
static time_t tzero, lasttime, thistime;	/* Wall */

void
initclock(void)
{
    tzero = time(&tzero);
}

void
clocker(char *file, int line)
{
    struct funky {
	long tv_sec;		       /* seconds since Jan. 1, 1970 */
	long tv_usec;		       /* and microseconds */
    };



    getrusage(RUSAGE_SELF, &rusage);
    getrusage(RUSAGE_CHILDREN, &rusagec);
    thistime = time(&thistime) - tzero;
    fprintf(stderr, "%s:%4d: %ld.%3.3ldu %ld.%3.3lds %ldMflt %ldSwp %ldIO", file, line,
	    rusage.ru_utime.tv_sec, rusage.ru_utime.tv_usec / 1000,
	    rusage.ru_stime.tv_sec, rusage.ru_stime.tv_usec / 1000,
	    rusage.ru_majflt, rusage.ru_nswap, rusage.ru_inblock + rusage.ru_oublock);
    fprintf(stderr, "\tChildren: %ld.%3.3ldu %ld.%3.3lds %ldMflt %ldSwp %ldIO\tWall: %u\n",
	    rusagec.ru_utime.tv_sec, rusagec.ru_utime.tv_usec / 1000,
	    rusagec.ru_stime.tv_sec, rusagec.ru_stime.tv_usec / 1000,
	    rusagec.ru_majflt, rusagec.ru_nswap, rusagec.ru_inblock + rusagec.ru_oublock,
	    thistime);
    lasttime = thistime;
    return;
}

int
compr(const void *a, const void *b)
{
    if (((struct data *) a)->r < ((struct data *) b)->r)
	return -1;
    if (((struct data *) a)->r > ((struct data *) b)->r)
	return 1;
    return 0;
}

int
comptheta(const void *a, const void *b)
{
/*
 * theta is an angle, and there should be some recognition
 * of 2pi here.
 */
    if (((struct data *) a)->theta < ((struct data *) b)->theta)
	return -1;
    if (((struct data *) a)->theta > ((struct data *) b)->theta)
	return 1;
    return 0;
}

int
compq(const void *a, const void *b)
{
/* sort largest to top of array */
    if (((struct qtype *) a)->dsq < ((struct qtype *) b)->dsq)
	return 1;
    if (((struct qtype *) a)->dsq > ((struct qtype *) b)->dsq)
	return -1;
    return 0;
}

int 
fpshellhdr(FILE * f)
{
    return fprintf(f, "#\t\tib, n, r, dr, theta, dtheta, bound, volume");
}


int
fpshell(FILE * f, const char *title, struct shell_t * s)
{
    int r;

    r = fprintf(f, "%12.12s: %12u %12u %12.7e %12.7e %12.7e %12.7e %12.7e %12.7e\n",
		title,
		s->ib, s->n, s->r, s->dr, s->theta, s->dtheta, s->bound,
		s->volume);
    return r;
}



int
binary_search(struct qtype U[], int low, int high, double target)
{
/*
 * Binary search -- returns an index into U[].
 *
 * Note, this is not the typical bsearch code, which demands that
 * target not exist in the list.
 *
 * REALLY note: this code wants U sorted in DECREASING ORDER.
 *
 * If target is found at U[k], return k.
 * If target is not found, return k such that U[k] > target
 *    and U[k+1] < target,  if U[k+1] exists.  If target > U[0],
 *    return -1, note that U[0 = -1+1] is still the first entry
 *    < target.
 * SO: target can be inserted after the returned index and not
 * disturb the sort-order of U.
 */

    int middle;

    middle = low + (high - low) / 2;
    while (low <= high) {
	middle = low + (high - low) / 2;
	if (target > U[middle].dsq)
	    high = middle - 1;
	else if (target < U[middle].dsq)
	    low = middle + 1;
	else
	    return middle;
    }
    if (target < U[middle].dsq)
	return middle;
    return middle - 1;
}

/* Dead code, held in reserve against a better day?

int
searchandinsert(struct qtype U[], struct qtype qtemp, int jmax)
{
    int r;
    int k;

    r = binary_search(U, 0, jmax - 1, qtemp.dsq);
    if (r == jmax)
	r--;
    if (r == 0)
	U[0] = qtemp;
    else {
	for (k = 0; k < r; k++) {
	    U[k] = U[k + 1];
	}
	/ * memcpy(U, U+1, (r)*sizeof(*U)); * /
	U[r] = qtemp;
    }
    return 0;
}

*/
