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

/* ISO 31-11 */


#include <err.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "rot.h"


static struct shell_t *mktorus_v(struct data *, int, struct shell_t, int);


struct shell_t *
mkshell_p(struct data * d, int nd, int nshell)
{
/*
 * partition the data in <d>, (<nd> points), into subdomains based
 * on rho, the distance from (0,0,0).  There will be <nshell> of these
 * subdomains, henceforth called "shells".
 *
 * The input data is assumed sorted on increasing rho.
 *
 * How the partitioning is done might vary.  Plausible methods:
 *
 * each shell has the same drho.
 * each shell has the same volume.
 * each shell has an equal (ignoring remainders) number of particles.
 * each shell has arbitrary boundaries.
 *
 */

/*
 * This version uses an equal number of particles per shell.
 * Remainder particles go in the last (outermost) shells.
 */

    int extra;
    int i;
    int ntot;
    int pershell;
    struct shell_t *s;


    pershell = nd / nshell;
    if (pershell == 0)
	errx(1, "Too many shells (%d) for %d data points", nshell, nd);

    if ((s = calloc(nshell, sizeof(struct shell_t))) == NULL)
	err(1, "mkshell");
    memset(s, 0, nshell * sizeof(struct shell_t));

    /*
     * s[i].n is the number of particles in the i-th shell.
     * s[i].ib is the first index (into <d>) in the i-th shell.
     * s[i].bound is the outer boundary of the i-th shell.
     * s[i].drho is the thickness of the i-th shell.
     */

    extra = nd % nshell;
    extra = nshell - extra;	       /* distribute these over the last
				        * shells */
    ntot = 0;
    for (i = 0; i < nshell; i++) {
	s[i].ib = ntot;
	s[i].n = pershell;
	if (i >= extra)
	    s[i].n++;
	ntot += s[i].n;
	s[i].r = i == 0 ? 0. : (s[i - 1].r + s[i - 1].dr);
	s[i].bound = (d[s[i].ib + s[i].n - 1].r + d[s[i].ib + s[i].n].r) / 2.;
	s[i].dr = s[i].bound - s[i].r;
	s[i].volume = powi(s[i].r + s[i].dr, 3) - powi(s[i].r, 3);
	s[i].volume *= 4. / 3. * M_PI;
    }

    /*
     * The last shell gets special treatment
     *
     * Leave the ternaries in this block alone, there is a special case that
     * they fit: one big shell.
     */

    i--;

    s[i].bound = (3. * d[nd - 1].r - d[nd - 2].r) / 2.;
    s[i].dr = s[i].bound - (i == 0 ? 0. : s[i - 1].bound);
    s[i].volume = 4. / 3. * M_PI * (powi(s[i].bound, 3) -
				    (i == 0 ? 0. : powi(s[i - 1].bound, 3)));

    return s;
}

struct shell_t *
mktorus_p(struct data * d, int nd, struct shell_t s, int nt)
{
/*
 * Given a single shell, s, of data d, divide it into nt tori.
 *
 * <d[s.ib]> is the first datum in the shell <s>.
 */

/*
 * the partitioning is based on equal number of particles
 * per torus.
 */

/*
 * V(i) = 2*PI/3*(r(i+1)^3-r(i)^3)*(cos(theta(i)-cos(theta(i+1))))
 * source:  derived from the triple integral
 */

    int extra;
    int i;
    int ntot;
    int pertorus;
    struct shell_t *t;


    pertorus = s.n / nt;
    if (pertorus == 0) {
	warnx("Too many tori (%d) for %d data points in this shell", nt, s.n);
	warnx("Changing to equivolume method for this shell at r= %e\n", s.r);
	return mktorus_v(d, nd, s, nt);
    }
    /*
     * We might try reducing nt and trying again until nt == 0 The shape of
     * the cluster is going to have a big impact. if it is globular, then
     * we're ok.  But suppose it is a disc? Probably, for non-globulars, we
     * should go to an equi-volume partitioning for both the shells and the
     * tori.
     */
    if ((t = calloc(nt, sizeof(struct shell_t))) == NULL)
	err(1, "mktorus");
    memset(t, 0, nt * sizeof(struct shell_t));

    /* First sort the shell on theta */

    qsort(d + s.ib, s.n, sizeof(struct data), comptheta);

    extra = s.n % nt;
    extra = nt - extra;		       /* distribute these over the last
				        * tori. */
    ntot = 0;
    for (i = 0; i < nt; i++) {
	t[i].n = pertorus;
	t[i].ib = s.ib + ntot;
	if (i >= extra)
	    t[i].n++;
	ntot += t[i].n;
	/* "bound" could be improved: */
	t[i].bound = (d[t[i].ib + t[i].n - 1].theta + d[t[i].ib + t[i].n].theta) / 2.;

	t[i].r = s.r;
	t[i].dr = s.dr;
	t[i].theta = i == 0 ? 0. : (t[i - 1].theta + t[i - 1].dtheta);
	t[i].dtheta = t[i].bound - t[i].theta;

	t[i].volume = powi(t[i].r + t[i].dr, 3) - powi(t[i].r, 3);
	t[i].volume *= cos(t[i].theta) - cos(t[i].theta + t[i].dtheta);
	t[i].volume *= 2. * M_PI / 3.;
    }

    i--;

    t[i].bound = M_PI;
    t[i].dtheta = M_PI - t[i].theta;
    t[i].volume = powi(t[i].r + t[i].dr, 3) - powi(t[i].r, 3);
    t[i].volume *= cos(t[i].theta) - cos(t[i].theta + t[i].dtheta);
    t[i].volume *= 2. * M_PI / 3.;

    return t;
}

struct shell_t *
mktorus_v(struct data * d, int nd, struct shell_t s, int nt)
{
/*
 * Given a single shell, s, of data d, divide it into nt tori.
 *
 * <d[s.ib]> is the first datum in the shell <s>.
 */

/*
 * the partitioning is based on equal volumes per torus.
 * the number of tori is specified in argument <nt>.
 */

/*
 * V(i) = 2*PI/3*(r(i+1)^3-r(i)^3)*(cos(theta(i)-cos(theta(i+1))))
 * source:  derived from the triple integral
 */

    int i;
    int ntot;
    struct shell_t *t;
    double theta;
    double cstheta;

    if ((t = calloc(nt, sizeof(struct shell_t))) == NULL)
	err(1, "mktorus");
    memset(t, 0, nt * sizeof(struct shell_t));

    /* First sort the shell on theta */

    qsort(d + s.ib, s.n, sizeof(struct data), comptheta);

    ntot = 0;
    theta = 0.;
    cstheta = 1.;

    for (i = 0; i < nt; i++) {
	t[i].r = s.r;
	t[i].dr = s.dr;
	t[i].theta = theta;
	t[i].dtheta = -theta;
	cstheta -= 2. / (double) nt;
	if (cstheta < -1.)
	    cstheta = -1.;	       /* Rounding error has been observed */
	theta = acos(cstheta);

	/*
	 * printf("i: %d cstheta: %e acos: %e theta: %e  cstheta+1: %e\n", i,
	 * cstheta, acos(cstheta), theta, cstheta+1.);
	 */
	t[i].dtheta += theta;
	t[i].bound = theta;
	t[i].volume = s.volume / (double) nt;
	t[i].ib = s.ib + ntot;
	t[i].n = 0;
	while (ntot < s.n) {
	    /*
        printf("mktorus_v:i=%d Test arg: [%d].theta %e theta %e test: %s diff: %e csth: %e\n", i,
	       t[i].ib + t[i].n, d[t[i].ib + t[i].n].theta, theta,
	       (d[t[i].ib + t[i].n].theta > theta) ? "TRUE" : "FALSE",
	       (d[t[i].ib + t[i].n].theta - theta),
	       cstheta);
	       */
	    if (d[t[i].ib + t[i].n].theta > theta)
		break;
	    t[i].n++;
	    ntot++;
	}
    }
    return t;
}

struct shell_t *
mkshell(struct data * d, int nd, int nshell)
{
    return mkshell_p(d, nd, nshell);
}

struct shell_t *
mktorus(struct data * d, int nd, struct shell_t s, int nt)
{
    return mktorus_p(d, nd, s, nt);
}
