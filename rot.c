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

/*
 * Todo: write a safe "realloc" with the calling convention of calloc(3),
 * checking for overflow, etc etc.
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

#ifdef LINUX
#include <values.h>
#include "linuxlib.c"
#endif

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

#undef DOCLOCK

#include "rot.h"

extern int comprho(const void *, const void *);

/************************************************************************/

int inserts_at_0 = 0;
int inserts_elsewhere = 0;
int non_inserts = 0;
int easy_inserts = 0;

int
neighbors(int ikid, int istart, int icount, va_list vargs)
{
    int i, j, k, r, nd, jmax;
    double radius, density;/* of "neighborhoods" */
    float mass;
    float *tx, *ty, *tz, *tradius, *tdens;
    struct data *d;
    struct qtype *Q;
    struct qtype qtemp;

 
    d = va_arg(vargs, struct data *);
    Q = va_arg(vargs, struct qtype *);
    nd = va_arg(vargs, int);
    jmax = va_arg(vargs, int);

    tx = (float *) malloc(nd*(sizeof( float)));
    ty = (float *) malloc(nd*(sizeof( float)));
    tz = (float *) malloc(nd*(sizeof( float)));
    tdens = (float *) malloc(nd*(sizeof( float)));
    tradius = (float *) malloc(nd*(sizeof( float)));

    inserts_at_0 = inserts_elsewhere = non_inserts = easy_inserts = 0;

    //set mass to some mass, assume all masses are equal
    mass = d[0].m;
    obtain_densities(tx,ty,tz,tradius,tdens,istart,icount,ikid,nd,mass);

    for (i = istart; i < istart + icount; i++) {
      d[i].radius = tradius[i];
      d[i].density = tdens[i];
    }

    /*
     * printf("Inserts at 0 %u, elsewhere %u, easy_inserts %u, non-inserts
     * %u\n", inserts_at_0, inserts_elsewhere, easy_inserts, non_inserts);
     */
    return 0;
}

int
neighbor_dispersion(int ikid, int istart, int icount, va_list vargs)
{
/*  Calculate the dispersion of radial velocity for each "neighborhood" */

    struct data *d;
    int jmax;

    int i;
    int j;
    double v, vsq;
    double tm;

    d = va_arg(vargs, struct data *);
    jmax = va_arg(vargs, int);

    for (i = istart; i < istart + icount; i++) {
	v = vsq = tm = 0.;
	for (j = 0; j < jmax; j++) {
	    v += d[d[i].nn[j].idx].vphi;
	    vsq += d[d[i].nn[j].idx].vphi * d[d[i].nn[j].idx].vphi;
	    tm += d[d[i].nn[j].idx].m;
	}
	v += d[i].vphi;

	vsq += d[i].vphi * d[i].vphi;
	tm += d[i].m;
	d[i].sigma2neigh = 1. / (jmax + 1) * (vsq - v * v / (jmax + 1));
	d[i].vphineigh = v / (jmax + 1);
	d[i].Rneigh = 0.5 * tm * vsq;
	d[i].Tneigh = 0.5 * tm * d[i].sigma2neigh;
    }
    return 0;
}

void
sortmem(struct data * d, int n, int jmax)
{
/*
 * <d> has been sorted.  Renumber the nearest neigbors (d[].nn)
 * to account for that.
 */

    /* old element i maps to to[i] */

    int *to;
    int i, j;

    if ((to = calloc(n, sizeof(int))) == NULL)
	err(1, "sortmem");

    for (i = 0; i < n; i++)
	to[d[i].i] = i;

    for (i = 0; i < n; i++)	       /* put the complexity into the
				        * indexing scheme */
	for (j = 0; j < jmax; j++)
	    d[i].nn[j].idx = to[d[i].nn[j].idx];

/*  These checks are not generally useful; they
 *  were present during development.
 *
 *    for (i = 0; i < n; i++)
 *	if (memarr[i].radius != d[i].radius)
 *	    errx(1, "sortmem failed at %d %d %e %e", i, d[i].i, memarr[i].radius, d[i].radius);
 *
 *    for (i = 0; i < n; i++) {
 *	double dsq;
 *	int k;
 *
 *	for (j = 0; j < jmax; j++) {
 *
 *	    k = memarr[i].nn[j].idx;
 *
 *	    dsq = SQD(d[i].x, d[k].x);
 *	    dsq += SQD(d[i].y, d[k].y);
 *	    dsq += SQD(d[i].z, d[k].z);
 *	    if (dsq != memarr[i].nn[j].dsq)
 * The reason this check fails by small amounts is the roundoff/truncation
 * error introduced by translation of d[x,y,z] to the center of mass.
 *		if (fabs((memarr[i].nn[j].dsq - dsq) / dsq) > 3.*DBL_EPSILON)
 *		    fprintf(logfile, "distance discrepancy: %6d %6d %6d %.12e %.12e %.12e %.12e %9e epsils\n",
 *			    i, j, k, memarr[i].nn[j].dsq, dsq,
 *			    memarr[i].nn[j].dsq / dsq,
 *			    memarr[i].nn[j].dsq - dsq,
 *			    (memarr[i].nn[j].dsq - dsq) / dsq /DBL_EPSILON);
 *	}
 *    }
 */
    free(to);
    return;
}

void
renumber(struct data * d, int nd)
{
/*
 * Does what it does.  This .i preserves the index of d before a
 * sort.  I.e. after a sort, .i can be used to adjust the indices
 * of "other things", chiefly the nearest-neighbors in d[].nn.
 */
    int i;

    for (i = 0; i < nd; i++)
	d[i].i = i;
    return;
}


int
main(int argc, char **argv)
{
    int i, k;
    int istart, iend, incr;
    int printneighbors = 0;
    int nd;
    int ndread;

    size_t areasize;

    struct data *d;
    struct qtype *dv;		       /* variable part of d */
    struct qtype *Q;

    double xcen, ycen, zcen, densitysum;
    double xm, ym, zm, tm;
    double xnm, ynm, znm;
    double rden;
    double dwaden;

    FILE *patrickoutfile;
    
    patrickoutfile = fopen("rot_output_for_td.dat","w");

    nd = procopt(&argc, &argv);
    initclock();
    clockit();
    /*
     * Allocate shared memory area
     *
     * This is highly dependent on OS and physical architecture.
     * Developed on the byte-addressable i386, and OpenBSD, which
     * returns mmap areas mapped on page boundaries.  Other architectures
     * or OSes will need to look explicitly at the size of the structs
     * involved, the needs and results of mmap, whether certain internal
     * data of the structs are optimally aligned in memory.  Presumably,
     * at this late date in history, there will be no problems.
     */

    areasize = (size_t) nd *sizeof(struct data);

    fprintf(logfile, "Areasize (d): %luKB\n", areasize / 1024);
    if ((d = (struct data *) mmap(NULL, areasize, PROT_READ | PROT_WRITE,
		  MAP_ANONYMOUS | MAP_SHARED, -1, (off_t) 0)) == MAP_FAILED)
	err(1, "d-memory mapping failed");
    memset(d, 0, areasize);	       /* not needed on Linux, not specified
				        * on BSD */
    areasize = (size_t) nd *jmax * sizeof(struct qtype);

    /*
     * For some perverse reason, I use two mmapped areas.  Why we can't just
     * use one, and adjust dv accordingly, is not obvious. This should be
     * implemented in future, in the interest of cleanliness.
     */

    fprintf(logfile, "Areasize (Q): %luKB\n", areasize / 1024);
    if ((dv = (struct qtype *) mmap(NULL, areasize, PROT_READ | PROT_WRITE,
		  MAP_ANONYMOUS | MAP_SHARED, -1, (off_t) 0)) == MAP_FAILED)
	memset(dv, 0, areasize);

    for (i = 0; i < nd; i++)
	d[i].nn = dv + i * jmax;       /* this is a bit of possibly i386
				        * specific cruft */
    /* Read data into prepared d array */

    ndread = readdblm(in, d, nd, 1, 7);
    fprintf(logfile, "read in %u records, %u comments, %u data\n",
	    nd, nd - ndread, ndread);
    nd = ndread;
    derive(d, nd);		       /* Calculates "derived" quantities.
				        * See the routine */

    if (nd <= (jmax + 1))
	errx(1, "too few data (%d)", nd);
    fprintf(logfile, "Initial data, the first few records:\n");
    for (i = 0; i < (10 > nd ? nd : 10); i++) {
	fprintf(logfile,
		"%6u:\tm=%19.12e X=(%19.12e,%19.12e,%19.12e) V=(%19.12e,%19.12e,%19.12e)\n"
		"\trho=%19.12e r=%15.8e ph=%19.12e th=%19.12e\n", i,
		d[i].m, d[i].x, d[i].y, d[i].z, d[i].vx, d[i].vy, d[i].vz,
		d[i].rho, d[i].r, d[i].phi, d[i].theta);
    }

    qsort(d, nd, sizeof(struct data), compr);	/* sort on r, radial distance */

    renumber(d, nd);

    if ((Q = calloc((size_t) (jmax), sizeof(struct qtype))) == NULL)
	err(1, "could not allocate near-neighbors queue");

    incr = nd / nproc;
    istart = 0;
    iend = istart + incr;

    clockit();

    /* Do see the comments for spawnit() */

    spawnit(nproc, nd, neighbors, d, Q, nd, jmax);

    clockit();

    /*
     * This loop could go in the child processes, and the file-reassembly
     * code from earlier versions of this program restored.  The argument at
     * this time is that it is not known what sort of calculations and output
     * are wanted at this point in the program.
     */

    if (printneighbors) {
	fprintf(out, "Nearest Neighbors:\n");
	for (i = 0; i < nd; i++) {
	    fprintf(out, "%6u %u %.12e: %.12e %.12e %.12e\n",
		    i, d[i].nn[0].idx, d[i].nn[0].dsq,
		    d[d[i].nn[0].idx].r, d[i].radius, d[i].density);
	    for (k = 0; k < jmax; k++)
		fprintf(out, " %u %.5e ", d[i].nn[k].idx,
			sqrt(d[i].nn[k].dsq));
	    fprintf(out, " | %u %u", d[i].nn[0].idx,
		    d[i].nn[jmax - 1].idx);
	    putc('\n', out);
	}
    }
    clockit();

/* Calculation of center of density
* Eqn II.3
*
*
*                                  (i)              (i)
*		X    =  SUM  (X RHO   )  /  SUM  RHO
*                d,j       i   i   j           i    j
*
*    RHO is density, X is position
*/
    xcen = ycen = zcen = 0.;
    tm = xm = ym = zm = 0.;
    xnm = ynm = znm = 0.;
    densitysum = 0.;

    for (i = 0; i < nd; i++) {
      fprintf(patrickoutfile,"%10d %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e %15.7e\n", d[i].i, d[i].density, d[i].x, d[i].y
	      ,d[i].z,d[i].vx,d[i].vy,d[i].vz);
        densitysum += d[i].density;
	xcen += d[i].x * d[i].density;
	ycen += d[i].y * d[i].density;
	zcen += d[i].z * d[i].density;
	tm += d[i].m;
	xm += d[i].m * d[i].x;
	ym += d[i].m * d[i].y;
	zm += d[i].m * d[i].z;
	xnm += d[i].x;
	ynm += d[i].y;
	znm += d[i].z;
    }
    if (densitysum == 0.)
	errx(1, "density sum is zero");
    xcen /= densitysum;
    ycen /= densitysum;
    zcen /= densitysum;
    xm /= tm;
    ym /= tm;
    zm /= tm;
    xnm /= nd;
    ynm /= nd;
    znm /= nd;

/* Calculation of density radius, and density weighted average density.
*
* EQN II.4 and II.5
*
*/
    rden = 0.;
    dwaden = 0.;

    for (i = 0; i < nd; i++) {
	rden +=
	    sqrt(SQD(d[i].x, xcen) + SQD(d[i].y, ycen) +
		 SQD(d[i].z, zcen)) * d[i].density;
	dwaden += d[i].density * d[i].density;
    }
    if (densitysum == 0.)
	errx(1, "Density sum (rho sum) is zero");
    rden /= densitysum;
    dwaden /= densitysum;

    fprintf(out, "Points %u\n", nd);
    fprintf(out, "Center of density:   %.12e %.12e %.12e\n",
	    xcen, ycen, zcen);
    fprintf(out, "Center of mass:      %.12e %.12e %.12e\n",
	    xm, ym, zm);
    fprintf(out, "Center of mass (m=1):%.12e %.12e %.12e\n",
	    xnm, ynm, znm);
    fprintf(out, "Density radius:      %.12e\n",
	    rden);
    fprintf(out, "Density-weighted density: %.12e\n",
	    dwaden);

/* Rotational energy and related quantities */

    xlate(d, nd, xcen, ycen, zcen);
    derive(d, nd);
    qsort(d, (size_t) nd, sizeof(struct data), compr);
    sortmem(d, nd, jmax);
    renumber(d, nd);		       /* needed? */

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
 *      e_rho:  (            cos(phi),            sin(phi),          0 )
 *	e_phi:	(           -sin(phi),            cos(phi),          0 )
 *      e_z:	(                   0,                   0,          1 )
 *
 */


    {
	struct shell_t *s;
	struct shell_t *t;
	int is;
	int it;
	double totv;
	double T, Ts, Tt;

	double v_rot;
	double mass;
	double omega;
	double sigma_r, sigma_theta, sigma_phi;
	double sigma;
	double sigma_phiomega;
	double density;
	double center_rho;
	double center_z;
	double sig2wc, mass2;

	double ke(struct data * d);

	s = mkshell(d, nd, nshell);

	for (is = 0; is < nshell; is++) {
	    fprintf(out, "%6u: ", is);
	    fpshell(out, "Shell", s + is);
	}
	fpshellhdr(out);
	for (is = 0; is < nshell; is++) {
	    fputc('\n', out);
	    t = mktorus(d, nd, s[is], ntorus);
	    s[is].sub = t;
	    fprintf(out, "        ");
	    fpshell(out, "Shell", s + is);
	    totv = 0.;
	    for (it = 0; it < ntorus; it++) {
		fprintf(out, "%6u: ", it);
		fpshell(out, "Torus", t + it);
		totv += t[it].volume;
	    }
	    fprintf(out, "#\t\tTotal volume: %12.7e %12.7e (error)\n", totv, s[is].volume - totv);
	}

	/* sample calculation -- total kinetic energy */

	if (0) {

	    T = 0.;

	    for (is = 0; is < nshell; is++) {
		struct shell_t *tw;

		Ts = 0.;
		tw = s[is].sub;
		for (it = 0; it < ntorus; it++, tw++) {
		    Tt = 0.;
		    for (i = tw->ib; i < tw->ib + tw->n; i++) {
			Tt += ke(d + i);
		    }
		    printf("torus[%u]: Tt= %e\n", it, Tt);
		    Ts += Tt;
		}
		printf("shell[%u]: Ts= %e\n\n", is, Ts);
		T += Ts;
	    }
	    Ts = Tt = T;
	    printf("Total:       %e\n\n", T);
	    /* Check */
	    for (i = 0; i < nd; i++) {
		T -= ke(d + i);
	    }
	    printf("Error:       %e\n\n", T);
	    Tt = 0;
	    for (i = 0; i < nd; i++) {
		Tt += ke(d + i);
	    }
	    printf("Error:       %e\n\n", Tt - Ts);

	}
/* Marc/Pau stuff: */

	for (is = 0; is < nshell; is++) {
	    struct shell_t *tw;

	    tw = s[is].sub;

	    for (it = 0; it < ntorus; it++, tw++) {
		double centroid_theta, centroid_r, centroid_rho, centroid_z, centroid_y;
		double wR;

		wR = tw->r + tw->dr;
		centroid_theta = tw->theta + 0.5 * tw->dtheta;

		centroid_r = 2. / 3. * (powi(wR, 3) - powi(tw->r, 3));
		centroid_r /= (powi(wR, 2) - powi(tw->r, 2));
		centroid_r *= sin(tw->dtheta / 2.) / (tw->dtheta / 2.);

		centroid_z = centroid_r * cos(centroid_theta);
		centroid_rho = centroid_r * sin(centroid_theta);
		centroid_y = centroid_rho;

		v_rot = 0.;
		mass = 0.;
		omega = 0.;
		sigma_r = 0.;
		sigma_theta = 0.;
		sigma_phi = 0.;
		sigma_phiomega = 0.;
		center_rho = 0.;
		center_z = 0.;

		sig2wc = 0.;
		mass2 = 0.;

		if (tw->n) {

		    for (k = tw->ib; k < tw->ib + tw->n; k++) {
			double mass_k;

			mass_k = usemass ? d[k].m : 1.;
			mass2 += mass_k * mass_k;
			v_rot += mass_k * d[k].vphi;
			omega += mass_k * d[k].vphi / d[k].rho;
			mass += mass_k;
			center_rho += mass_k * d[k].rho;
			center_z += mass_k * d[k].z;

			sig2wc += mass_k * d[k].vphi * d[k].vphi;

		    }

		    sig2wc -= v_rot * v_rot / mass;

		    /*
	             * biased: sig2wc /= mass;		/ * Biased weighted
	             * variance
	             */

		    sig2wc *= mass / (mass * mass - mass2);	/* Unbiased weighted
								 * variance */
		    sig2wc = sqrt(sig2wc);

		    v_rot /= mass;
		    omega /= mass;
		    density = mass / tw->volume;
		    center_rho /= mass;
		    center_z /= mass;

		    for (k = tw->ib; k < tw->ib + tw->n; k++) {
			double mass_k;

			mass_k = usemass ? d[k].m : 1.;
			sigma_r += mass_k * d[k].vr * d[k].vr;
			sigma_theta += mass_k * d[k].vtheta * d[k].vtheta;
			sigma_phi += mass_k *
			    (d[k].vphi - v_rot) * (d[k].vphi - v_rot);
			sigma_phiomega += mass_k *
			    (d[k].vphi - omega * d[k].rho) *
			    (d[k].vphi - omega * d[k].rho);
		    }
		    sigma_r /= mass;
		    sigma_theta /= mass;
		    sigma_phi /= mass;
		    sigma_phiomega /= mass;

		    sigma = sqrt(sigma_r + sigma_theta + sigma_phi);
		    sigma_r = sqrt(sigma_r);
		    sigma_theta = sqrt(sigma_theta);
		    sigma_phi = sqrt(sigma_phi);
		    sigma_phiomega = sqrt(sigma_phiomega);
		}
		if (it == 0) {
		    fprintf(out, "# Torus:         ");
		    fprintf(out, "r          ");
		    fprintf(out, "theta      ");
		    fprintf(out, "   (centroid z,rho)   ");
		    fprintf(out, "mass       ");
		    fprintf(out, "density    ");
		    fprintf(out, "vrot       ");
		    fprintf(out, "omega      ");
		    fprintf(out, "sigma_r    ");
		    fprintf(out, "sigma_thet ");
		    fprintf(out, "sigma_phi  ");
		    fprintf(out, "sigma_phi(omega)\n");
		}
		fprintf(out, "%3u,%3u %6u: ", is, it, tw->n);
		fprintf(out, "%10.3e %10.3e %10.3e %10.3e %10.3e ",
			tw->r, tw->theta, centroid_z, centroid_rho, mass);
		fprintf(out, "%10.3e %10.3e %10.3e %10.3e %10.3e ",
			density, v_rot, omega, sigma_r, sigma_theta);
		fprintf(out, "%10.3e %10.3e ",
			sigma_phi, /* sig2wc, */ sigma_phiomega);

		/*
		 * fprintf(out, "cen: %10.3e %10.3e ", centroid_r,
		 * centroid_theta);
		 */
		fputc('\n', out);

	    }
	}
    }

    exit(0);



/* * * * * * * * * * * * * * * * * * * * * * * * * * */

    /* Dump d[] to "out" */

    for (i = 0; i < nd; i++) {
	fprintf(out, "%6u ", i);
	fprintf(out, "%19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e ",
		d[i].m, d[i].x, d[i].y, d[i].z, d[i].vx, d[i].vy, d[i].vz);
	fprintf(out, "%19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e %19.12e ",
		d[i].rho, d[i].r, d[i].phi, d[i].theta, d[i].csphi, d[i].snphi, d[i].cstheta, d[i].sntheta);
	fprintf(out, "%19.12e %19.12e %19.12e %19.12e",
		d[i].radius, d[i].density, d[i].vphi, d[i].vphineigh);
	fprintf(out, "%19.12e %19.12e %19.12e\n",
		d[i].sigma2neigh, d[i].Rneigh, d[i].Tneigh);
    }
    exit(0);
}

double
ke(struct data * dd)
{
    double r;

    r = dd->vx * dd->vx + dd->vy * dd->vy + dd->vz * dd->vz;
    r = r * dd->m / 2.;
    return r;
}
