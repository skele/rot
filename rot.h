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

/* ISO 31-11 Notation */

#include <stdio.h>
#include <stdarg.h>

#define NPROCMAX	128
#define ND 60000
#define SQD(x,y) (((x)-(y))*((x)-(y)))


#ifdef DOCLOCK
#define clockit() clocker(__FILE__, __LINE__)
#else
#define clockit()
#endif


struct qtype {
    int idx;
    double dsq;
};

struct data {
    int i;
    double m;
    double x;
    double y;
    double z;
    double vx;
    double vy;
    double vz;
    double rho;			       /* cylindrical distance */
    double r;			       /* spherical coordinates */
    double theta;
    double phi;
    double snphi;
    double csphi;
    double sntheta;
    double cstheta;
    double radius;
    double density;
    double vr;
    double vrho;
    double vtheta;
    double vphi;
    double vphineigh;
    double sigma2neigh;
    double Rneigh;
    double Tneigh;
    struct qtype *nn;		       /* nn[JMAX]; */
};


struct shell_t {
    int n;
    int ib;
    double r;
    double dr;
    double theta;
    double dtheta;
    double bound;
    double volume;
    struct shell_t *sub;
};

/* Globals */

FILE *in;
FILE *out;
FILE *logfile;

int usemass;
int nproc;
int jmax;
int nshell;
int ntorus;


/*************/


/* misc.c */
extern void xlate(struct data * d, int nd, double x0, double y0, double z0);
extern void derive(struct data * d, int nd);
extern void cross(double *c1, double *c2, double *c3, double a1, double a2, double a3, double b1, double b2, double b3, double f);
extern double dot(double a1, double a2, double a3, double b1, double b2, double b3);
extern double powi(double x, int i);
extern void initclock(void);
extern void clocker(char *file, int line);
extern int compr(const void *a, const void *b);
extern int comptheta(const void *a, const void *b);
extern int compq(const void *a, const void *b);
extern int binary_search(struct qtype U[], int low, int high, double target);
extern int fpshellhdr(FILE * f);
extern int fpshell(FILE * f, const char *title, struct shell_t * s);

/* procopt.c */
extern int procopt(int *argc, char ***argv);

/* readdbl.c */
extern struct data *readdbl(FILE * in, int *np, int nskip, int nread);
int readdblm(FILE * in, struct data * d, int nd, int nskip, int nread);

/* spawnit.c */
extern int spawnit(int nfork, int nitem, int kidf(int ikid, int istart, int icount, va_list),...);

/* usage.c */
extern void usage(int arg);

/* mkshell.c */
struct shell_t *mkshell(struct data * d, int nd, int nshell);
struct shell_t *mktorus(struct data * d, int nd, struct shell_t s, int nt);

/*distances.cu*/
extern void obtain_densities(float *mx, float *my, float *mz, float *mradius, float *mdensity, int istart, int icount, int ikid, int N, float mass);
