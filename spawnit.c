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

#undef VERBOSE

int
 spawnit(int nfork, int nitem, int kidf(int ikid, int istart, int icount, va_list),...) {
/*
 * fork()s n processes.  Each child calls kidf.
 * the first three arguments to kidf are ints to be determined here.
 *
 *   ikid: an integer [0,nfork)  Unique to each instance of kidf.
 *   istart: begin here, calculated, see below
 *   icount: do this many
 *   vargs:  do it to these things.  The variable list passed to
 *           spawnit is passed on the kidf.
 *   What kidf does to the stuff in vargs, how it interprets istart and
 *   icount, is a matter between kidf and the caller of spawnint.
 *
 * this function waits for all the fork()ed children to terminate.
 *
 * an example usage:
 *
 * a calling routine wishes to fork() the work for multiplying two
 * arrays,  i.e., to calculate z(i) = x(i)*y(i);  It wishes to
 * partition this work among 4 children, each to handle a fourth
 * of the array.
 *
 * in the caller:
 * 	double x[100], y[100], z[100];
 *      ...
 *      spawnit( 4, 100, arrmult, x, y, z);
 *
 * the kidf routine:
 *
 * int arrmult(int ikid, int i0, int count, va_list args)
 * {
 *      double *x,*y,*z;
 *      int i;
 *
 *      x = va_arg(args, double *);
 *      y = va_arg(args, double *);
 *      z = va_arg(args, double *);
 *      for(i=i0; i< i0+count;i++)
 *          z[i]=x[i]*y[i];
 *      return;
 * }
 */


    pid_t *kids;
    pid_t pid;
    int status;
    int i;
    int r;
    int istart, icount, imod;

    va_list kidargs;

    va_start(kidargs, kidf);

    if (!(kids = calloc(nfork, sizeof(pid_t))))
	err(1, "kids");

    icount = nitem / nfork;
    imod = nitem % nfork;
    istart = -icount;

    for (i = 0; i < nfork; i++) {

	istart += icount;
	if (i == nfork - 1)
	    icount += imod;
refork:
	fflush(NULL);
	if ((pid = fork()) < 0) {
	    if (errno == EAGAIN) {
		errno = 0;
		pid = wait(&status);
		if (pid == -1)
		    err(1, "bizarre fork/wait, line %d", __LINE__);
		if (!(WIFEXITED(status) || WIFSIGNALED(status)))
		    goto refork;
#ifdef VERBOSE
		fprintf(stderr, "Child %d complete, status 0x%X\n", pid, status);
#endif
		goto refork;
	    } else
		err(1, "fork at %d, errno %d", i, errno);
	} else if (pid) {
	    kids[i] = pid;	       /* is acually a do-nought at present */
	} else {
	    /* Execution in one of the children: */
	    r = kidf(i, istart, icount, kidargs);
	    exit(r);
	}
#ifdef VERBOSE
	fprintf(stderr, "forked: %d\n", pid);
#endif
    }
    for (i = 0; i < nfork; i++) {      /* the do limits avoid an infinite
				        * loop */
#ifdef VERBOSE
	fprintf(stderr, "Waiting\n");
#endif
	pid = wait(&status);
	if (pid == -1) {
	    if (errno == ECHILD)
		break;
	    else
		warn("collected -1, %d", errno);
	}
#ifdef VERBOSE
	fprintf(stderr, "Collected %d, status 0x%X\n", pid, status);
#endif
	/* should probably check status with WIFxxx() macros. */
    }
    free(kids);
    return 0;
}
