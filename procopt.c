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
#include <sys/param.h>
#include <sys/utsname.h>

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

int
procopt(int *argc, char ***argv)
{
    int c;
    int i;
    int nrec;
    int r;
    char cmd[MAXPATHLEN + 100];
    char cwd[MAXPATHLEN];
    time_t t;
    FILE *wc;
    struct utsname unam;

    in = NULL;
    out = stdout;
    logfile = NULL;
    usemass = 0;
    jmax = 10;
    nproc = 4;
    nshell = 20;
    ntorus = 9;
    nrec = 0;

    while ((c = getopt(*argc, *argv, "i:j:k:l:o:s:t:hm")) != -1) {
	switch (c) {
	case 'i':
	    if ((in = fopen(optarg, "r")) == NULL)
		err(1, "input file: %s", optarg);

	    snprintf(cmd, sizeof(cmd), "/usr/bin/wc -l %s", optarg);
	    if (!(wc = popen(cmd, "r")))
		err(1, "could not popen %s", cmd);
	    r = fscanf(wc, "%d", &nrec);
	    if (r != 1 || ferror(wc))
		err(1, "trouble exectuting %s", cmd);
	    pclose(wc);
	    break;
	case 'j':
	    jmax = atoi(optarg);
	    if (jmax < 3)
		errx(1, "irrational j: %d", jmax);
	    break;
	case 'k':
	    nproc = atoi(optarg);
	    if (nproc < 1 || nproc > 1000)
		errx(1, "odd number of subprocesses (%d)", nproc);
	    break;
	case 'l':
	    if ((logfile = fopen(optarg, "w")) == NULL)
		err(1, "logfile file: %s", optarg);
	    setbuf(logfile, NULL);     /* unbuffered */
	    break;
	case 'm':
	    usemass = 1;
	    break;
	case 'o':
	    if (optarg[0] == '-') {
		out = stdout;
	    } else {
		if ((out = fopen(optarg, "w")) == NULL)
		    err(1, "out file: %s", optarg);
	    }
	    break;
	case 's':
	    nshell = atoi(optarg);
	    break;
	case 't':
	    ntorus = atoi(optarg);
	    break;
	case 'h':
	case 'H':
	    usage(0);
	    /* NOTREACHED */
	default:
	    usage(2);
	    /* NOTREACHED */
	}
    }

    if (in == NULL) {
	warnx("No input file specified.\n");
	usage(1);
    }
    if (logfile == NULL)
	if ((logfile = fopen("/dev/null", "w")) == NULL)
	    err(1, "could not open log on /dev/null");

    t = time(&t);
    if (getcwd(cwd, sizeof(cwd)) == NULL)
	err(1, "cwd");
    if (uname(&unam) == -1)
	err(1, "uname");

    fprintf(logfile, "# Command line: ");
    for (i = 0; i < *argc; i++)
	fprintf(logfile, "%s ", (*argv)[i]);
    fputc('\n', logfile);
    fprintf(logfile, "# Run commenced at %s", ctime(&t));
    fprintf(logfile, "# Current working directory: %s\n", cwd);
    fprintf(logfile, "# Run on %s : %s %s %s %s\n",
     unam.nodename, unam.sysname, unam.release, unam.machine, unam.version);

    fprintf(out, "# Command line: ");
    for (i = 0; i < *argc; i++)
	fprintf(out, "%s ", (*argv)[i]);
    fputc('\n', out);
    fprintf(out, "# Run commenced at %s", ctime(&t));
    fprintf(out, "# Current working directory: %s\n", cwd);
    fprintf(out, "# Run on %s : %s %s %s %s\n",
     unam.nodename, unam.sysname, unam.release, unam.machine, unam.version);

    *argc -= optind;
    *argv += optind;
    return nrec;
}
