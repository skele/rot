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

#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>

#include "rot.h"

extern char *__progname;

void
usage(int arg)
{
    fprintf(stderr,
	    "usage: %s -i file [-o file] [-l file] [-j int] [-s int] [-t int] [-hm]\n",
	    __progname);
    fprintf(stderr,
      "\t-i file\tTake input from <file>. No default, mandatory option.\n");
    fprintf(stderr,
	  "\t-o file\tWrite output on <file>. Default, standard output.\n");
    fprintf(stderr,
	    "\t-l file\tWrite log info on <file>. Default, /dev/null.\n");
    fprintf(stderr,
	    "\t-j jneigh\tUse <jneigh> nearest neighbors. Default, 10.\n");
    fprintf(stderr,
	    "\t-k nproc\tNumber of subprocesses to spawn. Default, 4.\n");
    fprintf(stderr,
      "\t\tUseful range: 1 to a small multiple of available processors.\n");
    fprintf(stderr, "\t-h\tPrint this help and exit.\n");
    fprintf(stderr,
    "\t-m\tUse mass in calculating density. By default, do not use mass.\n");
    fprintf(stderr,
	    "\t-s nshell\tNumber of spherical shells. Default, 20.\n");
    fprintf(stderr,
	  "\t-t ntorus\tNumber of tori per spherical shell. Default, 9.\n");
    exit(arg);
}
