/* TABLIX, PGA general timetable solver                              */
/* Copyright (C) 2002-2006 Tomaz Solc                                      */

/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU General Public License as published by    */
/* the Free Software Foundation; either version 2 of the License, or       */
/* (at your option) any later version.                                     */

/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU General Public License for more details.                            */

/* You should have received a copy of the GNU General Public License       */
/* along with this program; if not, write to the Free Software             */
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA */

/* $Id: error-pvm.c,v 1.3 2006-02-04 14:16:12 avian Exp $ */

#ifdef USE_MPI
    #include "mpi.h"
/*
#elifdef HAVE_LIBPVM3
  #include <pvm3.h>
*/
#endif
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "main.h"
#include "gettext.h"
#include "error.h"

int verbosity;
char *curmodule;

#ifdef USE_MPI
void msg_vsend(char *module, int msgtag, const char *fmt, va_list ap)
{
	int strsize, mpipack, mpisend;
	int position = 0;
	char msg[LINEBUFFSIZE];
    char buff[LINEBUFFSIZE];

	if(msgtag>verbosity&&msgtag!=MSG_FATAL) return;

	vsnprintf(msg, LINEBUFFSIZE, fmt, ap);
	msg[LINEBUFFSIZE-1]=0;

	strsize=strlen(module);
	mpipack=MPI_Pack(strsize, 1, MPI_INT, *buff, LINEBUFFSIZE, &position, MPI_COMM_WORLD);
    mpipack=MPI_Pack(module, strsize, MPI_CHAR, *buff, LINEBUFFSIZE, &position, MPI_COMM_WORLD);
	strsize=strlen(msg);
	mpipack=MPI_Pack(strsize, 1, MPI_INT, *buff, LINEBUFFSIZE, &position, MPI_COMM_WORLD);
    mpipack=MPI_Pack(msg, strsize, MPI_CHAR, *buff, LINEBUFFSIZE, &position, MPI_COMM_WORLD);

	mpisend=MPI_Send(buff, position, MPI_PACKED, parent, msgtag, MPI_COMM_WORLD);

	if(msgtag==MSG_FATAL) {
		MPI_Finalize();
		exit(1);
	}
}

#elifdef HAVE_LIBPVM3

void msg_vsend(char *module, int msgtag, const char *fmt, va_list ap)
{
	char msg[LINEBUFFSIZE];

	if(msgtag>verbosity&&msgtag!=MSG_FATAL) return;

	vsnprintf(msg, LINEBUFFSIZE, fmt, ap);
	msg[LINEBUFFSIZE-1]=0;

	pvm_initsend(0);
	pvm_pkstr(module);
	pvm_pkstr(msg);
	pvm_send(pvm_parent(), msgtag);

	if(msgtag==MSG_FATAL) {
		pvm_exit();
		exit(1);
	}
}
#endif

void msg_send(char *module, int msgtag, const char *fmt, ...)
{
	va_list ap;

	va_start(ap, fmt);

	msg_vsend(module, msgtag, fmt, ap);

	va_end(ap);
}
