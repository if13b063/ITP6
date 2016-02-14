/* TABLIX, PGA general timetable solver                              */
/* Copyright (C) 2002-2004 Tomaz Solc                                      */

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

/* $Id: transfer.c,v 1.6 2006-02-04 14:16:12 avian Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef USE_MPI
#include "mpi.h"

#include "main.h"
#include "error.h"
#include "assert.h"
#include "chromo.h"
#include "transfer.h"

/** @file
 * @brief File transfer through PVM.  */

/** @brief Wait for file sent from sender.
 *
 * @param filename File will be saved under this name.
 * @param sender Task ID of the sender task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int file_recv(char *filename, int sender, int msgtag)
{
	FILE *handle;
	char *buff;
	int notend;
	int len;
  MPI_Status status;
  MPI_Comm comm, comm2, comm3;

	buff=malloc(LINEBUFFSIZE);
	if (buff==NULL) return(-1);

	handle=fopen(filename, "w");
	if (handle==NULL) {
		free(buff);
		return(-1);
	}

	do {
    result= MPI_Recv(&notend, 1, MPI_INT, sender, msgtag, comm, &status);
			if(result!=MPI_SUCCESS){
        free(buff);
        fclose(handle);
        return(-1);
      }
    result= MPI_Recv(&len, 1, MPI_INT, sender, msgtag, comm2, &status);
			if(result!=MPI_SUCCESS){
        free(buff);
        fclose(handle);
        return(-1);
      }
		result= MPI_Recv(buff, len, MPI_CHAR, sender, msgtag, comm3, &status);
			if(result!=MPI_SUCCESS){
        free(buff);
        fclose(handle);
        return(-1);
      }

		fprintf(handle, "%s", buff);
	} while (notend);

	fclose(handle);
	free(buff);
	return(0);
}

/** @brief Multicast file to all recipients.
 *
 * @param filename File to multicast.
 * @param recipients Array of \a num task IDs.
 * @param num Number of recipients.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int file_mcast(char *filename, int *recipients, int num, int msgtag)
{
	FILE *handle;
	char *buff;
	int notend;
	int len;
	int result;
	MPI_Comm comm, comm2, comm3, comm4, comm5, comm6;

        buff=malloc(LINEBUFFSIZE);
        if (buff==NULL) return(-1);

	handle=fopen(filename, "r");
	if (handle==NULL) {
		free(buff);
		return(-1);
	}

        notend=1;
        while (fgets(buff, LINEBUFFSIZE, handle)!=NULL) {
                result=MPI_Bcast(&notend, 1, MPI_INT, 0, msgtag, comm);
                if(result!=MPI_SUCCESS){
                  fclose(handle);
                  free(buff);
                  return(-1);
                }
                len=strlen(buff);
                result=MPI_Bcast(&len, 1, MPI_INT, 0, msgtag, comm2);
                if(result!=MPI_SUCCESS){
                  fclose(handle);
                  free(buff);
                  return(-1);
                }
                result=MPI_Bcast(buff, len, MPI_CHAR, 0, msgtag, comm3);
                if(result!=MPI_SUCCESS){
                  fclose(handle);
                  free(buff);
                  return(-1);
                }
        }

        notend=0;
        strcpy(buff, "<!-- End of file -->\n");
        result=MPI_Bcast(&notend, 1, MPI_INT, 0, msgtag, comm4);
        if(result!=MPI_SUCCESS){
          fclose(handle);
          free(buff);
          return(-1);
        }
        len=strlen(buff);
        result=MPI_Bcast(&len, 1, MPI_INT, 0, msgtag, comm5);
        if(result!=MPI_SUCCESS){
          fclose(handle);
          free(buff);
          return(-1);
        }
        result=MPI_Bcast(buff, len, MPI_CHAR, 0, msgtag, comm6);
        if(result!=MPI_SUCCESS){
          fclose(handle);
          free(buff);
          return(-1);
        }

        fclose(handle);
	free(buff);
	return(0);
}

/** @brief Send file to one recipient.
 *
 * @param filename File to multicast.
 * @param recipient Task ID of the recipient task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int file_send(char *filename, int recipient, int msgtag)
{
	FILE *handle;
	char *buff;
	int mpisend, mpisend2, mpisend3, mpisend4, mpisend5, mpisend6;
  MPI_Comm comm, comm2, comm3, comm4, comm5, comm6;
	int notend;
	int len;

	//MPI_Init(&argc, &argv);?

        buff=malloc(LINEBUFFSIZE);
        if (buff==NULL) {
		return(-1);
	}

	handle=fopen(filename, "r");
	if (handle==NULL) {
		free(buff);
		return(-1);
	}

        notend=1;
        while (fgets(buff, LINEBUFFSIZE, handle)!=NULL) {
                mpisend=MPI_Send(&notend, 1, MPI_INT, recipient, msgtag, comm);
                if(mpisend!=MPI_SUCCESS){
                  fclose(handle);
                  free(buff);
                  return(-1);
                }
                len=strlen(buff);
                mpisend2=MPI_Send(&len, 1, MPI_INT, recipient, msgtag, comm2);
                if(mpisend2!=MPI_SUCCESS){
                  fclose(handle);
                  free(buff);
                  return(-1);
                }
                mpisend3=MPI_Send(buff, len, MPI_CHAR, recipient, msgtag, comm3);
                if(mpisend3!=MPI_SUCCESS){
                  fclose(handle);
                  free(buff);
                  return(-1);
                }
        }

        notend=0;
        strcpy(buff, "<!-- End of file -->\n");
        mpisend4=MPI_Send(&notend, 1, MPI_INT, recipient, msgtag, comm4);
        if(mpisend4!=MPI_SUCCESS){
          fclose(handle);
          free(buff);
          return(-1);
        }
        len=strlen(buff);
        mpisend5=MPI_Send(&len, 1, MPI_INT, recipient, msgtag, comm5);
        if(mpisend5!=MPI_SUCCESS){
          fclose(handle);
          free(buff);
          return(-1);
        }
        mpisend6=MPI_Send(buff, len, MPI_CHAR, recipient, msgtag, comm6);
        if(mpisend6!=MPI_SUCCESS){
          fclose(handle);
          free(buff);
          return(-1);
        }

        fclose(handle);
	free(buff);
	return(0);
}

/** @brief Wait for timetables from sender.
 *
 * Received timetables are stored at the end of the population.
 *
 * @param pop Pointer to the population struct.
 * @param num Number of chromosomes to receive.
 * @param sender Task ID of the sender task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int table_recv(population *pop, int num, int sender, int msgtag)
{
	int n,m;
	int result;
	int size;
	MPI_Status status
	MPI_Comm comm;
	table *tab;

	assert(num>0);
	assert(num<=pop->size);

	for(n=pop->size-num;n<pop->size;n++) {
		tab=pop->tables[n];
		for(m=0;m<tab->typenum;m++) {
			size=tab->chr[m].gennum;
			result= MPI_Recv(tab->chr[m].gen, 1, MPI_PACKED, sender, msgtag, comm, &status);
			if(result!=MPI_SUCCESS) return(-1);
		}

		tab->fitness=-1;
	}

	return(0);
}

/** @brief Send timetables to a recipient.
 *
 * Timetables are read from the start of the population.
 *
 * @param pop Pointer to the population struct.
 * @param num Number of chromosomes to send.
 * @param recipient Task ID of the recipient task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int table_send(population *pop, int num, int recipient, int msgtag)
{
	int n,m;
	int result;
	int size;
	table *tab;
	MPI_Comm comm;

	assert(num>0);
	assert(num<=pop->size);

	for(n=0;n<num;n++) {
		tab=pop->tables[n];
		for(m=0;m<tab->typenum;m++) {
			size=tab->chr[m].gennum;
      result = MPI_Send(tab->chr[m].gen, 1, MPI_PACKED, recipient, msgtag, comm);
			if(result!=MPI_SUCCESS) return(-1);
		}
	}

	return(0);
}

/** @brief Receive the entire population.
 *
 * @param sender Task ID of the sender task.
 * @param msgtag Message tag.
 * @return Pointer to the allocated population struct or NULL on error. */
population *population_recv(int sender, int msgtag)
{
	int result;
	int size, gencnt;
	int typenum;
	int tuplenum;
  MPI_Status status;
  MPI_Comm comm;
  int recv[4];

	population *pop;

  result=MPI_Recv(&recv, 4, MPI_INT, sender, msgtag, comm, &status);
  if(result!=MPI_SUCCESS) return(-1);
	size=recv[0];
	gencnt=recv[1];
	typenum=recv[2];
	tuplenum=recv[3];

	pop=population_new(size, typenum, tuplenum);
	if(pop==NULL) return(NULL);

	pop->gencnt=gencnt;

	result=table_recv(pop, pop->size, sender, msgtag);
	if(result<0) {
		population_free(pop);
		return(NULL);
	}

	return(pop);
}

/** @brief Send the entire population.
 *
 * @param pop Pointer to the population struct.
 * @param recipient Task ID of the recipient task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int population_send(population *pop, int recipient, int msgtag)
{
	int result;
  int send[4];
  MPI_Comm comm;

	assert(pop!=NULL);
	assert(pop->size>0);
	assert(pop->tables[0]->typenum>0);
	assert(pop->tables[0]->chr[0].gennum>0);

//	result=pvm_pkint(&pop->size, 1, 1);
//	if(result<0) return(-1);
//	result=pvm_pkint(&pop->gencnt, 1, 1);
//	if(result<0) return(-1);
//	result=pvm_pkint(&pop->tables[0]->typenum, 1, 1);
//	if(result<0) return(-1);
//	result=pvm_pkint(&pop->tables[0]->chr[0].gennum, 1, 1);
//	if(result<0) return(-1);

  send[0]=pop->size;
  send[1]=pop->dencnt;
  send[2]=pop->tables[0]->typenum;
  send[3]=pop->tables[0]->chr[0].gennum;
  mpisend=MPI_Send(&send, 4, MPI_INT, recipient, msgtag, comm);
  if(mpisend!=MPI_SUCCESS) return(-1);

	result=table_send(pop, pop->size, recipient, msgtag);
	if(result<0) return(-1);

	return(0);
}

<<<<<<< HEAD
//#elifdef HAVE_LIBPVM3
//#include <pvm3.h>
//
//#include "main.h"
//#include "error.h"
//#include "assert.h"
//#include "chromo.h"
//#include "transfer.h"
//
///** @file
// * @brief File transfer through PVM.  */
//
///** @brief Wait for file sent from sender.
// *
// * @param filename File will be saved under this name.
// * @param sender Task ID of the sender task.
// * @param msgtag Message tag.
// * @return 0 on success and -1 on error. */
//int file_recv(char *filename, int sender, int msgtag)
//{
//	FILE *handle;
//	char *buff;
//	int notend;
//
//	buff=malloc(LINEBUFFSIZE);
//	if (buff==NULL) return(-1);
//
//	handle=fopen(filename, "w");
//	if (handle==NULL) {
//		free(buff);
//		return(-1);
//	}
//
//	do {
//		pvm_recv(sender, msgtag);
//
//		pvm_upkint(&notend, 1, 1);
//		pvm_upkstr(buff);
//
//		fprintf(handle, "%s", buff);
//	} while (notend);
//
//	fclose(handle);
//	free(buff);
//	return(0);
//}
//
///** @brief Multicast file to all recipients.
// *
// * @param filename File to multicast.
// * @param recipients Array of \a num task IDs.
// * @param num Number of recipients.
// * @param msgtag Message tag.
// * @return 0 on success and -1 on error. */
//int file_mcast(char *filename, int *recipients, int num, int msgtag)
//{
//	FILE *handle;
//	char *buff;
//
//	int notend;
//
//        buff=malloc(LINEBUFFSIZE);
//        if (buff==NULL) return(-1);
//
//	handle=fopen(filename, "r");
//	if (handle==NULL) {
//		free(buff);
//		return(-1);
//	}
//
//        notend=1;
//        while (fgets(buff, LINEBUFFSIZE, handle)!=NULL) {
//                pvm_initsend(0);
//                pvm_pkint(&notend, 1, 1);
//                pvm_pkstr(buff);
//                pvm_mcast(recipients, num, msgtag);
//        }
//
//        pvm_initsend(0);
//        notend=0;
//        strcpy(buff, "<!-- End of file -->\n");
//        pvm_pkint(&notend, 1, 1);
//        pvm_pkstr(buff);
//        pvm_mcast(recipients, num, msgtag);
//
//        fclose(handle);
//	free(buff);
//	return(0);
//}
//
///** @brief Send file to one recipient.
// *
// * @param filename File to multicast.
// * @param recipient Task ID of the recipient task.
// * @param msgtag Message tag.
// * @return 0 on success and -1 on error. */
//int file_send(char *filename, int recipient, int msgtag)
//{
//	FILE *handle;
//	char *buff;
//
//	int notend;
//
//        buff=malloc(LINEBUFFSIZE);
//        if (buff==NULL) {
//		return(-1);
//	}
//
//	handle=fopen(filename, "r");
//	if (handle==NULL) {
//		free(buff);
//		return(-1);
//	}
//
//        notend=1;
//        while (fgets(buff, LINEBUFFSIZE, handle)!=NULL) {
//                pvm_initsend(0);
//                pvm_pkint(&notend, 1, 1);
//                pvm_pkstr(buff);
//                pvm_send(recipient, msgtag);
//        }
//
//        pvm_initsend(0);
//        notend=0;
//        strcpy(buff, "<!-- End of file -->\n");
//        pvm_pkint(&notend, 1, 1);
//        pvm_pkstr(buff);
//        pvm_send(recipient, msgtag);
//
//        fclose(handle);
//	free(buff);
//	return(0);
//}
//
///** @brief Wait for timetables from sender.
// *
// * Received timetables are stored at the end of the population.
// *
// * @param pop Pointer to the population struct.
// * @param num Number of chromosomes to receive.
// * @param sender Task ID of the sender task.
// * @param msgtag Message tag.
// * @return 0 on success and -1 on error. */
//int table_recv(population *pop, int num, int sender, int msgtag)
//{
//	int n,m;
//	int result;
//	int size;
//	table *tab;
//
//	assert(num>0);
//	assert(num<=pop->size);
//
//	pvm_recv(sender, msgtag);
//
//	for(n=pop->size-num;n<pop->size;n++) {
//		tab=pop->tables[n];
//		for(m=0;m<tab->typenum;m++) {
//			size=tab->chr[m].gennum;
//			result=pvm_upkint(tab->chr[m].gen, size, 1);
//
//			if(result<0) return(-1);
//		}
//
//		tab->fitness=-1;
//	}
//
//	return(0);
//}
//
///** @brief Send timetables to a recipient.
// *
// * Timetables are read from the start of the population.
// *
// * @param pop Pointer to the population struct.
// * @param num Number of chromosomes to send.
// * @param recipient Task ID of the recipient task.
// * @param msgtag Message tag.
// * @return 0 on success and -1 on error. */
//int table_send(population *pop, int num, int recipient, int msgtag)
//{
//	int n,m;
//	int result;
//	int size;
//	table *tab;
//
//	assert(num>0);
//	assert(num<=pop->size);
//
//	pvm_initsend(0);
//
//	for(n=0;n<num;n++) {
//		tab=pop->tables[n];
//		for(m=0;m<tab->typenum;m++) {
//			size=tab->chr[m].gennum;
//			result=pvm_pkint(tab->chr[m].gen, size, 1);
//
//			if(result<0) return(-1);
//		}
//	}
//
//	result=pvm_send(recipient, msgtag);
//	if(result<0) return(-1);
//
//	return(0);
//}
//
///** @brief Receive the entire population.
// *
// * @param sender Task ID of the sender task.
// * @param msgtag Message tag.
// * @return Pointer to the allocated population struct or NULL on error. */
//population *population_recv(int sender, int msgtag)
//{
//	int result;
//	int size, gencnt;
//	int typenum;
//	int tuplenum;
//
//	population *pop;
//
//	pvm_recv(sender, msgtag);
//
//	result=pvm_upkint(&size, 1, 1);
//	if(result<0) return(NULL);
//	result=pvm_upkint(&gencnt, 1, 1);
//	if(result<0) return(NULL);
//	result=pvm_upkint(&typenum, 1, 1);
//	if(result<0) return(NULL);
//	result=pvm_upkint(&tuplenum, 1, 1);
//	if(result<0) return(NULL);
//
//	pop=population_new(size, typenum, tuplenum);
//	if(pop==NULL) return(NULL);
//
//	pop->gencnt=gencnt;
//
//	result=table_recv(pop, pop->size, sender, msgtag);
//	if(result<0) {
//		population_free(pop);
//		return(NULL);
//	}
//
//	return(pop);
//}
//
///** @brief Send the entire population.
// *
// * @param pop Pointer to the population struct.
// * @param recipient Task ID of the recipient task.
// * @param msgtag Message tag.
// * @return 0 on success and -1 on error. */
//int population_send(population *pop, int recipient, int msgtag)
//{
//	int result;
//
//	assert(pop!=NULL);
//	assert(pop->size>0);
//	assert(pop->tables[0]->typenum>0);
//	assert(pop->tables[0]->chr[0].gennum>0);
//
//	pvm_initsend(0);
//
//	result=pvm_pkint(&pop->size, 1, 1);
//	if(result<0) return(-1);
//	result=pvm_pkint(&pop->gencnt, 1, 1);
//	if(result<0) return(-1);
//	result=pvm_pkint(&pop->tables[0]->typenum, 1, 1);
//	if(result<0) return(-1);
//	result=pvm_pkint(&pop->tables[0]->chr[0].gennum, 1, 1);
//	if(result<0) return(-1);
//
//	result=pvm_send(recipient, msgtag);
//	if(result<0) return(-1);
//
//	result=table_send(pop, pop->size, recipient, msgtag);
//	if(result<0) return(-1);
//
//	return(0);
//}
=======
/*
#elifdef HAVE_LIBPVM3
#include <pvm3.h>

#include "main.h"
#include "error.h"
#include "assert.h"
#include "chromo.h"
#include "transfer.h"

/** @file
 * @brief File transfer through PVM.  */

/** @brief Wait for file sent from sender.
 *
 * @param filename File will be saved under this name.
 * @param sender Task ID of the sender task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int file_recv(char *filename, int sender, int msgtag)
{
	FILE *handle;
	char *buff;
	int notend;

	buff=malloc(LINEBUFFSIZE);
	if (buff==NULL) return(-1);

	handle=fopen(filename, "w");
	if (handle==NULL) {
		free(buff);
		return(-1);
	}

	do {
		pvm_recv(sender, msgtag);

		pvm_upkint(&notend, 1, 1);
		pvm_upkstr(buff);

		fprintf(handle, "%s", buff);
	} while (notend);

	fclose(handle);
	free(buff);
	return(0);
}

/** @brief Multicast file to all recipients.
 *
 * @param filename File to multicast.
 * @param recipients Array of \a num task IDs.
 * @param num Number of recipients.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int file_mcast(char *filename, int *recipients, int num, int msgtag)
{
	FILE *handle;
	char *buff;

	int notend;

        buff=malloc(LINEBUFFSIZE);
        if (buff==NULL) return(-1);

	handle=fopen(filename, "r");
	if (handle==NULL) {
		free(buff);
		return(-1);
	}

        notend=1;
        while (fgets(buff, LINEBUFFSIZE, handle)!=NULL) {
                pvm_initsend(0);
                pvm_pkint(&notend, 1, 1);
                pvm_pkstr(buff);
                pvm_mcast(recipients, num, msgtag);
        }

        pvm_initsend(0);
        notend=0;
        strcpy(buff, "<!-- End of file -->\n");
        pvm_pkint(&notend, 1, 1);
        pvm_pkstr(buff);
        pvm_mcast(recipients, num, msgtag);

        fclose(handle);
	free(buff);
	return(0);
}

/** @brief Send file to one recipient.
 *
 * @param filename File to multicast.
 * @param recipient Task ID of the recipient task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int file_send(char *filename, int recipient, int msgtag)
{
	FILE *handle;
	char *buff;

	int notend;

        buff=malloc(LINEBUFFSIZE);
        if (buff==NULL) {
		return(-1);
	}

	handle=fopen(filename, "r");
	if (handle==NULL) {
		free(buff);
		return(-1);
	}

        notend=1;
        while (fgets(buff, LINEBUFFSIZE, handle)!=NULL) {
                pvm_initsend(0);
                pvm_pkint(&notend, 1, 1);
                pvm_pkstr(buff);
                pvm_send(recipient, msgtag);
        }

        pvm_initsend(0);
        notend=0;
        strcpy(buff, "<!-- End of file -->\n");
        pvm_pkint(&notend, 1, 1);
        pvm_pkstr(buff);
        pvm_send(recipient, msgtag);

        fclose(handle);
	free(buff);
	return(0);
}

/** @brief Wait for timetables from sender.
 *
 * Received timetables are stored at the end of the population.
 *
 * @param pop Pointer to the population struct.
 * @param num Number of chromosomes to receive.
 * @param sender Task ID of the sender task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int table_recv(population *pop, int num, int sender, int msgtag)
{
	int n,m;
	int result;
	int size;
	table *tab;

	assert(num>0);
	assert(num<=pop->size);

	pvm_recv(sender, msgtag);

	for(n=pop->size-num;n<pop->size;n++) {
		tab=pop->tables[n];
		for(m=0;m<tab->typenum;m++) {
			size=tab->chr[m].gennum;
			result=pvm_upkint(tab->chr[m].gen, size, 1);

			if(result<0) return(-1);
		}

		tab->fitness=-1;
	}

	return(0);
}

/** @brief Send timetables to a recipient.
 *
 * Timetables are read from the start of the population.
 *
 * @param pop Pointer to the population struct.
 * @param num Number of chromosomes to send.
 * @param recipient Task ID of the recipient task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int table_send(population *pop, int num, int recipient, int msgtag)
{
	int n,m;
	int result;
	int size;
	table *tab;

	assert(num>0);
	assert(num<=pop->size);

	pvm_initsend(0);

	for(n=0;n<num;n++) {
		tab=pop->tables[n];
		for(m=0;m<tab->typenum;m++) {
			size=tab->chr[m].gennum;
			result=pvm_pkint(tab->chr[m].gen, size, 1);

			if(result<0) return(-1);
		}
	}

	result=pvm_send(recipient, msgtag);
	if(result<0) return(-1);

	return(0);
}

/** @brief Receive the entire population.
 *
 * @param sender Task ID of the sender task.
 * @param msgtag Message tag.
 * @return Pointer to the allocated population struct or NULL on error. */
population *population_recv(int sender, int msgtag)
{
	int result;
	int size, gencnt;
	int typenum;
	int tuplenum;

	population *pop;

	pvm_recv(sender, msgtag);

	result=pvm_upkint(&size, 1, 1);
	if(result<0) return(NULL);
	result=pvm_upkint(&gencnt, 1, 1);
	if(result<0) return(NULL);
	result=pvm_upkint(&typenum, 1, 1);
	if(result<0) return(NULL);
	result=pvm_upkint(&tuplenum, 1, 1);
	if(result<0) return(NULL);

	pop=population_new(size, typenum, tuplenum);
	if(pop==NULL) return(NULL);

	pop->gencnt=gencnt;

	result=table_recv(pop, pop->size, sender, msgtag);
	if(result<0) {
		population_free(pop);
		return(NULL);
	}

	return(pop);
}

/** @brief Send the entire population.
 *
 * @param pop Pointer to the population struct.
 * @param recipient Task ID of the recipient task.
 * @param msgtag Message tag.
 * @return 0 on success and -1 on error. */
int population_send(population *pop, int recipient, int msgtag)
{
	int result;

	assert(pop!=NULL);
	assert(pop->size>0);
	assert(pop->tables[0]->typenum>0);
	assert(pop->tables[0]->chr[0].gennum>0);

	pvm_initsend(0);

	result=pvm_pkint(&pop->size, 1, 1);
	if(result<0) return(-1);
	result=pvm_pkint(&pop->gencnt, 1, 1);
	if(result<0) return(-1);
	result=pvm_pkint(&pop->tables[0]->typenum, 1, 1);
	if(result<0) return(-1);
	result=pvm_pkint(&pop->tables[0]->chr[0].gennum, 1, 1);
	if(result<0) return(-1);

	result=pvm_send(recipient, msgtag);
	if(result<0) return(-1);

	result=table_send(pop, pop->size, recipient, msgtag);
	if(result<0) return(-1);

	return(0);
}
*/
>>>>>>> ba75a5cb2033ab1db145fdbb71a4c6904b5d0e72
#endif
