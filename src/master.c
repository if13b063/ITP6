/* TABLIX, PGA general timetable solver                                    */
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

/* $Id: master.c,v 1.58 2006-08-29 14:32:24 avian Exp $ */

/* This whole thing needs to be cleaned up and documented. Any volunteers? */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_LOCALE_H
#include <locale.h>
#endif

#ifndef HAVE_LIBPVM3
    #include "mpi.h"
    #include <signal.h>
/*
#elifdef HAVE_LIBPVM3
  #include <pvm3.h>
  #include <signal.h>
*/
#endif


#include "chromo.h"
#include "main.h"
#include "counter.h"
#include "transfer.h"
#include "gettext.h"
#include "error.h"
#include "assert.h"

#ifndef HAVE_LIBPVM3
#include "nodes.h"

//#elifdef HAVE_LIBPVM3
//  #include "nodes.h"
#endif
int count=0, flag, msend2, msend;

char *cmd;
char prefix[256];

#ifndef HAVE_LIBPVM3
int ctrlc;
int timeout_reached;

int numlocals;
int maxlocals;
struct timeval t;
int mpimaker;

//#elifdef HAVE_LIBPVM3
//  int ctrlc;
//  int timeout_reached;
//
//  int numlocals;
//  int maxlocals;
//  struct timeval t;
#endif

void sighandler(int num);

void print_copyright()
{
        printf("\n");
        printf(_("\
This program is free software; you can redistribute it and/or modify\n\
it under the terms of the GNU General Public License as published by\n\
the Free Software Foundation; either version 2 of the License, or\n\
(at your option) any later version.\n"));
        printf("\n");
        printf(_("\
This program is distributed in the hope that it will be useful,\n\
but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
GNU General Public License for more details.\n"));
        printf("\n");
        printf(_("\
You should have received a copy of the GNU General Public License\n\
along with this program; if not, write to the Free Software\n\
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA\n"));
        printf("\n");

        printf("$Id: master.c,v 1.58 2006-08-29 14:32:24 avian Exp $\n");
        printf(_("Compile time options:\n"));

        #ifndef HAVE_LIBPVM3
        printf(_("\tParallel genetic algorithm (using MPI)\n"));

//        #elifdef HAVE_LIBPVM3
//          printf(_("\tParallel genetic algorithm (using PVM3)\n"));
        #else
        printf(_("\tLinear genetic algorith (debug mode)\n"));
        #endif
        #ifdef HAVE_CONV
        printf(_("\tSaving convergence info\n"));
        #endif
	printf(_("\tSearching for modules in %s\n"), HAVE_MODULE_PATH);
	printf(_("\tModule documentation available in %s\n"), HAVE_DATA_PATH);
	#ifdef DEBUG
	printf(_("\tAssertion checks enabled (debug mode)\n"));
	#endif
	#if ENABLE_NLS
	printf(_("\tNational language support enabled\n"));
	#endif
        exit(0);
}

void print_syntax()
{
        printf(_("Usage: %s [OPTION]... [FILE]\n"), cmd);
        printf("\n");

        #ifndef HAVE_LIBPVM3
        printf(_("\
  -n NODES              set number of computing nodes (default 4)\n"));

//        #elifdef HAVE_LIBPVM3
//          printf(_("\
//    -n NODES              set number of computing nodes (default 4)\n"));
        #endif

        printf(_("\
  -o PREFIX             use PREFIX when naming output files\n\
  -r                    restore saved populations\n\
  -h                    display this help\n\
  -v                    display version, compile time options and\n\
                        copyright information\n"));

  #ifndef HAVE_LIBPVM3
	printf(_("\
  -l NODES              set number of computing nodes that are allowed to\n\
                        simultaneously perform local search (set to -1 for \n\
                        no limit or 0 to disable local search. Default 1)\n"));

//	#elifdef HAVE_LIBPVM3
//    printf(_("\
//    -l NODES              set number of computing nodes that are allowed to\n\
//                          simultaneously perform local search (set to -1 for \n\
//                          no limit or 0 to disable local search. Default 1)\n"));
	#endif

	printf(_("\
  -d LEVEL              verbosity level (default 2)\n\
  -t N                  stop if no solution is found after N minutes (set\n\
                        to 0 for no time limit. Default 0)\n\
  -p PARAMETERS         set algorithm parameters. PARAMETERS is a comma\n\
                        separated string of options.\n"));
	printf(_("\
  -i PATH		set path to fitness modules\n\
                        (Default %s)\n"),
  			HAVE_MODULE_PATH);

        printf("\n");
        printf(_("Please report bugs to <tomaz.solc@tablix.org>.\n"));
        exit(0);
}



int main(int argc, char *argv[])
{
	char *loc=NULL;

    int c;
    
    

    #ifndef HAVE_LIBPVM3
    int mpimaker = MPI_Init(&argc, &argv);
    int d, a, msgtag, sender, sendernum;
	int n;

    double sum;
    int running;
	int nodereq;
    int restore;

    int *subtotals;
	int gnum;

	struct timeval start,end;
    char *buff;
	char *module;
    char fn[256];

	int timeout;

	population *pop=NULL;

    
    int strsize, mpibcast, mpibcast2;
    
    int recvi, recvi2, recvi3, recvi4, recvi5, recvi6, recvi7, recvi8, info;
    int unpack1, unpack2, unpack3, unpack4, unpack5, unpack6, unpack7, unpack8, unpack9, unpack10, unpack11, unpack12, unpack13;
    int *position;
    int intarr[3];
	MPI_Status st1, st2;

    #ifdef HAVE_CONV
        FILE **convfiles;
    #endif

//        #elifdef HAVE_LIBPVM3
//        int d, a, msgtag, sender, sendernum;
//	int n;
//
//        double sum;
//        int running;
//	int nodereq;
//        int restore;
//
//        int *subtotals;
//	int gnum;
//
//	struct timeval start,end;
//        char *buff;
//	char *module;
//        char fn[256];
//
//	int timeout;
//
//	population *pop=NULL;
//
//        #ifdef HAVE_CONV
//        FILE **convfiles;
//        #endif
        #endif

        cmd=argv[0];
	curmodule="tablix";

        strcpy(prefix, "./");
	verbosity=102;

        #ifndef HAVE_LIBPVM3
        restore=0;
	maxlocals=1;
	nodereq=4;

	timeout=0;

//        #elifdef HAVE_LIBPVM3
//        restore=0;
//	maxlocals=1;
//	nodereq=4;
//
//	timeout=0;
        #endif

	#ifdef HAVE_SETLOCALE
	loc=setlocale(LC_ALL, "");
	#else
	loc="C";
	#endif

	#if ENABLE_NLS && !defined DEBUG
	/* This won't compile without -O2. */
	bindtextdomain(PACKAGE, LOCALEDIR);
	textdomain(PACKAGE);
	#endif

        printf(_("TABLIX version %s, PGA general timetable solver\n"), VERSION);
        printf("Copyright (C) 2002-2006 Tomaz Solc\n");

        while ((c=getopt(argc, argv, "hn:vo:rl:d:t:p:i:"))!=-1) {
                switch (c) {
                        case 'v': print_copyright();
                        case 'h':
                        case '?': print_syntax();
                                  exit(0);
                        case 'o': strncpy(prefix, optarg, 256);
                                  break;
                        #ifndef HAVE_LIBPVM3
                        case 'l': sscanf(optarg, "%d", &maxlocals);
                                  break;
                        case 'n': sscanf(optarg, "%d", &nodereq);
                                  break;
                        case 'r': restore=1;
                                  break;
			case 'd': sscanf(optarg, "%d", &verbosity);
				  verbosity+=100;
				  break;
			case 't': sscanf(optarg, "%d", &timeout);
				  break;

//                        #elifdef HAVE_LIBPVM3
//                        case 'l': sscanf(optarg, "%d", &maxlocals);
//                                  break;
//                        case 'n': sscanf(optarg, "%d", &nodereq);
//                                  break;
//                        case 'r': restore=1;
//                                  break;
//			case 'd': sscanf(optarg, "%d", &verbosity);
//				  verbosity+=100;
//				  break;
//			case 't': sscanf(optarg, "%d", &timeout);
//				  break;
                        #endif
                }
        }

        printf("\n");

	if(loc==NULL) {
		//info(_("Locale not supported by C library. "
		//			"Using the fallback 'C' locale."));
		loc="C";
	}

        if (!(optind<argc)) {
                fatal(_("Missing file name. Try '%s -h' for more information."), cmd);
        }

        #ifdef HAVE_LIBPVM3
        execvp("tablix2_kernel", argv);

        perror(cmd);
        exit(1);

//        #elifndef HAVE_LIBPVM3
//        execvp("tablix2_kernel", argv);
//
//        perror(cmd);
//        exit(1);
        #else

        /*** Start nodes ***/

	assert(argv[argc]==NULL);
    
	node_startall(nodereq, &argv[1]);

        if(nodenum ==0) {
                //pvm_exit();
                MPI_Finalize();
                fatal(_("all nodes failed."));
        }

	snprintf(fn, 256, "%s %s",
		ngettext("PGA using %d node", "PGA using %d nodes", nodenum),
		ngettext("on %d host", "on %d hosts", hostnum));

        //info(fn, nodenum, hostnum);

	if (maxlocals>0) {
		debug(ngettext("maximum %d node will do local search",
			       "maximum %d nodes will do local search",
			       maxlocals),
		      maxlocals);
	} else if (maxlocals==0) {
		debug(_("local search disabled"));
	} else {
		debug(_("local search enabled on all nodes"));
	}
        debug(_("multicasting XML data"));

        /*** Send nodes XML data ***/

        buff=malloc(LINEBUFFSIZE);
        module=malloc(LINEBUFFSIZE);
        if (buff==NULL||module==NULL) {
                perror(cmd);
                //for(c=0;c<nodenum;c++) pvm_kill(nodetid[c]);
                //pvm_exit();
                MPI_Finalize();
                exit(1);
        }

        
        
        strsize=strlen(loc);
        
        for(count=0; count<nodenum; count++)
        {
            mpibcast=MPI_Send(&strsize, 1, MPI_INT, count, MSG_PARAMS, MPI_COMM_WORLD);
            mpibcast2=MPI_Send(&loc, strsize, MPI_CHAR, count, MSG_PARAMS, MPI_COMM_WORLD);
        }
        
	if (file_mcast(argv[optind], nodetid, nodenum, MSG_XMLDATA)) {
                perror(cmd);
		error(_("Failed to send problem description '%s' "
					"to computing nodes"), argv[optind]);
		free(buff);
		free(module);
                MPI_Finalize();
                //pvm_exit();
                exit(1);
        }

        /*** Restore populations from a file if requested ***/
        if (restore) {
                debug(_("Restoring saved populations"));

                running=0;
                for(n=0;n<nodenum;n++) {
                        snprintf(fn, 256, "%ssave%d.txt", prefix, n);

			pop=population_load(fn);

			if(pop==NULL) {

				//info(_("Failed to load %s"), fn);
                a=0;

                msend=MPI_Send(&a, 1, MPI_INT, nodetid[n], MSG_RESTOREPOP, MPI_COMM_WORLD);
			} else {
			a=1;
			msend=MPI_Send(&a, 1, MPI_INT, nodetid[n], MSG_RESTOREPOP, MPI_COMM_WORLD);

				population_send(pop, nodetid[n], MSG_RESTOREPOP);
				population_free(pop);

	                        running++;
			}
                }

		if (running<nodenum) {
                       	debug(_("%d nodes randomized"), nodenum-running);
		}
        } else {
        a=0;
        int count;
                for(count=0; count<nodenum; count++)
                {
                    msend2=MPI_Send(&a, 1, MPI_INT, count, MSG_RESTOREPOP, MPI_COMM_WORLD);
                }
        }

        debug(_("Initializing nodes"));

        /*** Send nodes TID to send migration to ***/

        for(c=0; c<nodenum; c++) {
		node_update(c);

        	if(verbosity>=MSG_DEBUG) {
	                printf("%x ", nodetid[c]);
	                fflush(stdout);
		}
        }

       	//if(verbosity>=MSG_DEBUG) printf("\n");
        //debug(_("parent ( %x ) listening\n"), pvm_mytid());

	/*** Set up files for saving convergence info ***/

        #ifdef HAVE_CONV
        convfiles=malloc(nodenum*sizeof(*convfiles));
        if (convfiles==NULL) {
                error(_("Can't allocate memory"));
                perror(cmd);
                free(buff);
                /* We don't need to kill the nodes here. */
                /* They should detect this. */
                //pvm_exit();
                MPI_Finalize();
                exit(1);
        }

        for(c=0;c<nodenum;c++) {
                snprintf(fn, 256, "%sconv%d.txt", prefix, c);
                convfiles[c]=fopen(fn, restore?"a":"w");

                if (convfiles[c]==NULL) {
                	error(_("Can't open file %s for writing"), fn);
                        perror(cmd);
                        free(buff);
                        for(a=0;a<c;a++) fclose(convfiles[a]);
                        free(convfiles);
                        MPI_Finalize();
                        //pvm_exit();
                        exit(1);
                }
        }
        #endif

	/*** Set up timer and signal handlers ***/

    ctrlc=0;
	timeout_reached=0;
    signal(SIGINT, sighandler);
	signal(SIGALRM, sighandler);
	alarm(timeout*60);

	gnum=-1;
	subtotals=NULL;		/* to keep gcc happy */

	gettimeofday(&start, NULL);

	/*** Check if the standard output is on a tty ***/

	c=fileno(stdout);
	if(isatty(c)) {
		cnt_newline=0;
	} else {
		cnt_newline=1;
	}

	cnt_recv=0;
	cnt_stopped=0;
	cnt_total=nodenum;
	counter_init();

	/*** MAIN LOOP STARTS HERE ***/

	numlocals=0;
	

    while (cnt_stopped<nodenum) {
        recvi=MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &st1);
        sender=st1.MPI_SOURCE;
        msgtag=st1.MPI_TAG;
        
		if(msgtag!=MSG_REPORT&&msgtag!=MSG_LOCALSYN) {
			counter_clear();
		}

        sendernum=node_find(sender);

        if(sendernum==-1&&msgtag!=MSG_NODEKILL) {
            error(_("received message from unknown node %x"), sender);
        }

        if (msgtag==MSG_MODINFO) {
            unpack1=MPI_Unpack(*fn, 256, *position, &a, 1, MPI_INT, MPI_COMM_WORLD);

			debug(_("%x has %d fitness functions"), sender, a);

			if(gnum<0) {
				gnum=a;
				subtotals=malloc(sizeof(*subtotals)*gnum);
			} else if(gnum!=a) {
				error(_("node %x has a different number of modules than others"), sender);
			}

            #ifdef HAVE_CONV
			if(!restore) {
                fprintf(convfiles[sendernum], "# Gen.\tFitness\tOK");
				for(c=0;c<gnum;c++) {
                    unpack2=MPI_Unpack(*fn, 256, *position, &strsize, 1, MPI_INT, MPI_COMM_WORLD);
                    unpack3=MPI_Unpack(*fn, 256, *position, &fn, strsize, MPI_CHAR, MPI_COMM_WORLD);
                    unpack4=MPI_Unpack(*fn, 256, *position, &a, 1, MPI_INT, MPI_COMM_WORLD);
					fprintf(convfiles[sendernum], "\t%s%s", fn, a?" (M)":"");
				}
				fprintf(convfiles[sendernum], "\n");

                fflush(convfiles[sendernum]);
			}
            #endif
		} else
		  if (msgtag==MSG_FATAL) {
            unpack5=MPI_Unpack(*buff, LINEBUFFSIZE, *position, &strsize, 1, MPI_INT, MPI_COMM_WORLD);
            unpack6=MPI_Unpack(*buff, LINEBUFFSIZE, *position, &module, strsize, MPI_CHAR, MPI_COMM_WORLD);
            unpack7=MPI_Unpack(*buff, LINEBUFFSIZE, *position, &strsize, 1, MPI_INT, MPI_COMM_WORLD);
            unpack8=MPI_Unpack(*buff, LINEBUFFSIZE, *position, &buff, strsize, MPI_CHAR, MPI_COMM_WORLD);

            msg_new(sender, _("%s: FATAL: %s"), module, buff);
            cnt_stopped++;
            node_stop(sendernum);
        } else
		  if (msgtag==MSG_NODEKILL) {
//              unpack9=MPI_Unpack(*fitbuff, 256, *position, &a, 1, MPI_INT, MPI_COMM_WORLD);
//			   pvm_upkint(&a, 1, 1);
//			   c=node_find(a);
//			   if(c==-1) {
//			       error(_("received quit notification "
//						"from unknown node %x"), a);
//			   } else
//              if(!nodeinfo[c].dead) {
//			      error(_("node %x has unexpectedly quit"), a);
//				  cnt_stopped++;
//				  node_stop(c);
//			   }
		  } else
            if (msgtag==MSG_INFO||msgtag==MSG_ERROR||msgtag==MSG_DEBUG) {
                unpack9=MPI_Unpack(*buff, LINEBUFFSIZE, *position, &strsize, 1, MPI_INT, MPI_COMM_WORLD);
                unpack10=MPI_Unpack(*buff, LINEBUFFSIZE, *position, &module, strsize, MPI_CHAR, MPI_COMM_WORLD);
                unpack11=MPI_Unpack(*buff, LINEBUFFSIZE, *position, &strsize, 1, MPI_INT, MPI_COMM_WORLD);
                unpack12=MPI_Unpack(*buff, LINEBUFFSIZE, *position, &buff, strsize, MPI_CHAR, MPI_COMM_WORLD);

			    msg_new(sender, "%s: %s", module, buff);
            } else
                if (msgtag==MSG_RESULTDATA) {
                    //info(_("node %x is uploading the result"), sender);
                    snprintf(fn, 256, "%sresult%d.xml", prefix, sendernum);

			         if (file_recv(fn, sender, MSG_RESULTDATA)) {
                        perror(cmd);
                    }

                    cnt_recv++;
	           		cnt_stopped++;
                       
                    node_stop(sendernum);
                } else
                    if (msgtag==MSG_REPORT&&gnum<0) {
                        error(_("node %x is reporting too soon"), sender);
                    } else
		              if (msgtag==MSG_REPORT) {
                        recvi5=MPI_Recv(&c, 1, MPI_INT, MPI_ANY_SOURCE, MSG_REPORT, MPI_COMM_WORLD, &st2);
                        recvi6=MPI_Recv(&a, 1, MPI_INT, MPI_ANY_SOURCE, MSG_REPORT, MPI_COMM_WORLD, &st2);
                        recvi7=MPI_Recv(&d, 1, MPI_INT, MPI_ANY_SOURCE, MSG_REPORT, MPI_COMM_WORLD, &st2);
                        recvi8=MPI_Recv(&subtotals, gnum, MPI_INT, MPI_ANY_SOURCE, MSG_REPORT, MPI_COMM_WORLD, &st2);

			            if(verbosity>=MSG_REPORT) {
                            counter_update(sender, c, d, a);
			            }
			            nodeinfo[sendernum].generations++;
                        
                        #ifdef HAVE_CONV
                            fprintf(convfiles[sendernum], "%d\t%d\t%d", a, c, d);
                            for(c=0;c<gnum;c++)
                                fprintf(convfiles[sendernum], "\t%d", subtotals[c]);
                            
                            fprintf(convfiles[sendernum], "\n");
                            fflush(convfiles[sendernum]);
                        #endif
                    } else
                    if (msgtag==MSG_POPDATA) {
                        //info(_("receiving population from %x"), sender);

                        snprintf(fn, 256, "%ssave%d.txt", prefix, sendernum);

                        pop=population_recv(sender, MSG_POPDATA);
                        population_save(fn, pop);
                        population_free(pop);

                        cnt_stopped++;
                        node_stop(sendernum);
		          } else
                    if (msgtag==MSG_LOCALSYN) {
                        recvi4=MPI_Recv(&c, 1, MPI_INT, MPI_ANY_SOURCE, MSG_REPORT, MPI_COMM_WORLD, &st2);

                        if(!c) {
                            numlocals--;
                            node_restart(sendernum);
                        } else {
                            if(numlocals!=maxlocals) {
                                c=1;
                                numlocals++;
                                node_stop(sendernum);
                            } 
                            else {
                                c=0;
                            }
                            
                            MPI_Send(&c, 1, MPI_INT, sender, MSG_LOCALACK, MPI_COMM_WORLD);
                        }
                    } else
                    {
                            error(_("unknown message received from %x. Perhaps a version mismatch?"), sender);
                    }
    }

	msg_flush();

    #ifdef HAVE_CONV
    for(c=0;c<nodenum;c++) fclose(convfiles[c]);
    
    free(convfiles);
    #endif

    if (ctrlc) {
        if(timeout_reached) {
            error(_("Time limit reached"));
        } else {
            error(_("Received interrupt"));
        }
	}

	/*if(cnt_recv<1) {
		info(_("No results were received"));
	} else if(cnt_recv<cnt_stopped) {
        info(_("Some results were received"));
        } else {
		info(_("All results were received"));
	   }*/

	if(verbosity>=MSG_INFO&&cnt_recv>0) {
		gettimeofday(&end, NULL);

		d=end.tv_sec-start.tv_sec;

		printf("\n");
		printf(_("node\tgenerations/minute\n"));
		a=0;
		for(c=0;c<nodenum;c++) {
			sum=((double) nodeinfo[c].generations*60)/d;
			printf("%x\t%.1f\n", nodetid[c], sum);
			a+=nodeinfo[c].generations;
		}
		sum=((double) a*60)/d;
		printf(_("\ntotal generations per minute: %.1f\n"), sum);
		printf(_("total time: %02d:%02d:%02d\n"), d/3600, (d/60)%60, d%60);
	}

    free(nodetid);
    free(nodeinfo);
    free(buff);

    //pvm_exit();
    MPI_Finalize();
    if (ctrlc) {
        exit(1);
	} else if(cnt_recv==0) {
		exit(2);
	   } else {
		  exit(0);
	   }
   #endif
}

#ifndef HAVE_LIBPVM3
void sighandler(int num)
{
    int mpisend;
    
    int c;

	if(num==SIGALRM) {
		timeout_reached=1;
	}
        if (ctrlc<1) {
                c=0;
                for(count=0; count<nodenum; count++)
                {
                    mpisend=MPI_Send(&c, 1, MPI_INT, count, MSG_SENDPOP, MPI_COMM_WORLD);
                }
                ctrlc++;
        } else {
                exit(1);
        }
}

//#elifdef HAVE_LIBPVM3
//  void sighandler(int num)
//  {
//          int c;
//
//    if(num==SIGALRM) {
//      timeout_reached=1;
//    }
//          if (ctrlc<1) {
//                  pvm_initsend(0);
//                  c=0;
//                  pvm_pkint(&c, 1, 1);
//                  pvm_mcast(nodetid, nodenum, MSG_SENDPOP);
//                  ctrlc++;
//          } else {
//                  exit(1);
//          }
//  }
#endif
