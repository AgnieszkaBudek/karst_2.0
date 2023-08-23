/*
 *
 *  This file is part of MUMPS 5.6.0, released
 *  on Wed Apr 19 15:50:57 UTC 2023
 *
 *
 *  Copyright 1991-2023 CERFACS, CNRS, ENS Lyon, INP Toulouse, Inria,
 *  Mumps Technologies, University of Bordeaux.
 *
 *  This version of MUMPS is provided to you free of charge. It is
 *  released under the CeCILL-C license 
 *  (see doc/CeCILL-C_V1-en.txt, doc/CeCILL-C_V1-fr.txt, and
 *  https://cecill.info/licences/Licence_CeCILL-C_V1-en.html)
 *
 */
#include "mpi.h"
#include "elapse.h"
MPI_Comm LIBSEQ_CALL MPI_Comm_f2c(LIBSEQ_INT comm)
{
  return 0;
}

LIBSEQ_INT LIBSEQ_CALL MPI_Init(LIBSEQ_INT *pargc, char ***pargv)
{
  return 0;
}

LIBSEQ_INT LIBSEQ_CALL MPI_Comm_rank(MPI_Comm comm, LIBSEQ_INT *rank)
{
  *rank=0;
  return 0;
}
LIBSEQ_INT LIBSEQ_CALL MPI_Finalize(void)
{
   return 0;
}

/* Internal: for MPI_IS_IN_PLACE tests from Fortran */

void LIBSEQ_CALL MUMPS_CHECKADDREQUAL(char *a, char*b, LIBSEQ_INT *i)
{
  if (a - b == 0)
   {
     *i=1;
   }
 else
   {
     *i=0;
   }
}

void LIBSEQ_CALL MUMPS_CHECKADDREQUAL_(char *a, char*b, LIBSEQ_INT *i)
 {
   MUMPS_CHECKADDREQUAL(a,b,i);
 }
void LIBSEQ_CALL mumps_checkaddrequal_(char *a, char*b, LIBSEQ_INT *i)
 {
   MUMPS_CHECKADDREQUAL(a,b,i);
 }
void LIBSEQ_CALL mumps_checkaddrequal__(char *a, char*b, LIBSEQ_INT *i)
 {
   MUMPS_CHECKADDREQUAL(a,b,i);
 }

double LIBSEQ_CALL MPI_Wtime()
{
  double val;
  mumps_elapse(&val);
  return val;
}
