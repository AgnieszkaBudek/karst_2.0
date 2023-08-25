include make.inc

#MUMPS_DVER = DROIDE_DBL_SEQ
#MUMPS_ZVER = DROIDE_cmplx16_SEQ


INCMUMPS = -I../include 
LIBMUMPS = -L../lib -ldmumps -lzmumps 

INCTOT = $(IORDERINGS) $(INCSEQ) $(INCMUMPS)
LIBTOT = $(LIBMUMPS) $(LORDERINGS) $(LIBSEQ) $(LIBFORT) -lpthread

all:
	$(MEX) $(OPTC) dmumpsmex.c $(INCTOT) $(LIBTOT)
	$(MEX) $(OPTC) zmumpsmex.c $(INCTOT) $(LIBTOT)

