MEX       = /opt/matlab/bin/mex -g
MUMPS_DIR = /home/jylexcel/MUMPS

LMETISDIR = /home/jylexcel/metis-4.0
LMETIS = -L$(LMETISDIR) -lmetis
LPORDDIR = $(MUMPS_DIR)/PORD/lib
LPORD = -L$(LPORDDIR) -lpord
LORDERINGS = $(LPORD) $(LMETIS)

LIBSEQ  = -L$(MUMPS_DIR)/libseq -lmpiseq

LIBFORTDIR= /home/afevre/compil/g95-install/lib/gcc-lib/i686-pc-linux-gnu/4.0.1
LIBFORT = $(LIBFORTDIR)/libf95.a $(LIBFORTDIR)/libgcc.a 

OPTC    = -O 

