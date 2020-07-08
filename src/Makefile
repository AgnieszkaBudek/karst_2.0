UNAME := $(shell uname)

ifeq ($(UNAME), Linux)
	CC		= clang        
	CXX 	= clang++                                            
	FL		= gfortran 
	MDIR	=  ~/.nix-profile/opt/mumps
	FLIBS	= -L ~/.nix-profile/lib -lgfortran -lquadmath -lblas
else 
	CC		= /usr/local/Cellar/gcc/9.3.0/bin/gcc-9        
	CXX 	= g++                                            
	FL		= /usr/local/Cellar/gcc/9.3.0/bin/gfortran-9 
	MDIR	= ~/Desktop/STUDIA/KRAS/MUMPS_4.7.3
	FLIBS	= -lblas -L/usr/local/gfortran/lib -lgfortran 
endif	

#MUMPS linking
MLIBS	= -ldmumps -lmpiseq -lpthread -lpord  -lstdc++ $(FLIBS)
ALGDIR  = ../algorithms

COPTS	= -Wall -pedantic -O2  
CINCLUDES	= -I ../include -I ../algorithms 
MINCLUDES	= -I ../include -I $(MDIR)/include -I $(MDIR)/libseq

CFLAGS	 = $(COPTS) $(CINCLUDES) $(MINCLUDES) 
CXXFLAGS = -O2 -std=c++11 -Weffc++ -Wall -Wextra -Wfatal-errors -Wunused -pedantic -w -Winline -O2
TMP      = -fno-rtti -ggdb -D_GLIBCXX_DEBUG -lsfml-graphics -lsfml-window -lsfml-system


OBJ	= network.o main.o printing.o printing_net.o node.o pore.o grain.o merging.o initial_modifications.o\
evolution.o network_maintenance.o setup.o dissolution.o control_and_verification.o import_export.o network_generation.o\
tests.o analysis.o pattern_analysis.o save_in_VTK.o transversal_diffussion_version.o transversal_diffussion_version_with_time.o
INC	= network.h node.h pore.h grain.h printing.h constants.h algorithms_cc.h


LINKOPTS  = $(COPTS) -lstdc++ $(FLIBS)  -lgfortran
LIBS      = -lm -L ../algorithms -lalgorithms
MFLAGS    = $(MLIBS) -L $(MDIR)/lib -L $(MDIR)/libseq
LFLAGS    = $(LINKOPTS) $(LIBS) $(MFLAGS) 

.c.o: 
	$(CC) $(CFLAGS) -c $< 

.cc.o: $(INC)
	$(CXX) $(CXXFLAGS) $(CINCLUDES) -c $<
   	 

karst: $(OBJ) 	
	$(CXX) -v -o ../bin/karst $^ $(LFLAGS) $(CXXFLAGS)

algorithms:
	$(MAKE) all -C $(ALGDIR)

clean:
	rm -f *.o

clean_alg:
	rm -f ../algorithms/*.o ../algorithms/*.a

all: 
	mkdir -p ../bin
	make algorithms
	make karst

