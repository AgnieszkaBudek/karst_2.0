
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
	FLIBS	= -lblas -L/usr/local/gfortran/lib 
endif

AR      = ar vr 


COPTS	= -Wall -pedantic -O2 
MINCLUDES	= -I $(MDIR)/include -I $(MDIR)/libseq 

CFLAGS	= $(COPTS) $(MINCLUDES) 


OBJ	= matrix_solvers.o triangularization.o main.o
INC	= algorithms.h 


LINKOPTS = $(COPTS)
MLIBS	= -lm -ldmumps -lmpiseq -lpthread -lpord -lstdc++  -lgfortran $(FLIBS)
MFLAGS	= $(LINKOPTS) $(MLIBS) -L $(MDIR)/lib -L $(MDIR)/libseq
CXXFLAGS= -std=c++11 -Weffc++ -Wall -Wextra -Wfatal-errors -Wunused -pedantic -w -Winline 
TMP     = -fno-rtti -ggdb -D_GLIBCXX_DEBUG -lsfml-graphics -lsfml-window -lsfml-system


.cc.o: 
	$(CXX) $(CXXFLAGS) -c $< 

.c.o: 
	$(CC) $(CFLAGS) -c $< 


algorithms: $(OBJ)	
#	$(FL) -o ../bin/algorithms $^ $(MFLAGS)
	$(AR) libalgorithms.a matrix_solvers.o triangularization.o 
	$(CXX) -v -o ../bin/algorithms $^ $(MFLAGS)  $(CXXFLAGS)


clean:
	rm -f *.o *.a


all: 
	mkdir -p ../bin
	make algorithms
