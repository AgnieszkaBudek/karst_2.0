
srcdir	=	./src

all: src_all
	bash -c "cp ./bin/karst ./tests"

#ifneq ("$(wildcard /Users/agnieszka/Desktop/KARST/symulacja/test)","")
#	bash -c "cp ./bin/karst /Users/agnieszka/Desktop/KARST/symulacja/test"
#else ifneq ("$(wildcard ~/KARST/symulacja/test)","")
#	bash -c "cp ./bin/karst ./tests"
#else
#	bash -c "echo "Warning: Directory ..../KARST/symulacja/test does not exists.""
#endif

# Default target
.DEFAULT:
    #@echo "No target specified. To build a program, use: make <program_name>"


src_all:
	$(MAKE) all -C $(srcdir)

clean:
	$(MAKE) clean  -C $(srcdir)
