
srcdir      = ./src

all: src_all


ifneq ("$(wildcard /Users/agnieszka/Desktop/KARST/symulacja/test)","") 
	bash -c "cp ./bin/karst /Users/agnieszka/Desktop/KARST/symulacja/test"
else ifneq ("$(wildcard ~/KARST/symulacja/test)","")
	bash -c "cp ./bin/karst ~/KARST/symulacja/test"
else  
	bash -c "echo "Warning: Directory ..../KARST/symulacja/test does not exists.""
endif		 

src_all:
	$(MAKE) all -C $(srcdir)
	
clean:
	$(MAKE) clean  -C $(srcdir)
