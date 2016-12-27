EXECDIR=../exec

all: clean hgt scripts

hgt: strsep.cpp newhgt.cpp fonctions.cpp utils_tree.cpp structures.h bioling.cpp
	@echo  "\n=> Building newhgt.cpp -> hgt"
	g++ newhgt.cpp -o hgt

scripts: 
	cp hgt $(EXECDIR)
	cp ageWBE.pl $(EXECDIR)
	cp run_hgt.pl $(EXECDIR)
	cp postTraitement.pl $(EXECDIR)

clean:
	@echo  "\n=> Cleaning executable and object files"
	rm -f hgt *.o
	mkdir -p $(EXECDIR) 
