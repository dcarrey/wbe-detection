EXECDIR=../exec

all: clean hgt scripts

hgt: strsep.cpp hgt.cpp fonctions.cpp utils_tree.cpp structures.h bioling.cpp
	@echo  "\n=> Building hgt.cpp -> hgt"
	g++ hgt.cpp -o hgt

scripts:
	-mkdir $(EXECDIR)
	cp hgt $(EXECDIR)
	cp ageWBE.pl $(EXECDIR)
	cp run_wbe.pl $(EXECDIR)
	cp PostTraitement.pm $(EXECDIR)
	cp ../examples/* $(EXECDIR)

clean:
	@echo  "\n=> Cleaning executable and object files"
	rm -f hgt *.o
	mkdir -p $(EXECDIR)
