COMPILER = g++
FLAGS = -O3
all: latticeeasy.h model.h parameters.h ffteasy.o latticeeasy.o evolution.o initialize.o output.o
	$(COMPILER) $(FLAGS) ffteasy.o latticeeasy.o evolution.o initialize.o output.o -lm -o latticeeasy

clean:
	rm -f *.o out output.txt out.txt latticeeasy *~ *.dat core*

cleaner:
	rm -f *.o out output.txt out.txt latticeeasy *~ *.dat *.img

ffteasy.o: ffteasy.cpp
	$(COMPILER) -c $(FLAGS) ffteasy.cpp

latticeeasy.o: latticeeasy.cpp latticeeasy.h model.h parameters.h
	$(COMPILER) -c $(FLAGS) latticeeasy.cpp

evolution.o: evolution.cpp latticeeasy.h model.h parameters.h
	$(COMPILER) -c $(FLAGS) evolution.cpp

initialize.o: initialize.cpp latticeeasy.h model.h parameters.h
	$(COMPILER) -c $(FLAGS) initialize.cpp

output.o: output.cpp latticeeasy.h model.h parameters.h
	$(COMPILER) -c $(FLAGS) output.cpp


