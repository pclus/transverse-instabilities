CC=gcc
CFLAGS=-Wall -Wno-maybe-uninitialized -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-result -g -pg -lm -lgsl -lgslcblas -O3 -fopenmp
SOURCES=visualization.c system.c netdynamics.c
DEPENDENCIES=$(SOURCES:.c=.h)
#GLOBAL_DEP=macros.h
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=../netdyn


$(EXECUTABLE): $(OBJECTS) $(SOURCES) $(DEPENDENCIES)
	$(CC) $(OBJECTS) -o $@ $(CFLAGS) 
        
build: $(SOURCES) $(DEPENDENCIES)
	$(CC) -c $(SOURCES) $(CFLAGS) -DVISUAL=true
        
visual: $(SOURCES) $(DEPENDENCIES)
	$(CC) -c $(SOURCES) $(CFLAGS) 

clean:
	rm -f *.o
