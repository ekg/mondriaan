CC=gcc
CFLAGS= -Wall -O3 -ffast-math -funroll-loops -fomit-frame-pointer -std=c99 -I../include
LFLAGS= -L../lib/ -lm -lMondriaan4

OBJECTS= MondriaanOpt.o errors.o matching.o matrix.o options.o solution.o solver.o
HEADERS= errors.h matching.h matrix.h options.h solution.h solver.h

all: MondriaanOpt

MondriaanOpt: $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) -o MondriaanOpt $(OBJECTS) $(LFLAGS)

:c.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o
	rm -f MondriaanOpt

veryclean:
	rm -f *.o
	rm -f MondriaanOpt
