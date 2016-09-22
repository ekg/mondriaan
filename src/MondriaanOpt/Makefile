
include ../../mondriaan.mk

LIBOBJECTS = matching.o matrix.o solution.o solver.o
OBJECTS= MondriaanOpt.o errors.o options.o ${LIBOBJECTS}
HEADERS= errors.h matching.h matrix.h options.h solution.h solver.h

all: MondriaanOpt lib/libMondriaan${MONDRIAANMAJORVERSION}Opt.a

MondriaanOpt: $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) -o MondriaanOpt $(OBJECTS) $(LIBS)
	cp MondriaanOpt ${MONDRIAANHOMEDIR}/tools/MondriaanOpt

:c.o:
	$(CC) -c $(CFLAGS) $<

clean:
	rm -f *.o
	rm -f MondriaanOpt ${MONDRIAANHOMEDIR}/tools/MondriaanOpt

veryclean:
	rm -f *.o
	rm -f MondriaanOpt ${MONDRIAANHOMEDIR}/tools/MondriaanOpt


#archive for linking
lib/libMondriaan${MONDRIAANMAJORVERSION}Opt.a: ${LIBOBJECTS}
	mkdir -p ../lib
	ar rcs libMondriaan${MONDRIAANMAJORVERSION}Opt.a $^
	mv libMondriaan${MONDRIAANMAJORVERSION}Opt.a ../lib/
