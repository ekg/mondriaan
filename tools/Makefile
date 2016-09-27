# Mondriaan Tools Makefile

include ../mondriaan.mk

#scrap default rule:
% :: %.c ;

#targets
TARGETS := Mondriaan MondriaanPlot Profile
LIB = ${MONDRIAANLIBDIR}/libMondriaan${MONDRIAANMAJORVERSION}.a
ifdef MATLABHOMEDIR
TARGETS := ${TARGETS} MatlabMondriaan MatlabMondriaanOpt
endif

all: ${TARGETS}

%.target: ${OBJDEPS}
	${CC} ${CFLAGS} -o ${@:%.target=%} ${OBJDEPS} ${LIBS}

%.${MEXSUFFIX}: ${OBJDEPS}
	${MEX} ${MEXFLAGS} -output ${@:%.mextarget=%.${MEXSUFFIX}} ${OBJDEPS} ${LIBS}

%.o: %.c %.h
	${CC} ${CFLAGS} -c -o $@ ${filter %.c,$^}

%.mex.o: %.c
	${CC} ${MEXCFLAGS} -c -o $@ $^

%.o: %.c
	${CC} ${CFLAGS} -c -o $@ $^

%.cstub:
	rm -f ${@:%.cstub=%}.o

%.vcstub:
	rm -f ${@:%.vcstub=%}

veryclean: clean ${TARGETS:%=%.vcstub}
	rm -f *.${MEXSUFFIX}

clean: ${TARGETS:%=%.cstub}
	rm -f `find . -name 'core'`

Profile: Profile.o ${LIB}
	make -r OBJDEPS='${^:%.a=}' $@.target

Mondriaan: Mondriaan.o ${LIB}
	make -r OBJDEPS='${^:%.a=}' $@.target

MondriaanPlot: MondriaanPlot.o ${LIB}
	make -r OBJDEPS='${^:%.a=}' $@.target

MatlabMondriaan: MatlabMondriaan.mex.o MatlabHelper.mex.o ${LIB}
	make -r OBJDEPS='${^:%.a=}' $@.${MEXSUFFIX}

MatlabMondriaanOpt: MatlabMondriaanOpt.mex.o MatlabHelper.mex.o ${LIB} ${MONDRIAANLIBDIR}/libMondriaan${MONDRIAANMAJORVERSION}Opt.a
	make -r OBJDEPS='${^:%.a=}' LIBS='-lMondriaan4Opt ${LIBS}' $@.${MEXSUFFIX}

