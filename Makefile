# Mondriaan Makefile

.SILENT:

all: lib lib_opt tool

include mondriaan.mk

clean:
	@echo ==== Cleaning Mondriaan package ====
	(cd src; make clean)
	(cd src/MondriaanOpt; make clean)
	(cd tools; make clean)
	(cd tests; make clean)

veryclean:
	@echo ==== Cleaning Mondriaan package thoroughly ====
	(cd src; make veryclean)
	(cd src/MondriaanOpt; make veryclean)
	(cd tools; make veryclean)
	(cd tests; make veryclean)

lib:
	@echo === Building Mondriaan library ===
	@echo ${MONDRIAANCOMPILEINFO}
	(cd src; make all)

lib_opt:
	@echo === Building MondriaanOpt library ===
	@echo ${MONDRIAANCOMPILEINFO}
	(cd src/MondriaanOpt; make all)

tool:
	@echo ==== Building Mondriaan tools ====
	@echo ${MONDRIAANCOMPILEINFO}
	(cd src; make all)
	(cd tools; make all)

test:
	@echo ==== Building Mondriaan tests ====
	@echo ${MONDRIAANCOMPILEINFO}
	(cd src; make all)
	(cd tools; make all)
	(cd tests; make all)

