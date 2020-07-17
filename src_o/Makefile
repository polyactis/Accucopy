
CXXTARGETS	=
#CXXLDFLAGS 	= $(BoostLib)
CTARGETS	=
SharedLibTargets	=
StaticLibTargets =


SRCS	= infer.cpp read_para.cpp BaseGADA.cc GADA.cc model_selection.cpp

ExtraTargets = infer GADA

infer:	%:	%.o read_para.o prob.o BaseGADA.o format.o model_selection.o
	$(CXXCOMPILER) $< read_para.o prob.o BaseGADA.o format.o model_selection.o $(CXXFLAGS) -o $@ $(CXXLDFLAGS) -lgsl -lgslcblas $(BoostLib)

GADA:   %:   %.o BaseGADA.o BaseGADA.h read_para.o format.o
	$(CXXCOMPILER) $< BaseGADA.o read_para.o format.o $(CXXFLAGS) -o $@ -lm $(CXXLDFLAGS) $(BoostLib)

recall_precision:	%:	%.o
	$(CXXCOMPILER) $< $(CXXFLAGS) -o $@ $(CXXLDFLAGS)

README.pdf: README.tex
	-pdflatex README.tex

#:= is different from =. The latter will cause the function evaluation every time the make variable is invoked.
#currentUnixTime:=$(shell echo "import time; print str(time.time()).replace('.', '_')"|python)
currentTime:=$(shell echo "from datetime import datetime; print str(datetime.now()).replace(' ', '_').replace(':', '')"|python)


debug: all ../src/main.rs
	-rm -rf debug
	-mkdir debug
	-mkdir -p ../target/debug/
	cp ./main.py ./debug/main.py
	chmod 755 ./debug/main.py
	cargo update
	cargo build
	-git checkout -- ../src/main.rs
	cp -r __init__.py ../LICENSE GADA ../target/debug/maestre configure infer plotCPandMCP.py plot_autocor_diff.py plot_model_select_result.py plot_coverage_after_normalization.py plot_tre.py plot_snp_maf_exp.py plot_snp_maf_peak.py debug/
	tar -cavf debug.$(currentTime).tar.gz debug/

release: all ../src/main.rs
	-rm -rf release
	-mkdir release
	-mkdir -p ../target/release/
	cp ./main.py ./release/main.py
	chmod 755 ./release/main.py
	cargo update
	cargo build --release
	-git checkout -- ../src/main.rs
	cp -r __init__.py ../LICENSE GADA ../target/release/maestre configure infer plotCPandMCP.py plot_autocor_diff.py plot_model_select_result.py plot_coverage_after_normalization.py plot_tre.py plot_snp_maf_exp.py plot_snp_maf_peak.py release/
	tar -cavf release.$(currentTime).tar.gz release/


include ../Makefile.common

# leave it at last to avoid being overridden by Makefile.common
clean: common_clean
	-rm README.aux README.log
	-rm -r release*
	-rm -r debug*
