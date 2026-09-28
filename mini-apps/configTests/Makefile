# configTests - unit tests for the config file / expression parsing pipeline.
#
# Standalone miniApp: builds its own test binaries against the (few) Vlasiator
# sources they need, without touching the toplevel build. The CI job
# "configTests_github" in .github/workflows/github-ci.yml runs `make test`.
#
#   make            build both test binaries
#   make test       build and run the whole suite
#   make clean      remove objects and binaries

ARCH = ${VLASIATOR_ARCH}
ifeq ($(strip $(ARCH)),)
ARCH := github_actions
endif

include ../../MAKE/Makefile.${ARCH}

# The arch file brings in heavy -march/-O3 settings for the solver binary;
# tests want warnings instead.
CXXFLAGS := -O1 -std=c++20 -fabi-version=0 -Wall -Wextra -I../../ $(filter-out -march=% -mavx -g -O3 -DUSE_JEMALLOC -DPAPI_MEM, $(CXXFLAGS))

TESTS = glossa_tests config_reader_test
OBJS  = readparameters.o version.o

default: $(TESTS)

readparameters.o: ../../readparameters.cpp ../../readparameters.h ../../version.h ../../qdparser.h ../../common.h ../../glossa.hpp
	$(CMP) ${CXXFLAGS} -c ../../readparameters.cpp -o $@

version.o: ../../version.cpp ../../version.h
	$(CMP) ${CXXFLAGS} -c ../../version.cpp -o $@

# Regenerate version.cpp the same way the toplevel Makefile does, if needed.
../../version.cpp:
	make -C../.. version.cpp

glossa_tests: glossa_tests.cpp ../../glossa.hpp
	$(CMP) ${CXXFLAGS} -c glossa_tests.cpp -o glossa_tests.o
	$(LNK) ${CXXFLAGS} glossa_tests.o -o $@ ${LDFLAGS}

config_reader_test: config_reader_test.cpp readparameters.o version.o ../../readparameters.h ../../glossa.hpp
	$(CMP) ${CXXFLAGS} -c config_reader_test.cpp -o config_reader_test.o
	${LNK} ${CXXFLAGS} config_reader_test.o $(OBJS) -o $@ ${LDFLAGS}

test: $(TESTS)
	./glossa_tests
	OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 ./config_reader_test

clean:
	rm -f *.o ${TESTS}

.PHONY: default test clean
