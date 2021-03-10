#CXX_OPTIONS = -O3 -DCHECK_DENSITY -DDEBUG_SOLVER -std=c++0x -W -Wall -Wextra -pedantic -Wno-missing-braces
#CXX_OPTIONS = -O3 -DCHECK_DENSITY -W -Wall -Wextra -pedantic -Wno-missing-braces
#CXX_OPTIONS = -O3 -DDEBUG_GUMICS -W -Wall -Wextra -pedantic -Wno-missing-braces
CXX_OPTIONS = -O3 -W -Wall -Wextra -pedantic -Wno-missing-braces -std=c++0x
#CXX_OPTIONS = -g -DDEBUG -W -Wall -Wextra -pedantic -Wno-missing-braces

# Uncomment one of the following:
include ../../MAKE/Makefile.${VLASIATOR_ARCH}

FLAGS = \
	$(INC_BOOST) \
	$(LIB_BOOST) \
	-lboost_program_options -lm

HEADERS = \
	../B0.hpp \
	../ode.hpp \
	../quadr.hpp

SOURCES = \
	../B0.cpp \
	../ode.cpp \
	../quadr.cpp

all: test1

test1: test1.cpp $(SOURCES) $(HEADERS) Makefile
	$(CMP) $(CXX_OPTIONS) $(SOURCES) test1.cpp $(FLAGS) -o test1

c: clean
clean:
	rm -f test1

