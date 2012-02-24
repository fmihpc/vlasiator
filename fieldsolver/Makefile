DEPS_LDZ = ../common.h ../parameters.h ../spatial_cell.hpp ../fieldsolver.h ../priorityqueue.h ../transferstencil.h londrillo_delzanna.cpp

MOVER = ldz

clean:
	rm -rf *.a *.o *~

ldz: ${DEPS_LDZ}
	${CMP} ${CXXFLAGS} ${FLAGS} -c londrillo_delzanna.cpp
	${AR} r libfieldsolver.a londrillo_delzanna.o

libfieldsolver.a: ${DEPS_LDZ}
	make ${MOVER}
