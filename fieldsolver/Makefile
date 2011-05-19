DEPS_LDZ = ../parameters.h ../cell_spatial.h ../fieldsolver.h londrillo_delzanna.cpp

MOVER = ldz

clean:
	rm -rf *.a *.o *~

ldz: ${DEPS_LDZ}
	${CMP} ${FLAGS} ${CXXFLAGS} -c londrillo_delzanna.cpp
	ar r libfieldsolver.a londrillo_delzanna.o

libfieldsolver.a: ${DEPS_LDZ}
	make ${MOVER}
