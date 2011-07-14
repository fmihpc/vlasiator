
DEPS_ASCIIBUILDER = ${INSTALL}/gridbuilder.h ${INSTALL}/gridbuilder.cpp asciireaders.h mpibuilder.h mpibuilder.cpp asciibuilder.h asciibuilder.cpp
DEPS_MPIBUILDER = ${INSTALL}/gridbuilder.h ${INSTALL}/gridbuilder.cpp ${INSTALL}/parameters.h mpibuilder.h mpibuilder.cpp
DEPS_RECTCUBOID = ${INSTALL}/gridbuilder.h ${INSTALL}/gridbuilder.cpp ${INSTALL}/parameters.h mpibuilder.h mpibuilder.cpp rect_cuboid_builder.h rect_cuboid_builder.cpp
DEPS_RESTARTBUILDER = ../gridbuilder.h ../gridbuilder.cpp ../parameters.h restartbuilder.h restartbuilder.cpp

asciibuilder.o: ${DEPS_ASCIIBUILDER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c asciibuilder.cpp
	ln -f -s ${INSTALL}/gridbuilders/asciibuilder.o ${INSTALL}/asciibuilder.o

mpibuilder.o: ${DEPS_MPIBUILDER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c mpibuilder.cpp
	ln -f -s ${INSTALL}/gridbuilders/mpibuilder.o ${INSTALL}/mpibuilder.o

rect_cuboid_builder.o: ${DEPS_RECTCUBOID}
	${CMP} ${CXXFLAGS} ${FLAGS} -c rect_cuboid_builder.cpp
	ln -f -s ${INSTALL}/gridbuilders/rect_cuboid_builder.o ${INSTALL}/rect_cuboid_builder.o

restartbuilder.o: ${DEPS_RESTARTBUILDER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c restartbuilder.cpp
	ln -f -s ${INSTALL}/gridbuilders/restartbuilder.o ${INSTALL}/restartbuilder.o

clean:
	rm -rf *.o *.a *~

