
DEPS_ASCIIREADER = ${INSTALL}/gridbuilder.h asciireaders.h asciibuilder.h asciibuilder.cpp
DEPS_RECTCUBOID = ${INSTALL}/gridbuilder.h ${INSTALL}/parameters.h rect_cuboid_builder.h rect_cuboid_builder.cpp

asciibuilder: ${DEPS_ASCIIREADER}
	${CMP} ${CXXFLAGS} ${FLAGS} -c asciibuilder.cpp
	ln -f -s ${INSTALL}/gridbuilders/asciibuilder.o ${INSTALL}/builder.o

rect_cuboid_builder: ${DEPS_RECTCUBOID}
	${CMP} ${CXXFLAGS} ${FLAGS} -c rect_cuboid_builder.cpp
	ln -f -s ${INSTALL}/gridbuilders/rect_cuboid_builder.o ${INSTALL}/builder.o

clean:
	rm -rf *.o *~

