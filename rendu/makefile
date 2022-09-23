
#########################
# Options du compilateur
#########################

CC = g++

########################
# Détection de l'OS
#########################

OS = $(shell uname)


C++FLAGS = -I. -I /usr/include \
               -I /usr/local/include \
               -I /opt/local/include/GL \
               -Wno-deprecated

ifeq (${OS},Darwin)
  C++FLAGS += -I /System/Library/Frameworks/OpenGL.framework/Headers
else
  C++FLAGS += -I /usr/include/GL -I /usr/include
endif

CCOPTIONS = -g -O0 -Wall

#########################
# Librairies
#########################

ifeq ($(OS),Darwin)
	LIBS = -lm -framework OpenGL -framework GLUT
else
	LIBS = -lm -lGLU -lGL -lglut
endif


EXEC = particleView


#--------------------------------------------------------------------------
# Compilation

all: $(EXEC)

SRC_C = ${wildcard *.cpp}

OBJS = $(patsubst %.cpp,%.o,${SRC_C})

.cpp.o:
	${CC} -c $*.cpp ${C++FLAGS} $(CCOPTIONS)

$(EXEC) : ${OBJS}
	${CC} -o ${EXEC} ${OBJS} ${LD_FLAGS} ${LIBS}

# Nettoyer
clean:
	rm *.o ${EXEC} *~

mem :
	valgrind --leak-check=yes $(EXEC) -I 800 800

toto:
	echo ${OS}
	echo ${OBJS}
