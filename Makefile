PROJECTROOT=.
INCLUDES=-I${PROJECTROOT}/include -I${PROJECTROOT}/lidar -I${PROJECTROOT}/dataAcquisition -I${PROJECTROOT}/utils
CFLAGS=-std=c++11 -pthread

all: bin/demo

bin/demo: demo/main.cpp lidar/lidar.cpp dataAcquisition/dataAcquisition.cpp utils/utils.cpp
	mkdir -p bin
	g++ ${CFLAGS} ${INCLUDES} demo/main.cpp lidar/lidar.cpp dataAcquisition/dataAcquisition.cpp utils/utils.cpp -o bin/demo
