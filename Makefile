PROJECTROOT=.
INCLUDES=-I${PROJECTROOT}/include -I${PROJECTROOT}/lidar -I${PROJECTROOT}/dataAcquisition
CFLAGS=-std=c++11 -pthread

all: bin/demo

bin/demo: demo/main.cpp lidar/lidar.cpp dataAcquisition/dataAcquisition.cpp
	mkdir -p bin
	g++ ${CFLAGS} ${INCLUDES} demo/main.cpp lidar/lidar.cpp dataAcquisition/dataAcquisition.cpp -o bin/demo
