cmake_minimum_required(VERSION 3.6)
project(stopsign)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(ROOTSYS /home/jsvirzi/root)
include_directories("${ROOTSYS}/include")

# -L/home/jsvirzi/root/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread
# -pthread -lm -ldl -rdynamic


set(SOURCE_FILES stopsign.cpp)
add_executable(stopsign ${SOURCE_FILES})
target_link_libraries(stopsign ${ROOTSYS}/lib/libCore.so;${ROOTSYS}/lib/libCint.so;${ROOTSYS}/lib/libRIO.so;${ROOTSYS}/lib/libNet.so;${ROOTSYS}/lib/libHist.so;${ROOTSYS}/lib/libGraf.so;${ROOTSYS}/lib/libGraf3d.so;${ROOTSYS}/lib/libGpad.so;${ROOTSYS}/lib/libTree.so;${ROOTSYS}/lib/libRint.so;${ROOTSYS}/lib/libPostscript.so;${ROOTSYS}/lib/libMatrix.so;
${ROOTSYS}/lib/libPhysics.so;${ROOTSYS}/lib/libMathCore.so;${ROOTSYS}/lib/libThread.so)
