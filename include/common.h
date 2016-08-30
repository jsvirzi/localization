#ifndef COMMON_H
#define COMMON_H

#include <math.h>

typedef unsigned short int uint16;
typedef unsigned char uint8;

typedef struct {
    int threadStarted;
    int run;
    int errorCode;
    const char *threadName;
    pthread_t tid;
    int loopWait;
} ThreadParams;

const double invTwoPi = 1.0 / (2.0 * M_PI);

#endif
