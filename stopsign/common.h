#ifndef COMMON_H
#define COMMON_H

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

#endif
