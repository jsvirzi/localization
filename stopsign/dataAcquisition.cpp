#include <signal.h>
#include <fcntl.h>
#include "geometry.h"
#include "lidar.h"
#include "dataAcquisition.h"
#include "common.h"

static int *run = 0;
static void stop(int sig) {
    if (run != 0) {
        *run = 0;
    }
    printf("lidar module: SIGINT received. shutting down application\n");
    // syslog(LOG_NOTICE, "imu: sig %d received", sig);
}

/* this will run in its own thread */
void *dataAcquisitionLoop(void *ptr) {
    DataAcquisitionParams *dataAcquisitionParams = (DataAcquisitionParams *)ptr;
    ThreadParams *threadParams = &dataAcquisitionParams->threadParams;
    run = &threadParams->run;

/* exception handling */
    signal(SIGINT, stop);
    signal(SIGTERM, stop);

    int fd = open(dataAcquisitionParams->filename, O_CREAT);

    threadParams->threadStarted = 1;

    while(threadParams->run == 1) {
//        if(dataAcquisitionParams->lidarPacketTail != )
        if(threadParams->loopWait != 0) {
            usleep(threadParams->loopWait);
        } else {
            usleep(10000); /* default pacing for loop */
        }
    }

    threadParams->threadStarted = 0;

}