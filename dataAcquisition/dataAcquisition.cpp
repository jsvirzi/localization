#include <signal.h>
#include <fcntl.h>
#include <cstdio>
#include <unistd.h>
#include <time.h>
#include <lidar.h>
#include <sys/stat.h>
#include <common.h>
#include "lidar.h"
#include "dataAcquisition.h"
#include "common.h"
#include "utils.h"

static int *run = 0;
static void stop(int sig) {
    if (run != 0) {
        *run = 0;
    }
    printf("lidar module: SIGINT received. shutting down application\n");
    // syslog(LOG_NOTICE, "imu: sig %d received", sig);
}

const int timeBetweenSyncs = 5; /* flush data every 5 seconds, at least */

/* this will run in its own thread */
void *dataAcquisitionLoop(void *ptr) {
    DataAcquisitionParams *dataAcquisitionParams = (DataAcquisitionParams *)ptr;
    ThreadParams *threadParams = &dataAcquisitionParams->threadParams;
    LidarServerParams *lidarServerParams = dataAcquisitionParams->lidarServerParams;
    run = &threadParams->run;

    DataAcquisitionPreamble preamble;

/* exception handling */
    signal(SIGINT, stop);
    signal(SIGTERM, stop);

    int fd = open(dataAcquisitionParams->filename, O_RDWR | O_CREAT | O_TRUNC, S_IWRITE | S_IREAD);
    if(fd < 0) {
        perror("file open");
        threadParams->errorCode = -1;
        return 0;
    }

    threadParams->threadStarted = 1;

    time_t startTime = time(0);
    time_t nextSyncTime = startTime + timeBetweenSyncs;

    while(threadParams->run == 1) {

        /* any lidar data to write? */
        int size = dataAcquisitionParams->lidarServerParams->packetBufferSize;
        preamble.structId = LidarPackedId;
        preamble.structLength = sizeof(DataAcquisitionPreamble) + sizeof(LidarPacket);
        preamble.syncWord = syncWord;
        preamble.time = elapsedTime();
        while(dataAcquisitionParams->lidarPacketTail != lidarServerParams->packetBufferHead) {
            LidarPacket *lidarPacket = &lidarServerParams->packetBuffer[dataAcquisitionParams->lidarPacketTail];
            write(fd, &preamble, sizeof(DataAcquisitionPreamble));
            write(fd, lidarPacket, sizeof(LidarPacket));
            dataAcquisitionParams->lidarPacketTail = (dataAcquisitionParams->lidarPacketTail + 1) % size;
        }

        /* general bookkeeping */
        time_t now = time(0);
        if(now > nextSyncTime) {
            fsync(fd);
            nextSyncTime += timeBetweenSyncs;
            printf("data sync occurred at %ld", now);
        }

        /* pace thread */
        if(threadParams->loopWait != 0) {
            usleep(threadParams->loopWait);
        } else {
            usleep(10000); /* default pacing for loop */
        }

    }

    close(fd);

    threadParams->threadStarted = 0;

    return 0;
}