#include <iostream>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <netdb.h>

#include <arpa/inet.h>
#include <poll.h>
#include <errno.h>
#include <fcntl.h>
#include <sys/file.h>
#include <math.h>
#include <csignal>

#include "geometry.h"
#include "lidar.h"
#include "common.h"
#include "dataAcquisition.h"

static int run = 1;

static void stop(int sig) {
    run = 0;
    printf("SIGINT received. shutting down application\n");
    // syslog(LOG_NOTICE, "imu: sig %d received", sig);
}

void test();

const int nLidarPackets = 1024;
//LidarPacket lidarPackets[NLIDARPACKETS];
//int lidarPacketBufferHead = 0;
//int lidarPacketBufferSize = NLIDARPACKETS;
//void *lidarLloop(void *p);

enum {
    LidarThreadId = 0,
    DataAcquisitionThreadId = 1,
    NThreads
};
pthread_t tid[NThreads];

int waitForThreadToStart(ThreadParams *threadParams) {
    while(run && threadParams->threadStarted == 0) {
        printf("waiting for %s service to start\n", threadParams->threadName);
        if(threadParams->errorCode != 0) {
            printf("error code %d encountered starting %s service", threadParams->errorCode, threadParams->threadName);
            break;
        }
    }

    if(threadParams->errorCode == 0) {
        printf("%s service started\n", threadParams->threadName);
    }

    return threadParams->errorCode;
}

int waitForThreadToStop(ThreadParams *threadParams) {
    char *threadResult;
    pthread_join(threadParams->tid, (void **)&threadResult);
}

int main(int argc, char **argv) {

    int i, err;
    std::string ofile = "dataAcquisition.jsv";

    for(i=1;i<argc;++i) {
        if(strcmp(argv[i], "-o") == 0) {
            ofile = argv[++i];
        }
    }

    /* important preliminaries */
    if ((sizeof(uint8) != 1) ||
        (sizeof(uint16) != 2) ||
        (sizeof(LidarPacket) != lidarPacketSize) ||
        (sizeof(LidarDataBlock) != lidarBlockSize) ||
        (sizeof(LidarChannelDatum) != lidarDatumSize)) {
        printf("sizeof(uint16) = %ld. sizeof(uint8) = %ld\n", sizeof(uint16), sizeof(uint8));
        printf("struct sizes = %ld %ld %ld\n", sizeof(LidarChannelDatum), sizeof(LidarDataBlock), sizeof(LidarPacket));
        return -1;
    }

    // test();

    signal(SIGINT, stop); /* exception handling */
    signal(SIGTERM, stop); /* exception handling */

    ThreadParams *threadParams;

    /* get the lidar thread up and running */
    LidarServerParams lidarServerParams;
    memset(&lidarServerParams, 0, sizeof(lidarServerParams));
    threadParams = &lidarServerParams.threadParams;
    threadParams->threadStarted = 0;
    threadParams->run = 1;
    threadParams->errorCode = 0;
    threadParams->threadName = "lidar";
    lidarServerParams.packetBufferSize = nLidarPackets;
    lidarServerParams.packetBuffer = new LidarPacket[nLidarPackets];
    lidarServerParams.packetBufferHead = 0;
    err = pthread_create(&tid[LidarThreadId], NULL, lidarLoop, &lidarServerParams);
    threadParams->tid = tid[LidarThreadId];

    if(waitForThreadToStart(threadParams) != 0) { printf("problems in paradise"); exit(1); }

    /* data acquisition loop up and running */
    DataAcquisitionParams dataAcquisitionParams;
    memset(&dataAcquisitionParams, 0, sizeof(dataAcquisitionParams));
    threadParams = &dataAcquisitionParams.threadParams;
    threadParams->threadStarted = 0;
    threadParams->run = 1;
    threadParams->errorCode = 0;
    threadParams->threadName = "data acquisition";
    dataAcquisitionParams.lidarPacketBuffer = lidarServerParams.packetBuffer;
    dataAcquisitionParams.lidarPacketTail = 0;
    dataAcquisitionParams.lidarPacketSize = nLidarPackets;
    err = pthread_create(&tid[DataAcquisitionThreadId], NULL, dataAcquisitionLoop, &dataAcquisitionParams);
    threadParams->tid = tid[DataAcquisitionThreadId];

    if(waitForThreadToStart(threadParams) != 0) { printf("problems in paradise"); exit(1); }

    while(run) {
        sleep(1);
        printf("running...\n");
    }

    waitForThreadToStop(&lidarServerParams.threadParams);
    waitForThreadToStop(&dataAcquisitionParams.threadParams);
    printf("goodbye\n");

}

void clearHistogram(int *histogram, int nX) {
    memset(histogram, 0, sizeof(int) * nX);
}

void fillHistogram(int *histogram, int x) {
    ++histogram[x];
}

void fillHistogram(int *histogram, int nX, int x, int y) {
    ++histogram[y * nX + x];
}

void analyzePointCloud(LidarData *data, int nPoints) {
    int iChannel, iPoint;
    LidarData *datum = data;
    for(iChannel=0;iChannel<nChannels;++iChannel) {
        clearHistogram(histogramAzimuth[iChannel], 256);
        clearHistogram(histogramIntensity[iChannel], 256);
        clearHistogram(histogramIntensityAzimuth[iChannel], 256 * 256);
    }
    for(iPoint=0;iPoint<nPoints;++iPoint,++datum) {
        int intensity = datum->intensity & 0xff;
        int azimuth = ((int)(datum->phi * invTwoPi * 255.0)) & 0xff;
        int channel = datum->channel;
        fillHistogram(histogramIntensity[channel], intensity);
        fillHistogram(histogramAzimuth[channel], azimuth);
        fillHistogram(histogramIntensityAzimuth[channel], 256, azimuth, intensity);
    }
}
