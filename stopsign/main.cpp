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

static int run = 1;

static void stop(int sig) {
    run = 0;
    printf("SIGINT received. shutting down application\n");
    // syslog(LOG_NOTICE, "imu: sig %d received", sig);
}

void test();

void *lidar_loop(void *p);

enum {
    LidarThreadId = 0,
    NThreads
};
pthread_t tid[NThreads];

int main() {

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

    LidarServerParams lidarServerParams;
    memset(&lidarServerParams, 0, sizeof(lidarServerParams));
    lidarServerParams.threadStarted = 0;
    lidarServerParams.run = 1;
    int err = pthread_create(&tid[LidarThreadId], NULL, lidar_loop, &lidarServerParams);

    while(run && lidarServerParams.threadStarted == 0) {
        printf("waiting for lidar service to start\n");
        if(lidarServerParams.errorCode != 0) {
            printf("error code %d encountered starting lidar service", lidarServerParams.errorCode);
            break;
        }
    }
    if(lidarServerParams.errorCode == 0) {
        printf("lidar service started\n");
    }

    while(run) {
        sleep(1);
        printf("running...\n");
    }

    char *thread_result;
    pthread_join(tid[LidarThreadId], (void**) &thread_result);
    printf("goodbye\n");

}

