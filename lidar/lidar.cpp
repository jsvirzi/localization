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
#include <time.h>

#include "lidar.h"

static int *run = 0;
static void stop(int sig) {
    if (run != 0) {
        *run = 0;
    }
    printf("lidar module: SIGINT received. shutting down application\n");
    // syslog(LOG_NOTICE, "imu: sig %d received", sig);
}

const double cycleTimeBetweenFirings = 2.304; // microseconds
const double rechargePeriod = 18.43; // microseconds
const int distanceUnit = 2.0; // millimeters
const int VLP16Device = 0x22;
const double laserPolarAngle[nChannels] = {
        -15.0 * M_PI / 180.0,
        1.0 * M_PI / 180.0,
        -13.0 * M_PI / 180.0,
        3.0 * M_PI / 180.0,
        -11.0 * M_PI / 180.0,
        5.0 * M_PI / 180.0,
        -9.0 * M_PI / 180.0,
        7.0 * M_PI / 180.0,
        -7.0 * M_PI / 180.0,
        9.0 * M_PI / 180.0,
        -5.0 * M_PI / 180.0,
        11.0 * M_PI / 180.0,
        -3.0 * M_PI / 180.0,
        13.0 * M_PI / 180.0,
        -1.0 * M_PI / 180.0,
        15.0 * M_PI / 180.0
};

int histogramIntensity[nChannels][256];
int histogramAzimuth[nChannels][256];
int histogramIntensityAzimuth[nChannels][256 * 256];
void analyzePointCloud(LidarData *data, int nPoints);

#define NPOINTS (3600 * 16)
LidarData pointCloud[NPOINTS];
int pointCloudIndex;

int convertLidarPacketToLidarData(LidarPacket *lidarPacket, LidarData *lidarData) {
    int iBlock, iChannel, nPoints = 0;
    double azimuth[NBLOCKS];
    double interpolatedDeltaAzimuth = 0.0;
    for(iBlock=1;iBlock<NBLOCKS;++iBlock) {
        LidarDataBlock *lidarDataBlock = &lidarPacket->dataBlock[iBlock];
        azimuth[iBlock] = getAzimuth(lidarDataBlock);
        double a = azimuth[iBlock] - azimuth[iBlock-1];
        interpolatedDeltaAzimuth += (a * a);
    }
    interpolatedDeltaAzimuth = 0.5 * sqrt(interpolatedDeltaAzimuth) / (NBLOCKS - 1);

    for(iBlock=0;iBlock<NBLOCKS;++iBlock) {
        LidarDataBlock *lidarDataBlock = &lidarPacket->dataBlock[iBlock];
        double azimuth0 = azimuth[iBlock];
        for(iChannel=0;iChannel<nChannels;++iChannel) {
            LidarChannelDatum *datum = &lidarDataBlock->data[iChannel];
            double a = distanceUnit * datum->distance[0];
            lidarData->R = 256.0 * a + distanceUnit * datum->distance[1];
            lidarData->phi = azimuth0 + iChannel * cycleTimeBetweenFirings;
            lidarData->theta = laserPolarAngle[iChannel];
            lidarData->intensity = datum->reflectivity;
            lidarData->channel = iChannel;
            ++lidarData;
            ++nPoints;
        }
        azimuth0 += interpolatedDeltaAzimuth;
        for(iChannel=0;iChannel<nChannels;++iChannel) {
            LidarChannelDatum *datum = &lidarDataBlock->data[nChannels + iChannel];
            double a = distanceUnit * datum->distance[0];
            lidarData->R = 256.0 * a + distanceUnit * datum->distance[1];
            lidarData->phi = azimuth0 + iChannel * cycleTimeBetweenFirings;
            lidarData->theta = laserPolarAngle[iChannel];
            lidarData->intensity = datum->reflectivity;
            lidarData->channel = iChannel;
            ++lidarData;
            ++nPoints;
        }
    }
    return nPoints;
}

int isBeginningPacket(LidarPacket *packet) {
    LidarDataBlock *lidarDataBlock = &packet->dataBlock[0];
    double azimuth = getAzimuth(lidarDataBlock);
//    printf("isBeginningPacket(LidarPacket packet=%p): azimuth = %.5f\n", packet, azimuth);
    double c = cos(azimuth), s = sin(azimuth);
    if ((fabs(c - 1.0) < 0.1) && (fabs(s - 0.0) < 0.1)) return 1;
    else return 0;
}

void dumpLidarPacket(LidarPacket *packet) {
    int iBlock;
    for(iBlock=0;iBlock<NBLOCKS;++iBlock) {
        LidarDataBlock *lidarDataBlock = &packet->dataBlock[iBlock];
        double azimuth = getAzimuth(lidarDataBlock);
        printf("flag(%d) = 0x%x. azimuth = %.3f(rad) = %.3f(deg)\n", iBlock, packet->dataBlock[iBlock].flag, azimuth, azimuth * 180.0 / M_PI);
//        double azimuth = lidarDataBlock->azimuth
//        printf("")
    }
}

double getAzimuth(LidarDataBlock *block) {
    unsigned char azimuthLo = block->azimuthLo;
    unsigned char azimuthHi = block->azimuthHi;
    int azimuthData = azimuthHi;
    azimuthData = (azimuthData << 8) | azimuthLo;
    double azimuth = azimuthData * M_PI / 18000.0;
    return azimuth;
}

/* this will run on its own thread */
void *lidarLoop(void *ptr) {

    LidarServerParams *lidarServerParams = (LidarServerParams *)ptr;
    ThreadParams *threadParams = &lidarServerParams->threadParams;
    run = &threadParams->run;

/* exception handling */
    signal(SIGINT, stop);
    signal(SIGTERM, stop);

    int i;
    int port = 2368;
    int sockfd = socket(PF_INET, SOCK_DGRAM, 0);
    if (sockfd == -1) { perror("socket"); threadParams->errorCode = -1; return 0; }

    struct sockaddr_in my_addr;                     // my address information
    memset(&my_addr, 0, sizeof(my_addr));    // initialize to zeros
    my_addr.sin_family = AF_INET;            // host byte order
    my_addr.sin_port = htons(port);          // port in network byte order
    my_addr.sin_addr.s_addr = INADDR_ANY;    // automatically fill in my IP

    if (bind(sockfd, (struct sockaddr *) &my_addr, sizeof(struct sockaddr)) == -1) { perror("bind"); }

    if (fcntl(sockfd, F_SETFL, O_NONBLOCK | FASYNC) < 0) { perror("non-block"); threadParams->errorCode = -2; return 0; }

    time_t t0 = time(0);
    int nCycles = 0, nPackets = 0, packet0Index = 0;

    threadParams->threadStarted = 1;

    while (threadParams->run == 1) {

//        int retval, pollTimeout = 1000;
//        struct pollfd fds;
//        memset(&fds, 0, sizeof(fds));
//        fds.fd = sockfd;
//        fds.events = POLLIN;
//        do {
//            retval = poll(&fds, 1, pollTimeout);
//            if (retval == 0) {
//                printf("Velodyne timed out\n");
//            }
//            int errorCondition = (fds.revents & POLLERR) || (fds.revents & POLLHUP) || (fds.revents & POLLNVAL);
//
//           if (errorCondition) {
//                printf("error condition exists in UDP packets from Velodyne Lidar\n");
//            }
//            usleep(10000);
//        } while((fds.revents & POLLIN) == 0);

        struct sockaddr_in sender_address;
        socklen_t sender_address_len = sizeof(sender_address);
        int nbytes = recvfrom(sockfd, &lidarServerParams->packetBuffer[lidarServerParams->packetBufferHead], sizeof(LidarPacket), 0, (struct sockaddr *) &sender_address, &sender_address_len);

        // printf("nbytes = %d %ld\n", (int)nbytes, sizeof(LidarPacket));

        if (nbytes < 0) {
            if (errno != EWOULDBLOCK) {
                perror("recvfail");
            }
        }

        if (nbytes == sizeof(LidarPacket)) {
            LidarPacket *packet = &lidarServerParams->packetBuffer[lidarServerParams->packetBufferHead];
            lidarServerParams->packetBufferHead = (lidarServerParams->packetBufferHead + 1) % lidarServerParams->packetBufferSize;
            if(isBeginningPacket(packet) && (nPackets >= (packet0Index + 10))) {
                packet0Index = nPackets;
                printf("found beginning packet. previous index = %d\n", pointCloudIndex);
                analyzePointCloud(pointCloud, pointCloudIndex);
                pointCloudIndex = 0;
//                dumpLidarPacket(packet);
                ++nCycles;
                time_t deltaTime = time(0) - t0;
                if (deltaTime > 0) {
                    printf("frequency = %d/%ld = %ld\n", nCycles, deltaTime, nCycles / deltaTime);
                }
            }
            ++nPackets;
            int n = convertLidarPacketToLidarData(packet, &pointCloud[pointCloudIndex]);
            pointCloudIndex += n;
        }

        if(threadParams->loopWait != 0) {
            usleep(threadParams->loopWait);
        }

    }

    threadParams->threadStarted = 0;

    return 0;
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