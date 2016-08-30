#ifndef DATAACQUISITION_H
#define DATAACQUISITION_H

#include "common.h"
#include "lidar.h"
#include "dataAcquisition.h"

typedef struct {
    ThreadParams threadParams;
    const char *filename;
    LidarPacket *lidarPacketBuffer;
    int lidarPacketTail;
    int lidarPacketSize;
    LidarServerParams *lidarServerParams;
} DataAcquisitionParams;

void *dataAcquisitionLoop(void *ptr);

typedef struct {
    long unsigned int syncWord;
    long int structId;
    long int structLength;
    long int time;
} DataAcquisitionPreamble;

#endif
