#ifndef LIDAR_H
#define LIDAR_H

#include "common.h"

#define NBLOCKS 12

typedef struct {
    uint8 distance[2];
    uint8 reflectivity;
} LidarChannelDatum;

typedef struct {
    uint16 flag;
    uint8 azimuthLo, azimuthHi; // these need rearrangement
    LidarChannelDatum data[32];
} LidarDataBlock;

typedef struct {
//    uint8 header[42];
    LidarDataBlock dataBlock[NBLOCKS];
    uint8 timestamp[4];
    uint8 factoryField[2];
} LidarPacket;

typedef struct {
    double R, theta, phi; // R = distance, phi = azimuth, theta = altitude
    int intensity, channel;
} LidarData;

enum {
    StrongestReturn = 0x37,
    LastReturn = 0x38,
    DualReturn = 0x39
} ReturnModes;

double getAzimuth(LidarDataBlock *block);

const int nChannels = 16;
const int lidarPacketSize = 1206;
const int lidarBlockSize = 100;
const int lidarDatumSize = 3;

typedef struct {
    ThreadParams threadParams;
    LidarPacket *packetBuffer;
    int packetBufferHead;
    int packetBufferSize;
    const char *ipAddress;
} LidarServerParams;

void *lidarLoop(void *ptr);

#endif
