#ifndef DATAACQUISITION_H
#define DATAACQUISITION_H

typedef struct {
    ThreadParams threadParams;
    char *filename;
    LidarPacket *lidarPacketBuffer;
    int lidarPacketTail;
    int lidarPacketSize;
} DataAcquisitionParams;

void *dataAcquisitionLoop(void *ptr);

#endif
