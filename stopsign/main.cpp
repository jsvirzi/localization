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

#define BUFSIZE 1300
unsigned char buf [ BUFSIZE ];

typedef unsigned short int uint16;
typedef unsigned char uint8;

typedef struct {
    uint8 distance[2];
    uint8 reflectivity;
} LidarChannelDatum;

typedef struct {
    uint16 flag;
    uint16 azimuth;
    LidarChannelDatum data[32];
} LidarDataBlock;

typedef struct {
//    uint8 header[42];
    LidarDataBlock dataBlock[12];
    uint8 timestamp[4];
    uint8 factoryField[2];
} LidarPacket;

void dumpLidarPacket(LidarPacket *packet) {
    int iBlock;
    for(iBlock=0;iBlock<12;++iBlock) {
        LidarDataBlock *lidarDataBlock = &packet->dataBlock[iBlock];
        printf("flag(%d) = 0x%x. azimuth = 0x%x\n", iBlock, packet->dataBlock[iBlock].flag, lidarDataBlock->azimuth);
//        double azimuth = lidarDataBlock->azimuth
//        printf("")
    }
}

int main() {

    int i;
    int port = 2368;
    int sockfd = socket(PF_INET, SOCK_DGRAM, 0);
    if (sockfd == -1) { perror("socket"); return -1; }

    printf("sizeof(uint16) = %ld. sizeof(uint8) = %ld\n", sizeof(uint16), sizeof(uint8));
    printf("struct sizes = %ld %ld %ld\n", sizeof(LidarChannelDatum), sizeof(LidarDataBlock), sizeof(LidarPacket));
    getchar();

    struct sockaddr_in my_addr;                     // my address information
    memset(&my_addr, 0, sizeof(my_addr));    // initialize to zeros
    my_addr.sin_family = AF_INET;            // host byte order
    my_addr.sin_port = htons(port);          // port in network byte order
    my_addr.sin_addr.s_addr = INADDR_ANY;    // automatically fill in my IP

    if (bind(sockfd, (struct sockaddr *) &my_addr, sizeof(struct sockaddr)) == -1) { perror("bind"); }

    if (fcntl(sockfd, F_SETFL, O_NONBLOCK | FASYNC) < 0) { perror("non-block"); return -1; }

    while (getchar) {
        int retval, pollTimeout = 1000;
        struct pollfd fds;
        memset(&fds, 0, sizeof(fds));
        fds.fd = sockfd;
        fds.events = POLLIN;
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
        int nbytes = recvfrom(sockfd, buf, BUFSIZE, 0, (struct sockaddr *) &sender_address, &sender_address_len);

        // printf("nbytes = %d %ld\n", (int)nbytes, sizeof(LidarPacket));

        if (nbytes < 0) {
            if (errno != EWOULDBLOCK) {
                perror("recvfail");
            }
        }

        if (nbytes == sizeof(LidarPacket)) {
            LidarPacket *packet = (LidarPacket *) buf;
            dumpLidarPacket(packet);
        }

    }
}