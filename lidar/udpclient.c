/* 
 * udpclient.c - A simple UDP client
 * usage: udpclient <host> <port>
 */
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

/* 
 * error - wrapper for perror
 */
void error(char *msg) {
    perror(msg);
    exit(0);
}

#if 0

int main(int argc, char **argv) {
    int sockfd, portno, n;
    int serverlen;
    struct sockaddr_in serveraddr;
    struct hostent *server;
    char *hostname;
    char buf[BUFSIZE];

    /* check command line arguments */
//    if (argc != 3) {
//       fprintf(stderr,"usage: %s <hostname> <port>\n", argv[0]);
//       exit(0);
//    }
//    hostname = argv[1];
//    portno = atoi(argv[2]);

    portno = 2368;

    /* socket: create the socket */
    sockfd = socket(PF_INET, SOCK_DGRAM, 0);
    if (sockfd < 0) 
        error("ERROR opening socket");

    /* gethostbyname: get the server's DNS entry */
    // server = gethostbyname(hostname);
    // if (server == NULL) {
        // fprintf(stderr,"ERROR, no such host as %s\n", hostname);
        // exit(0);
    // }

    /* build the server's Internet address */
    bzero((char *) &serveraddr, sizeof(serveraddr));
    serveraddr.sin_family = AF_INET;
    // serveraddr.sin_port = htons(portno);
    // serveraddr.sin_addr.s_addr = inet_addr("192.168.1.201");
    // serveraddr.sin_addr.s_addr = inet_addr("192.168.1.77");
    serveraddr.sin_addr.s_addr = INADDR_ANY;
    serveraddr.sin_port = htons(2368);
    // memset(serveraddr.sin_zero, '\0', sizeof serveraddr.sin_zero);  
    // bcopy((char *)server->h_addr, (char *)&serveraddr.sin_addr.s_addr, server->h_length);

    if(bind(sockfd, (struct sockaddr *)&serveraddr, sizeof(serveraddr)) == -1) {
	perror("bind");
	return -1;
    }

    /* get a message from the user */
    bzero(buf, BUFSIZE);
    // printf("Please enter msg: ");
    // fgets(buf, BUFSIZE, stdin);

    /* send the message to the server */
    // serverlen = sizeof(serveraddr);
    // n = sendto(sockfd, buf, strlen(buf), 0, &serveraddr, serverlen);
    // if (n < 0) 
      // error("ERROR in sendto");
    
    /* print the server's reply */
    while(1) {
        printf("recvfrom() in\n");
    // n = recvfrom(sockfd, buf, BUFSIZE, 0, &serveraddr, &serverlen);
    n = recvfrom(sockfd, buf, 1, 0, &serveraddr, &serverlen);
        printf("recvfrom() out\n");
    buf[serverlen] = 0;
    if (n < 0) 
      error("ERROR in recvfrom");
    printf("Echo from server: %s. length = %d", buf, serverlen);
    }
    return 0;
}

#endif

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
    int i;
    for(i=0;i<12;++i) {
        printf("flag(%d) = 0x%x\n", i, packet->dataBlock[i].flag); 
    }
}

int main() {

    int i;
    int port = 2368;
    int sockfd_ = socket(PF_INET, SOCK_DGRAM, 0);
    if (sockfd_ == -1) { perror("socket"); return; }

    printf("sizeof(uint16) = %ld. sizeof(uint8) = %ld\n", sizeof(uint16), sizeof(uint8));
    printf("struct sizes = %ld %ld %ld\n", sizeof(LidarChannelDatum), sizeof(LidarDataBlock), sizeof(LidarPacket));
    getchar();
  
    struct sockaddr_in my_addr;                     // my address information
    memset(&my_addr, 0, sizeof(my_addr));    // initialize to zeros
    my_addr.sin_family = AF_INET;            // host byte order
    my_addr.sin_port = htons(port);          // port in network byte order
    my_addr.sin_addr.s_addr = INADDR_ANY;    // automatically fill in my IP
  
    if (bind(sockfd_, (struct sockaddr *)&my_addr, sizeof(struct sockaddr)) == -1) { perror("bind"); }
  
    if (fcntl(sockfd_,F_SETFL, O_NONBLOCK|FASYNC) < 0) { perror("non-block"); return; }

    while(getchar) {
        int retval, pollTimeout = 1000;
        struct pollfd fds;
        //do {
          //  retval = poll(&fds, 1, pollTimeout);
        //} while((fds.revents & POLLIN) == 0);

        struct sockaddr_in sender_address;
        int sender_address_len = sizeof(sender_address);
        int nbytes = recvfrom(sockfd_, buf, BUFSIZE, 0, (struct sockaddr *)&sender_address, &sender_address_len);
        printf("nbytes = %d %ld\n", nbytes, sizeof(LidarPacket));

        if(nbytes == sizeof(LidarPacket)) {
            LidarPacket *packet = (LidarPacket *)buf;
            dumpLidarPacket(packet);
        }

//        for(i=0;i<nbytes;++i) {
//            printf("%c", buf[i]);
//            if(i && ((i % 32) == 0)) printf("\n");
//        }
    }


}
