#include <time.h>

unsigned long int elapsedTime() {
    struct timespec timeInfo;
    clock_gettime(CLOCK_REALTIME, &timeInfo);
    long timeValue = timeInfo.tv_sec;
    timeValue = timeValue * 1000000 + timeInfo.tv_nsec;
    return timeValue;
}
