#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <TMatrixD.h>

typedef struct {
    double x, y, z;
} Vector3;

typedef struct {
    double x1, y1, z1, x2, y2, z2;
} Interval;

enum {
    X_AXIS = 1,
    Y_AXIS = 2,
    Z_AXIS = 3
};

double getLine(Vector3 *points, int nPoints, int primaryAxis, double lineParameters[4], double endpoints[2]);
double getPlane(Vector3 *points, int nPoints, double *theta);
double solveLinearSystem(TMatrixD &matrix, double *b, double *x);
void projectPointOntoPlane(Vector3 &n, double d, Vector3 &p0, Vector3 &p);

#endif
