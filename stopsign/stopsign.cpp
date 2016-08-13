#include "math.h"

#include <TROOT.h>
#include <TMatrixD.h>
#include <TApplication.h>
#include <TRandom3.h>

double root2 = 1.414;

typedef struct {
	double x, y, z;
} Vector3;

typedef struct {
	double x1, y1, z1, x2, y2, z2;
} Interval;

double getPlane(Vector3 *points, int nPoints, double *theta);

int main(int argc, char **argv) {
	printf("hello world\n");
	double theta0[3];
	theta0[0] = 1.0;
	theta0[1] = 2.0;
	theta0[2] = 3.0;
	int nPoints = 1024;
	Vector3 *points = new Vector3 [ nPoints ];
	TRandom3 random;
	double a = 1.0, b = 1.5, c = 2.0, d = 1.0;
	for(int i=0;i<nPoints;++i) {
		double x = random.Uniform(1.0);
		double y = random.Uniform(1.0);
		double z = (-d - a * x - b * y) / c;
		points[i].x = x;
		points[i].y = y;
		points[i].z = z;
	}
	double theta[3];
	double det = getPlane(points, nPoints, theta);
	printf("det = %f\n", det);
}

double getPlane(Vector3 *points, int nPoints, double *theta) {
	double sumXX = 0.0, sumXY = 0.0, sumXZ = 0.0, sumYY = 0.0, sumYZ = 0.0, sumZZ = 0.0;
	double sumX = 0.0, sumY = 0.0, sumZ = 0.0;
	for(int i=0;i<nPoints;++i) {
		double x = points[i].x;
		double y = points[i].y;
		double z = points[i].z;
		sumX += x;
		sumY += y;
		sumZ += z;
		sumXX += x * x;
		sumXY += x * y;
		sumXZ += x * z;
		sumYY += y * y;
		sumYZ += y * z;
		sumZZ += z * z;
	}
	TMatrixD matrix(3, 3);
	matrix[0][0] = sumXX;
	matrix[0][1] = sumXY;
	matrix[0][2] = sumXZ;
	matrix[1][0] = sumXY;
	matrix[1][1] = sumYY;
	matrix[1][2] = sumYZ;
	matrix[2][0] = sumXZ;
	matrix[2][1] = sumYZ;
	matrix[2][2] = sumZZ;
	double n[3];
	n[0] = -sumX;
	n[1] = -sumY;
	n[2] = -sumZ;
	double det0 = matrix.Determinant();
	for(int i=0;i<3;++i) {
		TMatrixD M(3, 3);
		for(int j=0;j<3;++j) { for(int k=0;k<3;++k) { M[j][k] = matrix[j][k]; } }
		for(int j=0;j<3;++j) { M[j][i] = n[j]; }
		double det = M.Determinant();
		theta[i] = det / det0;
	}
	return det0;
}

int crossProduct(Vector3 *u, Vector3 *v, Vector3 *n, int normalize) {
	if(normalize == 0) {
		n->x = u->y * v->z - v->y * u->z;
		n->y = u->z * v->x - v->z * u->x;
		n->z = u->x * v->y - u->y * v->x;
	} else {
		double nx = u->y * v->z - v->y * u->z;
		double ny = u->z * v->x - v->z * u->x;
		double nz = u->x * v->y - u->y * v->x;
		double mag = sqrt(nx * nx + ny * ny + nz * nz);
		if(mag == 0.0) return 1;
		double invmag = 1.0 / mag;
		n->x = nx * invmag;
		n->y = ny * invmag;
		n->z = nz * invmag;
	}
	return 0;
}

int getNormalToPlane(Vector3 *vec1, Vector3 *vec2, Vector3 *vec3, Vector3 *n) {
	Vector3 u;
	u.x = vec2->x - vec1->x;
	u.y = vec2->y - vec1->y;
	u.z = vec2->z - vec1->z;
	Vector3 v;
	v.x = vec3->x - vec1->x;
	v.y = vec3->y - vec1->y;
	v.z = vec3->z - vec1->z;
	crossProduct(&u, &v, n, 1);
	// double nx = u.y * v.z - u.z * v.y;
	// double ny = v.x * u.z - u.x * v.z;
	// double nz = u.x * v.z - v.x * u.y;
	// double mag = sqrt(nx * nx + ny * ny + nz * nz);
	// if(mag == 0.0) { return -1; }
	// double invmag = 1.0 / mag;
	// n->x = nx * invmag;
	// n->y = ny * invmag;
	// n->z = nz * invmag;
	return 0;
}

double intervalLength(Interval *interval) {
	double ax = interval->x2 - interval->x1;
	double ay = interval->y2 - interval->y1;
	double az = interval->z2 - interval->z1;
	return sqrt(ax * ax + ay * ay + az * az);
}

double vectorLength(Vector3 *vector) {
	double ax = vector->x;
	double ay = vector->y;
	double az = vector->z;
	return sqrt(ax * ax + ay * ay + az * az);
} 

int getPlane(Interval *interval, int nIntervals, Vector3 *p, Vector3 *n) {
	double scale;
	double z1, z2, l1, l2;
	double x0 = 1.0 / (1.0 + root2);
	double z0 = (2.0 + root2) * x0;
	Vector3 u, v;
	if(nIntervals == 3) {

		u.x = interval[0].x2 - interval[0].x1;
		u.y = interval[0].y2 - interval[0].y1;
		u.z = interval[0].z2 - interval[0].z1;
		v.x = interval[1].x1 - interval[0].x1;
		v.y = interval[1].y1 - interval[0].y1;
		v.z = interval[1].z1 - interval[0].z1;
		Vector3 n1;
		crossProduct(&u, &v, &n1, 0);

		scale = intervalLength(&interval[1]);
		l1 = intervalLength(&interval[0]);
		l2 = intervalLength(&interval[2]);
		z1 = z0 - 0.5 * l1; 
		z2 = z0 - 0.5 * l2;
//		estimateCenter(&interval[0],
	}
}

int estimateCenter(Interval *interval, Vector3 *planeNormal, double length, Vector3 *center) {
	Vector3 v;
	v.x = interval->x2 - interval->x1;
	v.y = interval->y2 - interval->y1;
	v.z = interval->z2 - interval->z1;
	Vector3 n;
	crossProduct(planeNormal, &v, &n, 1); 
	center->x = 0.5 * (interval->x1 + interval->x2) + n.x * length;
	center->y = 0.5 * (interval->y1 + interval->y2) + n.y * length;
	center->z = 0.5 * (interval->z1 + interval->z2) + n.z * length;
}

