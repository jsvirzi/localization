#include "math.h"

#include <TROOT.h>
#include <TMatrixD.h>
#include <TApplication.h>
#include <TRandom3.h>

#include "geometry.h"

double root2 = 1.414;

// int generateGrate(double *m, double *b, double *delta, double *l, int nGrates, Vector3 points, int nPoints);
int generateGrate(Interval *grates, int nGrates, Vector3 *points, int nPoints);
void generateLine(Vector3 *points, int nPoints, int primaryAxis, double lineParameters[4], double endPoints[2], double sigma);

TRandom3 rndm;

void test() {
	/* y is primary axis.
	 * x = m1 * y + b1
	 * z = m2 * y + b2
	 * put the stop sign center at x0 = 5m, y0 = -3m, z0 = 2m
	 * the stop sign is slightly "yawed" by 4 degrees = 0.07 radians
	 * the stop sign is slightly tipped (pitched) by 3 degrees = 0.05 radians
	 * m1 = tan(0.07) = 0.07
	 * b1 = 5.0 - m1 * y0
	 * m2 = tan(0.05) = 0.05
	 * b2 = 2.0 - m2 * y0
	 */
	double x0 = 5.0, y0 = -3.0, z0 = 2.0;
	double m1 = 0.07;
	double b1 = x0 - m1 * y0;
	double m2 = 0.05;
	double b2 = z0 - m2 * y0;
	double inputLineParameters[4] = { m1, b1, m2, b2 };
	int nPoints = 100;
	Vector3 *points = new Vector3[nPoints];
	int iAxis;
	for(iAxis=X_AXIS;iAxis<=Z_AXIS;++iAxis) {
		double endPoints[2] = {-3.5, -2.5};
		double sigma = 0.0002;
		generateLine(points, nPoints, iAxis, inputLineParameters, endPoints, sigma);
		double lineParameters[4] = {0.0, 0.0, 0.0, 0.0};
		getLine(points, nPoints, iAxis, lineParameters, endPoints);
		printf("input line parameters = (%f, %f, %f, %f\n",
			   inputLineParameters[0], inputLineParameters[1], inputLineParameters[2], inputLineParameters[3]);
		printf("final line parameters = (%f, %f, %f, %f\n",
			   lineParameters[0], lineParameters[1], lineParameters[2], lineParameters[3]);
	}
	delete [] points;
	return;
}

int mainX(int argc, char **argv) {
	printf("hello world\n");

	double theta0[3];
	theta0[0] = 1.0;
	theta0[1] = 2.0;
	theta0[2] = 3.0;

	int nPoints = 3 * 250;
	Vector3 *points = new Vector3 [ nPoints ];

    rndm = TRandom3();

#if 0
	/*
	 * generate random points on a plane
	 */
	double a = 1.0, b = 1.5, c = 2.0, d = 1.0;
	for(int i=0;i<nPoints;++i) {
		double x = random.Uniform(1.0);
		double y = random.Uniform(1.0);
		double z = (-d - a * x - b * y) / c;
		points[i].x = x;
		points[i].y = y;
		points[i].z = z;
	}
#endif

	/*
	 * generate a grate
	 */
	double alpha = 5.0 * M_PI / 180.0;
	double sigma = 2.0; /* 2 cm uncertainty */
	double sinAlpha = sin(alpha), cosAlpha = cos(alpha);
	double m[3], b[3], length[3];
	Interval grateInterval[3];
	Vector3 origin = { 0.0, 0.0, 0.0 };
	// first work out the coordinates in the stop sign plane with origin at the center
	m[0] = 0.01, m[1] = 0.012, m[2] = 0.08; // nominally 0.01
	b[0] = -20.0, b[1] = -0.1, b[2] = 20.1; // nominally -20, 0, 20
	length[0] = 30.0, length[1] = 42.0, length[2] = 28; // TODO more consistent numbers
	int iGrate, nGrates = 3;
	for(iGrate=0;iGrate<nGrates;++iGrate) {
		Interval *interval = &grateInterval[iGrate];
		double a = 0.5 * length[iGrate];
		interval->x1 = origin.x - a;
		interval->x2 = origin.x + a;
		interval->y1 = origin.y + b[iGrate] - a * sinAlpha; // only y picks up a component
		interval->y2 = origin.y + b[iGrate] + a * sinAlpha;
		interval->z1 = origin.z;
		interval->z2 = origin.z;
	}

    generateGrate(grateInterval, 3, points, nPoints);

	double unitNormal[3];
	double d = getPlane(points, nPoints, unitNormal);
	printf("(a=%.3f) * x + (b=%.3f) * y + (z=%.3f) * z + %f = 0\n", unitNormal[0], unitNormal[1], unitNormal[2], d);

    Vector3 *planePoints = new Vector3 [ nPoints ];

    Vector3 n;
    n.x = unitNormal[0];
    n.y = unitNormal[1];
    n.z = unitNormal[2];
    for(int iPoint=0;iPoint<nPoints;++iPoint) {
        projectPointOntoPlane(n, d, points[iPoint], planePoints[iPoint]);
    }

    /*
     * cluster in y
     */
    std::vector<double> yClusters;
    double y = planePoints[0].y, maxClusterDistance = 8.0;
    yClusters.push_back(y);
    for(int iPoint=1;iPoint<nPoints;++iPoint) {
        double newY = planePoints[iPoint].y;
        std::vector<double>::const_iterator it, last = yClusters.end();
        bool clusterFound = false;
        for(it=yClusters.begin();it!=last;++it) {
            double y = *it;
            if(fabs(y - newY) < maxClusterDistance) { clusterFound = true; break; }
        }
        if(!clusterFound) { yClusters.push_back(newY); }
    }

    int nClusters = yClusters.size();
    printf("%d clusters found", nClusters);

	/*
	 * for each cluster in y, fit a line
	 */

	Vector3 ***clusterPoints = new Vector3 ** [ nClusters ];
	for(int i=0;i<nClusters;++i) {
		clusterPoints[i] = new Vector3 * [ nPoints ]; // prepare for maximum size
	}
	int *clusterSize = new int [ nClusters ];
	double *sumX = new double [ nClusters ];
	double *sumY = new double [ nClusters ];
	double *sumZ = new double [ nClusters ];
	double *sumXX = new double [ nClusters ];
	double *sumXY = new double [ nClusters ];
	double *sumXZ = new double [ nClusters ];
	double *sumYY = new double [ nClusters ];
	double *sumYZ = new double [ nClusters ];
	double *sumZZ = new double [ nClusters ];
	for(int iCluster=0;iCluster<nClusters;++iCluster) { clusterSize[iCluster] = 0; }
	for(int iCluster=0;iCluster<nClusters;++iCluster) {
		for(int iPoint=0;iPoint<nPoints;++iPoint) {
			Vector3 *point = &planePoints[iPoint];
			double y = point->y;
			double yCluster = yClusters[iCluster];
			int size = 0;
			double sX = 0.0, sY = 0.0, sZ = 0.0, sXX = 0.0, sXY = 0.0, sXZ = 0.0, sYY = 0.0, sYZ = 0.0, sZZ = 0.0;
			if(fabs(y - yCluster) < maxClusterDistance) {
				clusterPoints[iCluster][size] = point;
				double x = point->x;
				double y = point->y;
				double z = point->z;
				sX += x;
				sY += y;
				sXX += x * x;
				sXY += x * y;
				sXZ += x * z;
				sYY += y * y;
				sYZ += y * z;
				sZZ += z * z;
				++size;
			}
			clusterSize[iCluster] = size;
			sumX[iCluster] = sX;
			sumY[iCluster] = sY;
			sumZ[iCluster] = sZ;
			sumXX[iCluster] = sXX;
			sumXY[iCluster] = sXY;
			sumXZ[iCluster] = sXZ;
			sumYY[iCluster] = sYY;
			sumYZ[iCluster] = sYZ;
			sumZZ[iCluster] = sZZ;
		}
	}

	for(int i=0;i<nClusters;++i) { delete [] clusterPoints[i]; }
	delete [] clusterPoints;
	delete [] clusterSize;
	delete [] sumX;
	delete [] sumY;
	delete [] sumZ;
	delete [] sumXX;
	delete [] sumXY;
	delete [] sumXZ;
	delete [] sumYY;
	delete [] sumYZ;
	delete [] sumZZ;
}

// TODO not tested
/*
 * plane is specified by nx * x + ny * y + nz * z + d = 0, where (nx, ny, nz) is a unit vector
 */
void projectPointOntoPlane(Vector3 &n, double d, Vector3 &p0, Vector3 &p) {

    double nx = n.x, ny = n.y, nz = n.z;
    double x0 = p0.x, y0 = p0.y, z0 = p0.z;
    p.x = (ny * ny + nz * nz) * x0 - nz * nx * (z0 + d) - nx * ny * y0;
    p.y = (nx * nx + nz * nz) * y0 - nz * ny * (z0 + d) - nx * ny * x0;
    p.z = -1.0 * (d + p.x * nx / nz + p.y * ny / nz);

    return;

	TMatrixD matrix(3, 3);
	matrix[0][0] = 0.0;
	matrix[0][1] = -n.z;
	matrix[0][2] = n.y;
	matrix[1][0] = n.z;
    matrix[1][1] = 0.0;
	matrix[1][2] = -n.x;
	matrix[2][0] = -n.y;
    matrix[2][1] = n.x;
    matrix[2][2] = 0.0;
	double b[3], x[3];
	b[0] = n.y * p0.z - n.z * p0.y;
	b[1] = n.z * p0.x - n.x * p0.z;
	b[2] = n.x * p0.y - n.y * p0.x;
	solveLinearSystem(matrix, b, x);
	p.x = x[0];
	p.y = x[1];
	p.z = x[2];
}

/*
 * solve Mx = b
 */
double solveLinearSystem(TMatrixD &matrix, Vector3 &b, Vector3 &x) {
	double bNew[3];
	double xNew[3];
	bNew[0] = b.x;
	bNew[1] = b.y;
	bNew[2] = b.z;
	double det = solveLinearSystem(matrix, bNew, xNew);
	x.x = xNew[0];
	x.y = xNew[1];
	x.z = xNew[2];
	return det;
}

double solveLinearSystem(TMatrixD &matrix, double *b, double *x) {
	double det0 = matrix.Determinant();
	if (det0 == 0.0) return det0;
	for(int i=0;i<3;++i) {
		TMatrixD M(3, 3);
		for(int j=0;j<3;++j) { for(int k=0;k<3;++k) { M[j][k] = matrix[j][k]; } }
		for(int j=0;j<3;++j) { M[j][i] = b[j]; }
		double det = M.Determinant();
		x[i] = det / det0;
	}
	return det0;
}

void generateLine(Vector3 *points, int nPoints, int primaryAxis, double lineParameters[4], double endPoints[2], double sigma) {
	int iPoint;
	double tBegin = endPoints[0], tFinal = endPoints[1];
	double deltaT = (tFinal - tBegin) / (nPoints - 1);
	double x, y, z;
	double m1 = lineParameters[0], b1 = lineParameters[1], m2 = lineParameters[2], b2 = lineParameters[3];
	Vector3 *point = points;
	for(iPoint=0;iPoint<nPoints;++iPoint,++point) {
		switch(primaryAxis) {
			case 1:
				x = tBegin + iPoint * deltaT;
				y = m1 * x + b1;
				z = m2 * x + b2;
				break;
			case 2:
				y = tBegin + iPoint * deltaT;
				x = m1 * y + b1;
				z = m2 * y + b2;
				break;
			case 3:
				z = tBegin + iPoint * deltaT;
				x = m1 * z + b1;
				y = m2 * z + b2;
				break;
		}
		point->x = rndm.Gaus(x, sigma);
		point->y = rndm.Gaus(y, sigma);
		point->z = rndm.Gaus(z, sigma);
		// printf("point(%d) = (%f, %f, %f)\n", iPoint, point->x, point->y, point->z);
	}
}

double getLine(Vector3 *points, int nPoints, int primaryAxis, double *lineParameters, double endpoints[2]) {
	double sumXX = 0.0, sumXY = 0.0, sumXZ = 0.0, sumYY = 0.0, sumYZ = 0.0, sumZZ = 0.0;
	double sumX = 0.0, sumY = 0.0, sumZ = 0.0, minOrdinate, maxOrdinate;
	switch(primaryAxis) {
		case X_AXIS:
			minOrdinate = points[0].x;
			maxOrdinate = points[0].x;
			break;
		case Y_AXIS:
			minOrdinate = points[0].y;
			maxOrdinate = points[0].y;
			break;
		case Z_AXIS:
			minOrdinate = points[0].z;
			maxOrdinate = points[0].z;
			break;
	}
	for (int i = 0; i < nPoints; ++i) {
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
		switch (primaryAxis) {
			case 1:
				if (x < minOrdinate) minOrdinate = x;
				if (x > maxOrdinate) maxOrdinate = x;
				break;
			case 2:
				if (y < minOrdinate) minOrdinate = y;
				if (y > maxOrdinate) maxOrdinate = y;
				break;
			case 3:
				if (z < minOrdinate) minOrdinate = z;
				if (z > maxOrdinate) maxOrdinate = z;
				break;
		}
	}
	double sum1 = nPoints;
	double a, aInv;
	switch (primaryAxis) {
		case 1:
			a = sumX * sumX - sumXX * sum1;
			aInv = 1.0 / a;
			lineParameters[0] = (sumX * sumY - sumXY * sum1) * aInv;
			lineParameters[1] = (sumX * sumXY - sumY * sumXX) * aInv;
			lineParameters[2] = (sumX * sumZ - sumXZ * sum1) * aInv;
			lineParameters[3] = (sumX * sumXZ - sumZ * sumXX) * aInv;
			break;
		case 2:
			a = sumY * sumY - sumYY * sum1;
			aInv = 1.0 / a;
			lineParameters[0] = (sumX * sumY - sumXY * sum1) * aInv;
			lineParameters[1] = (sumY * sumXY - sumX * sumYY) * aInv;
			lineParameters[2] = (sumY * sumZ - sumYZ * sum1) * aInv;
			lineParameters[3] = (sumY * sumYZ - sumZ * sumYY) * aInv;
			break;
		case 3:
			a = sumZ * sumZ - sumZZ * sum1;
			aInv = 1.0 / a;
			lineParameters[0] = (sumX * sumZ - sumXZ * sum1) * aInv;
			lineParameters[1] = (sumZ * sumXZ - sumX * sumZZ) * aInv;
			lineParameters[2] = (sumY * sumZ - sumYZ * sum1) * aInv;
			lineParameters[3] = (sumZ * sumYZ - sumY * sumZZ) * aInv;
			break;
	}
	endpoints[0] = minOrdinate;
	endpoints[1] = maxOrdinate;
	return a / nPoints;
}

double getPlane(Vector3 *points, int nPoints, double *unitNormal) {
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

    double normal[3];
	double det0 = solveLinearSystem(matrix, n, normal);
    double d = sqrt(normal[0] * normal[0] + normal[1] * normal[1] + normal[2] * normal[2]);

    unitNormal[0] = normal[0] / d;
    unitNormal[1] = normal[1] / d;
    unitNormal[2] = normal[2] / d;

	return d;
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

int generateGrate(Interval *grates, int nGrates, Vector3 *points, int nPoints) {
	double acc = 0.0, length = 0.0, sigma = 0.02;
	int iGrate, index = 0;
	for(iGrate=0;iGrate<nGrates;++iGrate) { length += intervalLength(&grates[iGrate]); }
	int *n = new int [ nGrates ];
	acc = 0.0;
	for(iGrate=0;iGrate<nGrates;++iGrate) {
		acc += intervalLength(&grates[iGrate]);
		double fraction = acc / length;
		n[iGrate] = (int)(fraction * nPoints);
		if(n[iGrate] > nPoints) {
			n[iGrate] = nPoints;
		}
	}
	for(iGrate=(nGrates-1);iGrate>0;--iGrate) {
		n[iGrate] = n[iGrate] - n[iGrate-1];
	}
	for(iGrate=0;iGrate<nGrates;++iGrate) {
        Interval *grate = &grates[iGrate];
		for (int j = 0; j < n[iGrate]; ++j, ++index) {
			double t = (double)(n[iGrate] - 1); /* -1 so range of j covers full range of interval */
			t = (double) j / t;
			double x = grate->x1 + t * (grate->x2 - grate->x1);
			double y = grate->y1 + t * (grate->y2 - grate->y1);
			double z = grate->z1 + t * (grate->z2 - grate->z1);
			points[index].x = rndm.Gaus(x, sigma);
			points[index].y = rndm.Gaus(y, sigma);
			points[index].z = rndm.Gaus(z, sigma);
			if(index >= nPoints) { break; } /* shouldn't really happen */
		}
	}
	return 0;
}

#if 0

int generateGrate(double *m, double *b, double *delta, double *l, int nGrates,
	Vector3 points, int nPoints) {
    int i, j, k, kPrev, m, n;
	double acc = 0.0, length = 0.0;
	for(i=0;i<nGrates;++i) { length += l[i]; }
	for(i=j=0;i<nGrates;++i) {
		acc += l[i];
		double fraction = acc / length;
		kPrev = k;
		k = fraction * nPoints;
		n = k - kPrev;
		if(k > nPoints) k = nPoints;
		double x1 = -l[i] / 2.0, x2 = l[i] / 2.0;
		double y1 = b[i] - 0.1, y2 = b[i] + 0.1;
		double z1 = -l[i] / 2.0 * tanAlpha, z2 = l[i] / 2.0 * tanAlpha;
		for (m = 0; j < k; ++j, ++m) {
			double t = (double) m / (double) n;
			double x = x1 + t * (x2 - x1);
			double y = y1 + t * (y2 - y1);
			double z = z1 + t * (z2 - z1);
			points[j].x = random.Gaus(x, sigma);
			points[j].y = random.Gaus(y, sigma);
			points[j].z = random.Gaus(z, sigma);
		}
	}
    for(j=0;j<3;++j) {
        int index = 0;
    }
}

#endif