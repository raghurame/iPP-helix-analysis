#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <assert.h>
#include "postprocess.h"

CURVE findCurveParameters (float *OOP)
{
	CURVE cosineParameters;
	float OOP_avg, OOP_sum = 0;

	for (int i = 0; i < 360; ++i)
	{
		OOP_sum += OOP[i];
	}
	cosineParameters.center = OOP_sum / 360;

	for (int i = 0; i < 360; ++i)
	{
		if (i == 0)
		{
			cosineParameters.localMinima_Xpos = i;
			cosineParameters.localMaxima_Xpos = i;
			cosineParameters.localMinima_Yval = OOP[i];
			cosineParameters.localMaxima_Yval = OOP[i];
		}

		if (cosineParameters.localMaxima_Yval < OOP[i])
		{
			cosineParameters.localMaxima_Yval = OOP[i];
			cosineParameters.localMaxima_Xpos = i;
		}
		if (cosineParameters.localMinima_Yval > OOP[i])
		{
			cosineParameters.localMinima_Yval = OOP[i];
			cosineParameters.localMinima_Xpos = i;
		}
	}

	return cosineParameters;
}

float computeCosineCurve (float c1, float c2, float c3, float c4, float theta)
{
	float cosineValue = c1 * cos (degreeToRadians((theta * c2) + c3)) + c4;
	return cosineValue;
}

float computeError (OOP_DATA *oopinfo, OOP_DATA *modelCurve, float precision)
{
	int nDataPoints = (360.0 / precision);
	float error = 0;

	for (int i = 0; i < nDataPoints; ++i)
	{
		error += fabs (oopinfo[i].OOP - modelCurve[i].OOP) / oopinfo[i].OOP;
	}
	error /= nDataPoints;
	error *= 100;
	return error;
}

COSINE findAverage_Cos (COSINE *cosParameters, float precision)
{
	int nDataPoints = (360.0 / precision);
	COSINE averageCosine;

	averageCosine.C1 = 0;
	averageCosine.C2 = 0;
	averageCosine.C3 = 0;
	averageCosine.C4 = 0;
	averageCosine.errorPercent = 0;
	averageCosine.covariance = 0;
	averageCosine.phaseShift = 0;

	for (int i = 0; i < nDataPoints; ++i)
	{
		averageCosine.C1 += cosParameters[i].C1;
		averageCosine.C2 += cosParameters[i].C2;
		averageCosine.C3 += cosParameters[i].C3;
		averageCosine.C4 += cosParameters[i].C4;

		averageCosine.errorPercent += cosParameters[i].errorPercent;
		averageCosine.covariance += cosParameters[i].covariance;
		averageCosine.phaseShift += cosParameters[i].phaseShift;
	}

	averageCosine.C1 /= nDataPoints;
	averageCosine.C2 /= nDataPoints;
	averageCosine.C3 /= nDataPoints;
	averageCosine.C4 /= nDataPoints;

	averageCosine.errorPercent /= nDataPoints;
	averageCosine.covariance /= nDataPoints;
	averageCosine.phaseShift /= nDataPoints;

	return averageCosine;
}

COSINE findSD_Cos (COSINE *cosParameters, float precision, COSINE averageCosine)
{
	int nDataPoints = (360.0 / precision);
	COSINE sdCosine;

	sdCosine.C1 = 0;
	sdCosine.C2 = 0;
	sdCosine.C3 = 0;
	sdCosine.C4 = 0;
	sdCosine.errorPercent = 0;
	sdCosine.covariance = 0;
	sdCosine.phaseShift = 0;

	for (int i = 0; i < nDataPoints; ++i)
	{
		sdCosine.C1 += pow ((cosParameters[i].C1 - averageCosine.C1), 2);
		sdCosine.C2 += pow ((cosParameters[i].C2 - averageCosine.C2), 2);
		sdCosine.C3 += pow ((cosParameters[i].C3 - averageCosine.C3), 2);
		sdCosine.C4 += pow ((cosParameters[i].C4 - averageCosine.C4), 2);

		sdCosine.errorPercent += pow ((cosParameters[i].errorPercent - averageCosine.errorPercent), 2);
		sdCosine.covariance += pow ((cosParameters[i].covariance - averageCosine.covariance), 2);
		sdCosine.phaseShift += pow ((cosParameters[i].phaseShift - averageCosine.phaseShift), 2);
	}

	sdCosine.C1 /= (nDataPoints - 1);
	sdCosine.C2 /= (nDataPoints - 1);
	sdCosine.C3 /= (nDataPoints - 1);
	sdCosine.C4 /= (nDataPoints - 1);

	sdCosine.errorPercent /= (nDataPoints - 1);
	sdCosine.covariance /= (nDataPoints - 1);
	sdCosine.phaseShift /= (nDataPoints - 1);

	sdCosine.C1 = sqrt (sdCosine.C1);
	sdCosine.C2 = sqrt (sdCosine.C2);
	sdCosine.C3 = sqrt (sdCosine.C3);
	sdCosine.C4 = sqrt (sdCosine.C4);

	sdCosine.errorPercent = sqrt (sdCosine.errorPercent);
	sdCosine.covariance = sqrt (sdCosine.covariance);
	sdCosine.phaseShift = sqrt (sdCosine.phaseShift);

	return sdCosine;
}

int findStats (COSINE *averageCosine, COSINE *sdCosine, COSINE *cosParameters, float precision)
{
	int nDataPoints = (360.0 / precision);

	for (int i = 0; i < nDataPoints; ++i)
	{
		(*averageCosine).C1 += cosParameters[i].C1;
		(*averageCosine).C2 += cosParameters[i].C2;
		(*averageCosine).C3 += cosParameters[i].C3;
		(*averageCosine).C4 += cosParameters[i].C4;

		(*averageCosine).errorPercent += cosParameters[i].errorPercent;
		(*averageCosine).covariance += cosParameters[i].covariance;
		(*averageCosine).phaseShift += cosParameters[i].phaseShift;
	}

	(*averageCosine).C1 /= nDataPoints;
	(*averageCosine).C2 /= nDataPoints;
	(*averageCosine).C3 /= nDataPoints;
	(*averageCosine).C4 /= nDataPoints;

	(*averageCosine).errorPercent /= nDataPoints;
	(*averageCosine).covariance /= nDataPoints;
	(*averageCosine).phaseShift /= nDataPoints;

	for (int i = 0; i < nDataPoints; ++i)
	{
		(*sdCosine).C1 += pow ((cosParameters[i].C1 - (*averageCosine).C1), 2);
		(*sdCosine).C2 += pow ((cosParameters[i].C2 - (*averageCosine).C2), 2);
		(*sdCosine).C3 += pow ((cosParameters[i].C3 - (*averageCosine).C3), 2);
		(*sdCosine).C4 += pow ((cosParameters[i].C4 - (*averageCosine).C4), 2);

		(*sdCosine).errorPercent += pow ((cosParameters[i].errorPercent - (*averageCosine).errorPercent), 2);
		(*sdCosine).covariance += pow ((cosParameters[i].covariance - (*averageCosine).covariance), 2);
		(*sdCosine).phaseShift += pow ((cosParameters[i].phaseShift - (*averageCosine).phaseShift), 2);
	}

	(*sdCosine).C1 /= (nDataPoints - 1);
	(*sdCosine).C2 /= (nDataPoints - 1);
	(*sdCosine).C3 /= (nDataPoints - 1);
	(*sdCosine).C4 /= (nDataPoints - 1);

	(*sdCosine).errorPercent /= (nDataPoints - 1);
	(*sdCosine).covariance /= (nDataPoints - 1);
	(*sdCosine).phaseShift /= (nDataPoints - 1);

	(*sdCosine).C1 = sqrt ((*sdCosine).C1);
	(*sdCosine).C2 = sqrt ((*sdCosine).C2);
	(*sdCosine).C3 = sqrt ((*sdCosine).C3);
	(*sdCosine).C4 = sqrt ((*sdCosine).C4);

	(*sdCosine).errorPercent = sqrt ((*sdCosine).errorPercent);
	(*sdCosine).covariance = sqrt ((*sdCosine).covariance);
	(*sdCosine).phaseShift = sqrt ((*sdCosine).phaseShift);

	return 1;

}

float computeCovariance (OOP_DATA *oopinfo, OOP_DATA *modelCurve, float precision)
{
	int nDataPoints = (360.0 / precision);

	float covariance = 0, originalAvg, modelAvg, originalSum = 0, modelSum = 0;

	for (int i = 0; i < nDataPoints; ++i)
	{
		originalAvg += oopinfo[i].OOP;
		modelAvg += modelCurve[i].OOP;
	}
	originalAvg /= nDataPoints;
	modelAvg /= nDataPoints;

	for (int i = 0; i < nDataPoints; ++i)
	{
		covariance += ((oopinfo[i].OOP - originalAvg) * (modelCurve[i].OOP - modelAvg));
	}
	covariance /= (nDataPoints - 1);

	return covariance;
}

void updateCurveParameter (COSINE **cosParameters, int nTimeframes, CURVE_HIGHPRECISION curveParameters, float precision, OOP_DATA *oopinfo)
{
	// Assigning constants
	// printf("maxima: %f; minima: %f\n", curveParameters.localMaxima_Xpos, curveParameters.localMinima_Xpos);

	(*cosParameters)[nTimeframes].C1 = curveParameters.localMaxima_Yval - curveParameters.center;

	if ((fabs (curveParameters.localMaxima_Xpos - curveParameters.localMinima_Xpos) > 100))
		(*cosParameters)[nTimeframes].C2 = 540 / fabs (curveParameters.localMaxima_Xpos - curveParameters.localMinima_Xpos);
	else if ((fabs (curveParameters.localMaxima_Xpos - curveParameters.localMinima_Xpos)) <= 100)
		(*cosParameters)[nTimeframes].C2 = 180 / fabs (curveParameters.localMaxima_Xpos - curveParameters.localMinima_Xpos);

	(*cosParameters)[nTimeframes].C3 = -curveParameters.localMaxima_Xpos * (*cosParameters)[nTimeframes].C2;

	(*cosParameters)[nTimeframes].phaseShift = curveParameters.localMaxima_Xpos;
	
	(*cosParameters)[nTimeframes].C4 = curveParameters.center;

	if ((*cosParameters)[nTimeframes].phaseShift < -360)
		(*cosParameters)[nTimeframes].phaseShift += 360;

	if ((*cosParameters)[nTimeframes].phaseShift > 180)
		(*cosParameters)[nTimeframes].phaseShift -= 180;

	if ((*cosParameters)[nTimeframes].phaseShift > 90)
		(*cosParameters)[nTimeframes].phaseShift -= 90;

	// Generating model cosine curve
	int nDataPoints = (360.0 / precision);
	float theta = 0.0;
	OOP_DATA *modelCurve;
	modelCurve = (OOP_DATA *) malloc (nDataPoints * sizeof (OOP_DATA));

	for (int i = 0; i < nDataPoints; ++i)
	{
		modelCurve[i].OOP = computeCosineCurve ((*cosParameters)[nTimeframes].C1, (*cosParameters)[nTimeframes].C2, (*cosParameters)[nTimeframes].C3, (*cosParameters)[nTimeframes].C4, theta);
		theta += precision;
	}

	// Finding errors
	(*cosParameters)[nTimeframes].errorPercent = computeError (oopinfo, modelCurve, precision);

	// Finding covariance
	(*cosParameters)[nTimeframes].covariance = computeCovariance (oopinfo, modelCurve, precision);
}

CURVE_HIGHPRECISION findCurveParameters_highPrecision (struct iPP_OOP *oopinfo, float precision)
{
	CURVE_HIGHPRECISION cosineParameters;

	float OOP_avg, OOP_sum = 0, currentTheta = 0.0;
	float nDataPoints = (int) (360.0 / precision);

	for (int i = 0; i < nDataPoints; ++i)
	{
		OOP_sum += oopinfo[i].OOP;
	}
	cosineParameters.center = OOP_sum / nDataPoints;

	for (int i = 0; i < nDataPoints; ++i)
	{
		if (i == 0)
		{
			cosineParameters.localMinima_Xpos = currentTheta;
			cosineParameters.localMaxima_Xpos = currentTheta;
			cosineParameters.localMinima_Yval = oopinfo[i].OOP;
			cosineParameters.localMaxima_Yval = oopinfo[i].OOP;
		}

		if (cosineParameters.localMaxima_Yval < oopinfo[i].OOP)
		{
			cosineParameters.localMaxima_Yval = oopinfo[i].OOP;
			cosineParameters.localMaxima_Xpos = currentTheta;
		}
		if (cosineParameters.localMinima_Yval > oopinfo[i].OOP)
		{
			cosineParameters.localMinima_Yval = oopinfo[i].OOP;
			cosineParameters.localMinima_Xpos = currentTheta;
		}
		currentTheta += precision;
	}

	return cosineParameters;
}

int makeDirectory (const char *folderName)
{
	struct stat checkDirectory = {0};
	if (stat (folderName, &checkDirectory) == -1)
	{
		mkdir (folderName, 0700);
		return 1;
	}
	else
	{
		return 0;
	}
}

void iPP_printOOP_TAVG (struct iPP_OOP_TAVG *oopinfo_tavg, const char *outputFolder)
{
	char *outputLocation;
	outputLocation = (char *) malloc (100 * sizeof (char));
	sprintf (outputLocation, "%s/%s",outputFolder, "avg_oop.csv");
	FILE *output;
	output = fopen (outputLocation, "w");

	fprintf(output, "%s, %s\n", "theta", "OOP");

	for (int Looptheta = 0; Looptheta < 360; ++Looptheta)
	{
		fprintf(output, "%d, %.4f\n", Looptheta, oopinfo_tavg[Looptheta].OOP_TAVG);
	}
}

void iPP_computeOOP_TAVG (struct iPP_OOP *oopinfo, struct iPP_OOP_TAVG **oopinfo_tavg, int nTimeframes)
{
	for (int Looptheta = 0; Looptheta < 360; ++Looptheta)
	{
		(*oopinfo_tavg)[Looptheta].oop_sum += oopinfo[Looptheta].OOP;
		(*oopinfo_tavg)[Looptheta].OOP_TAVG = (*oopinfo_tavg)[Looptheta].oop_sum / nTimeframes;
	}
}

float rotateX (float angle_radians)
{
	float sineVariable;
	sineVariable = sinf (angle_radians);
	return sineVariable;
}

float rotateZ (float angle_radians)
{
	float cosineVariable;
	cosineVariable = cosf (angle_radians);
	return cosineVariable;
}

struct iPP_OOP *iPP_computeOOP_findXYZ ()
{
	float theta_radians;
	struct iPP_OOP *oopinfo;
	oopinfo = (struct iPP_OOP *) malloc (360 * sizeof (struct iPP_OOP));

	for (int theta = 0; theta < 360; ++theta)
	{
		oopinfo[theta].x1 = 0;
		oopinfo[theta].y1 = 0;
		oopinfo[theta].z1 = 0;
		oopinfo[theta].y2 = 0;

		theta_radians = degreeToRadians(theta);
		oopinfo[theta].x2 = rotateX(theta_radians);
		oopinfo[theta].z2 = rotateZ(theta_radians);
	}

	return oopinfo;
}

struct iPP_OOP *iPP_computeOOP_findXYZ_highPrecision (float precision)
{
	float theta_radians, theta = 0;
	struct iPP_OOP *oopinfo;
	int nDataPoints_current = 0, nDataPoints = (int) (360.0 / precision);
	oopinfo = (struct iPP_OOP *) malloc (nDataPoints * sizeof (struct iPP_OOP));

	for (; nDataPoints_current < nDataPoints; )
	{
		oopinfo[nDataPoints_current].x1 = 0;
		oopinfo[nDataPoints_current].y1 = 0;
		oopinfo[nDataPoints_current].z1 = 0;
		oopinfo[nDataPoints_current].y2 = 0;

		theta_radians = degreeToRadians(theta);
		oopinfo[nDataPoints_current].x2 = rotateX(theta_radians);
		oopinfo[nDataPoints_current].z2 = rotateZ(theta_radians);

		nDataPoints_current += 1;
		theta += precision;
	}

	return oopinfo;
}


struct iPP_OOP *iPP_computeOOP (struct iPP_IOP *iopinfo, int nIOPs)
{
	float vx1, vx2, vy1, vy2, vz1, vz2, magnitude1, magnitude2, dotproduct, costheta, theta, oop_local, oop_sum = 0;
	// int nDataPoints = (int) (360.0 / 0.0001);
	// printf("nDataPoints: %d\n", nDataPoints);

	struct iPP_OOP *oopinfo;
	oopinfo = (struct iPP_OOP *) malloc (360 * sizeof (struct iPP_OOP));

	// int progressWidth = 30;

	/* Calculating x1, x2, y1, y2, z1, z2 */
	oopinfo = iPP_computeOOP_findXYZ ();

	for (int Looptheta = 0; Looptheta < 360; ++Looptheta)
	{
		vx2 = oopinfo[Looptheta].x2 - oopinfo[Looptheta].x1;
		vy2 = oopinfo[Looptheta].y2 - oopinfo[Looptheta].y1;
		vz2 = oopinfo[Looptheta].z2 - oopinfo[Looptheta].z1;

		// fprintf(stdout, "(%s: %.4f, %s: %.4f, %s: %.4f), (%s: %.4f, %s: %.4f, %s: %.4f)\n", "x1", oopinfo[Looptheta].x1, "y1", oopinfo[Looptheta].y1, "z1", oopinfo[Looptheta].z1, "x2", oopinfo[Looptheta].x2, "y2", oopinfo[Looptheta].y2, "z2", oopinfo[Looptheta].z2);

		oop_sum = 0;

		for (int nIOP = 0; nIOP < nIOPs; nIOP++)
		{
			vx1 = iopinfo[nIOP].x2 - iopinfo[nIOP].x1;
			vy1 = iopinfo[nIOP].y2 - iopinfo[nIOP].y1;
			vz1 = iopinfo[nIOP].z2 - iopinfo[nIOP].z1;

			dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
			magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
			magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
			costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
			theta = acosf (costheta);

			oop_local = ((3 * costheta * costheta) - 1) / 2;
			oop_sum += oop_local;
		}

		oopinfo[Looptheta].OOP = oop_sum / nIOPs;
		// fprintf(stdout, "%s: %.5f (%s: %d)\n", "oop_sum", oop_sum, "Looptheta", Looptheta);
		// fprintf(stdout, "%s: %.4f\n", "oopinfo[Looptheta].OOP", oopinfo[Looptheta].OOP);
		// if (Looptheta % 50 == 0)
		// 	sleep(1);

	}
	return oopinfo;
}

struct iPP_OOP *iPP_computeOOP_highPrecision (struct iPP_IOP *iopinfo, int nIOPs, float precision, struct iPP_OOP *oopinfo)
{
	float vx1, vx2, vy1, vy2, vz1, vz2, magnitude1, magnitude2, dotproduct, costheta, theta, oop_local, oop_sum = 0, Looptheta = 0;
	int nDataPoints = (int) (360.0 / precision), nDataPoints_current = 0;

	for (; nDataPoints_current < nDataPoints; )
	{
		// fprintf(stdout, "\r%f; %d", Looptheta, nDataPoints_current);

		vx2 = oopinfo[nDataPoints_current].x2 - oopinfo[nDataPoints_current].x1;
		vy2 = oopinfo[nDataPoints_current].y2 - oopinfo[nDataPoints_current].y1;
		vz2 = oopinfo[nDataPoints_current].z2 - oopinfo[nDataPoints_current].z1;

		oop_sum = 0;

		for (int nIOP = 0; nIOP < nIOPs; nIOP++)
		{
			vx1 = iopinfo[nIOP].x2 - iopinfo[nIOP].x1;
			vy1 = iopinfo[nIOP].y2 - iopinfo[nIOP].y1;
			vz1 = iopinfo[nIOP].z2 - iopinfo[nIOP].z1;

			dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
			magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
			magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
			costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
			theta = acosf (costheta);

			oop_local = ((3 * costheta * costheta) - 1) / 2;
			oop_sum += oop_local;
		}

		oopinfo[nDataPoints_current].OOP = oop_sum / nIOPs;
		Looptheta += precision;
		if (nDataPoints_current < nDataPoints)
			nDataPoints_current += 1;
	}
	return oopinfo;
}

struct iPP_OOP *iPP_computeOOP_highPrecision_2D (struct iPP_IOP *iopinfo, int nIOPs, float precision, struct iPP_OOP *oopinfo)
{
	float vx1, vx2, vy1, vy2, vz1, vz2, magnitude1, magnitude2, dotproduct, costheta, theta, oop_local, oop_sum = 0, Looptheta = 0;
	int nDataPoints = (int) (360.0 / precision), nDataPoints_current = 0;

	for (; nDataPoints_current < nDataPoints; )
	{
		// fprintf(stdout, "\r%f; %d", Looptheta, nDataPoints_current);

		vx2 = oopinfo[nDataPoints_current].x2 - oopinfo[nDataPoints_current].x1;
		vy2 = oopinfo[nDataPoints_current].y2 - oopinfo[nDataPoints_current].y1;
		vz2 = oopinfo[nDataPoints_current].z2 - oopinfo[nDataPoints_current].z1;

		oop_sum = 0;

		for (int nIOP = 0; nIOP < nIOPs; nIOP++)
		{
			vx1 = iopinfo[nIOP].x2 - iopinfo[nIOP].x1;
			vy1 = iopinfo[nIOP].y2 - iopinfo[nIOP].y1;
			vz1 = iopinfo[nIOP].z2 - iopinfo[nIOP].z1;

			dotproduct = (vx1 * vx2) + (vz1 * vz2);
			magnitude1 = (vx1 * vx1) + (vz1 * vz1);
			magnitude2 = (vx2 * vx2) + (vz2 * vz2);
			costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
			theta = acosf (costheta);

			oop_local = ((3 * costheta * costheta) - 1) / 2;
			oop_sum += oop_local;
		}

		oopinfo[nDataPoints_current].OOP = oop_sum / nIOPs;
		Looptheta += precision;
		if (nDataPoints_current < nDataPoints)
			nDataPoints_current += 1;
	}
	return oopinfo;
}

struct iPP_OOP *iPP_computeOOP_highPrecision_2D_bounds (struct iPP_IOP *iopinfo, int nIOPs, float precision, struct iPP_OOP *oopinfo)
{
	float vx1, vx2, vy1, vy2, vz1, vz2, magnitude1, magnitude2, dotproduct, costheta, theta, oop_local, oop_sum = 0, Looptheta = 0;
	int nDataPoints = (int) (360.0 / precision), nDataPoints_current = 0;
	float iop_com_x, iop_com_y, iop_com_z;
	float yhi_bounds = -13, ylo_bounds = -30;

	for (; nDataPoints_current < nDataPoints; )
	{
		// fprintf(stdout, "\r%f; %d", Looptheta, nDataPoints_current);

		vx2 = oopinfo[nDataPoints_current].x2 - oopinfo[nDataPoints_current].x1;
		vy2 = oopinfo[nDataPoints_current].y2 - oopinfo[nDataPoints_current].y1;
		vz2 = oopinfo[nDataPoints_current].z2 - oopinfo[nDataPoints_current].z1;

		oop_sum = 0;

		for (int nIOP = 0; nIOP < nIOPs; nIOP++)
		{
			iop_com_x = (iopinfo[nIOP].x2 + iopinfo[nIOP].x1) / 2;
			iop_com_y = (iopinfo[nIOP].y2 + iopinfo[nIOP].y1) / 2;
			iop_com_z = (iopinfo[nIOP].z2 + iopinfo[nIOP].z1) / 2;

			if (iop_com_y <= yhi_bounds && iop_com_y > ylo_bounds)
			{
				vx1 = iopinfo[nIOP].x2 - iopinfo[nIOP].x1;
				vy1 = iopinfo[nIOP].y2 - iopinfo[nIOP].y1;
				vz1 = iopinfo[nIOP].z2 - iopinfo[nIOP].z1;

				dotproduct = (vx1 * vx2) + (vz1 * vz2);
				magnitude1 = (vx1 * vx1) + (vz1 * vz1);
				magnitude2 = (vx2 * vx2) + (vz2 * vz2);
				costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
				theta = acosf (costheta);

				oop_local = ((3 * costheta * costheta) - 1) / 2;				
			}
			else
			{
				oop_local = 0;
			}
			oop_sum += oop_local;
		}

		oopinfo[nDataPoints_current].OOP = oop_sum / nIOPs;
		Looptheta += precision;
		if (nDataPoints_current < nDataPoints)
			nDataPoints_current += 1;
	}
	return oopinfo;
}

void *iPP_computeRIOP_TAVG (struct iPP_RIOP *riopinfo, struct iPP_RIOP_TAVG **riopinfo_tavg, int nTimeframes, float simBoxSize, float dr)
{
	int nRIOPs = (int) simBoxSize / dr;

	for (int nRIOP = 0; nRIOP < nRIOPs; ++nRIOP)
	{
		(*riopinfo_tavg)[nRIOP].r_min = riopinfo[nRIOP].r_min;
		(*riopinfo_tavg)[nRIOP].r_max = riopinfo[nRIOP].r_max;
		if (isnan(riopinfo[nRIOP].RIOP))
			(*riopinfo_tavg)[nRIOP].riop_sum += 0;
		else
			(*riopinfo_tavg)[nRIOP].riop_sum += riopinfo[nRIOP].RIOP;
		(*riopinfo_tavg)[nRIOP].RIOP_TAVG = (*riopinfo_tavg)[nRIOP].riop_sum / nTimeframes;
	}
	return 0;
}

void iPP_printRIOP_TAVG (struct iPP_RIOP_TAVG *riopinfo_tavg, int nRIOPs, const char *outputFolder)
{
	char *outputLocation;
	outputLocation = (char *) malloc (100 * sizeof (char));
	sprintf (outputLocation, "%s/%s",outputFolder, "avg_riop.csv");
	FILE *output;
	output = fopen (outputLocation, "w");

	fprintf(output, "%s, %s, %s\n", "r_min", "r_max", "RIOP_TAVG");

	for (int nRIOP = 0; nRIOP < nRIOPs; ++nRIOP)
	{
		// fprintf(stdout, "%s: %.2f; %s: %.2f; %s: %.4f\n", "r_min", riopinfo_tavg[nRIOP].r_min, "r_max", riopinfo_tavg[nRIOP].r_max, "RIOP", riopinfo_tavg[nRIOP].RIOP_TAVG);
		fprintf(output, "%.2f, %.2f, %.4f\n", riopinfo_tavg[nRIOP].r_min, riopinfo_tavg[nRIOP].r_max, riopinfo_tavg[nRIOP].RIOP_TAVG);
	}
}

struct iPP_IOP *iPP_computeIOP (struct lammpsdump *dumpinfo, struct isotactic_polypropylene *monomers, int nMonomers, int *nIOPs)
{
	int nIOPs_local = 0;
	assert (nIOPs);

	struct iPP_IOP *iopinfo;
	iopinfo = (struct iPP_IOP *) malloc (nMonomers * sizeof (struct iPP_IOP));

	for (int i = 0; i < nMonomers-3; ++i)
	{
		if (monomers[i].chain_id_number == monomers[i+3].chain_id_number)
		{
			iopinfo[nIOPs_local].a = monomers[i].CH2;
			iopinfo[nIOPs_local].b = monomers[i].CH;
			iopinfo[nIOPs_local].x1 = (dumpinfo[monomers[i].CH2 - 1].x_dump_unwrapped + dumpinfo[monomers[i].CH - 1].x_dump_unwrapped) / 2;
			iopinfo[nIOPs_local].y1 = (dumpinfo[monomers[i].CH2 - 1].y_dump_unwrapped + dumpinfo[monomers[i].CH - 1].y_dump_unwrapped) / 2;
			iopinfo[nIOPs_local].z1 = (dumpinfo[monomers[i].CH2 - 1].z_dump_unwrapped + dumpinfo[monomers[i].CH - 1].z_dump_unwrapped) / 2;
			
			iopinfo[nIOPs_local].c = monomers[i+3].CH2;
			iopinfo[nIOPs_local].d = monomers[i+3].CH;			
			iopinfo[nIOPs_local].x2 = (dumpinfo[monomers[i+3].CH2 - 1].x_dump_unwrapped + dumpinfo[monomers[i+3].CH - 1].x_dump_unwrapped) / 2;
			iopinfo[nIOPs_local].y2 = (dumpinfo[monomers[i+3].CH2 - 1].y_dump_unwrapped + dumpinfo[monomers[i+3].CH - 1].y_dump_unwrapped) / 2;
			iopinfo[nIOPs_local].z2 = (dumpinfo[monomers[i+3].CH2 - 1].z_dump_unwrapped + dumpinfo[monomers[i+3].CH - 1].z_dump_unwrapped) / 2;

			/* Finding center of IOPs */
			iopinfo[nIOPs_local].x_center = (iopinfo[i].x1 + iopinfo[i].x2) / 2;
			iopinfo[nIOPs_local].y_center = (iopinfo[i].y1 + iopinfo[i].y2) / 2;
			iopinfo[nIOPs_local].z_center = (iopinfo[i].z1 + iopinfo[i].z2) / 2;

			/* Checking */
			// fprintf(stdout, "%d ===> (%s: %.2f, %s: %.2f, %s: %.2f), (%s: %.2f, %s: %.2f, %s: %.2f)\n", nIOPs_local, "x1", iopinfo[nIOPs_local].x1, "y1", iopinfo[nIOPs_local].y1, "z1", iopinfo[nIOPs_local].z1, "x2", iopinfo[nIOPs_local].x2, "y2", iopinfo[nIOPs_local].y2, "z2", iopinfo[nIOPs_local].z2);
			// fprintf(stdout, "%s: %d; %s: %d; %s: %d; %s: %d;\n", "a", iopinfo[nIOPs_local].a, "b", iopinfo[nIOPs_local].b, "c", iopinfo[nIOPs_local].c, "d", iopinfo[nIOPs_local].d);
			// fflush(stdout);
			// sleep(1);

			nIOPs_local++;
		}
	}

	*nIOPs = nIOPs_local;
	return iopinfo;
}

void displayProgressbar (int dx, int x, int progressWidth, const char *processName)
{
	float progressBar_float = (float) (dx + 1) * progressWidth / x;
	int isFirst = 1;

	fprintf(stdout, "\rProcessing %s --> %s", processName, "[");
	for (int i = 0; i < progressWidth; ++i)
	{
		if (i < progressBar_float)
			fprintf(stdout, "%s", "=");
		else
		{
			if (isFirst)
			{
				fprintf(stdout, "%s", ">");
				// fflush (stdout);
				isFirst = 0;
			}
			else
			{
				fprintf(stdout, "%s", " ");
				// fflush (stdout);
			}
		}
	}
	fprintf(stdout, "%s [%.0f%%]", "]", (float) progressBar_float * 100 / progressWidth);
	// fflush(stdout);

}

struct iPP_RIOP *iPP_computeRIOP (struct iPP_IOP *iopinfo, struct lammpsdump *dumpinfo, struct isotactic_polypropylene *monomers, int nIOPs, float dr, float simBoxSize)
{
	struct iPP_RIOP *riopinfo;
	int nRIOPs = simBoxSize / dr;
	riopinfo = (struct iPP_RIOP *) malloc (nRIOPs * sizeof (struct iPP_RIOP));

	float distance_iop, riop_local;
	int denom = 0;
	int progressBar = 0, progressWidth = 30;
	float progressBar_float = 0;
	float vx1, vx2, vy1, vy2, vz1, vz2, dotproduct, magnitude1, magnitude2, costheta, theta;

	fprintf(stdout, "\n");

	/* Checking iopinfo */
	// for (int nIOP = 0; nIOP < nIOPs; ++nIOP)
	// {
	// 	fprintf(stdout, "(%s: %.2f, %s: %.2f, %s: %.2f), (%s: %.2f, %s: %.2f, %s: %.2f)\n", "x1", iopinfo[nIOP].x1, "y1", iopinfo[nIOP].y1, "z1", iopinfo[nIOP].z1, "x2", iopinfo[nIOP].x2, "y2", iopinfo[nIOP].y2, "z2", iopinfo[nIOP].z2);
	// 	fflush(stdout);
	// 	sleep(1);
	// }
	// exit(1);

	for (int nRIOP = 0; nRIOP < nRIOPs; nRIOP++)
	{

		/* Displaying progress bar */
		displayProgressbar (nRIOP, nRIOPs, progressWidth, "RIOP");

		riopinfo[nRIOP].RIOP = 0;
		denom = 0;

		for (int i = 0; i < nIOPs; ++i)
		{
			for (int j = 0; j < nIOPs; ++j)
			{
				distance_iop = sqrt (
					((iopinfo[i].x_center - iopinfo[j].x_center) * (iopinfo[i].x_center - iopinfo[j].x_center)) + 
					((iopinfo[i].y_center - iopinfo[j].y_center) * (iopinfo[i].y_center - iopinfo[j].y_center)) + 
					((iopinfo[i].z_center - iopinfo[j].z_center) * (iopinfo[i].z_center - iopinfo[j].z_center))
					);

				riopinfo[nRIOP].r_min = nRIOP * dr;
				riopinfo[nRIOP].r_max = riopinfo[nRIOP].r_min + dr;

				if (i != j && distance_iop <= riopinfo[nRIOP].r_max && distance_iop > riopinfo[nRIOP].r_min)
				{
					vx1 = iopinfo[i].x2 - iopinfo[i].x1;
					vy1 = iopinfo[i].y2 - iopinfo[i].y1;
					vz1 = iopinfo[i].z2 - iopinfo[i].z1;
					vx2 = iopinfo[j].x2 - iopinfo[j].x1;
					vy2 = iopinfo[j].y2 - iopinfo[j].y1;
					vz2 = iopinfo[j].z2 - iopinfo[j].z1;

					dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
					magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
					magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
					costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
					theta = acosf (costheta);

					riop_local = ((3 * costheta * costheta) - 1) / 2;
					riopinfo[nRIOP].RIOP += riop_local;
					denom++;
				}
			}
		}

		riopinfo[nRIOP].RIOP /= denom;

		// if (isnan(riopinfo[nRIOP].RIOP))
		// {
		// 	fprintf(stdout, "\n");
		// 	return riopinfo;
		// }
	}
	fprintf(stdout, "\n");

	return riopinfo;
}

struct iPP_LOP *iPP_computeLOP (struct lammpsdump *dumpinfo, struct isotactic_polypropylene *monomers, int nMonomers, int *nLOPs, float *AVERAGE_LOCAL_ORDER_PARAMETER, const char *outputFolder)
{
	assert (nLOPs);
	assert (AVERAGE_LOCAL_ORDER_PARAMETER);

	/* Output file location */
	char *LOPdump_location, *LOPavg_location;
	LOPdump_location = (char *) malloc (100 * sizeof (char));
	sprintf (LOPdump_location, "%s/%s", outputFolder, "dump_lop.csv");
	LOPavg_location = (char *) malloc (100 * sizeof (char));
	sprintf (LOPavg_location, "%s/%s", outputFolder, "avg_lop.csv");

	/* Output file */
	FILE *LOPdump = fopen (LOPdump_location, "w"); 
	FILE *LOPavg = fopen (LOPavg_location, "a");

	int nLOPs_local = 0;
	float sum_lops = 0;

	struct iPP_LOP *lopinfo;
	lopinfo = (struct iPP_LOP *) malloc (nMonomers * sizeof (struct iPP_LOP));

	float vx1, vy1, vz1, vx2, vy2, vz2, dotproduct, magnitude1, magnitude2, costheta, theta, LOCAL_ORDER_PARAMETER;

	for (int i = 0; i < nMonomers; ++i)
	{
		if (i >= 3 && (i+3) < nMonomers && monomers[i].chain_id_number == monomers[i+3].chain_id_number && monomers[i].chain_id_number == monomers[i-3].chain_id_number)
		{
			lopinfo[i].a = monomers[i].CH3;
			lopinfo[i].b = monomers[i+3].CH3;
			lopinfo[i].c = monomers[i].CH2;
			lopinfo[i].d = monomers[i-3].CH2;

			vx1 = dumpinfo[monomers[i].CH3 - 1].x_dump_unwrapped - dumpinfo[monomers[i+3].CH3 - 1].x_dump_unwrapped;
			vy1 = dumpinfo[monomers[i].CH3 - 1].y_dump_unwrapped - dumpinfo[monomers[i+3].CH3 - 1].y_dump_unwrapped;
			vz1 = dumpinfo[monomers[i].CH3 - 1].z_dump_unwrapped - dumpinfo[monomers[i+3].CH3 - 1].z_dump_unwrapped;
			vx2 = dumpinfo[monomers[i].CH2 - 1].x_dump_unwrapped - dumpinfo[monomers[i-3].CH2 - 1].x_dump_unwrapped;
			vy2 = dumpinfo[monomers[i].CH2 - 1].y_dump_unwrapped - dumpinfo[monomers[i-3].CH2 - 1].y_dump_unwrapped;
			vz2 = dumpinfo[monomers[i].CH2 - 1].z_dump_unwrapped - dumpinfo[monomers[i-3].CH2 - 1].z_dump_unwrapped;

			dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
			magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
			magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
			costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
			theta = acosf (costheta);

			lopinfo[i].LOP1 = ((3 * costheta * costheta) - 1) / 2;
			nLOPs_local++;
			sum_lops += lopinfo[i].LOP1;

			fprintf(LOPdump, "%d,%.4f,%d,%d,%d,%d\n", nLOPs_local, lopinfo[i].LOP1, lopinfo[i].a, lopinfo[i].b, lopinfo[i].c, lopinfo[i].d);
		}
		if ((i+5) < nMonomers && monomers[i].chain_id_number == monomers[i+3].chain_id_number && monomers[i].chain_id_number == monomers[i+5].chain_id_number && monomers[i].chain_id_number == monomers[i+2].chain_id_number)
		{
			lopinfo[i].a = monomers[i].CH3;
			lopinfo[i].b = monomers[i+3].CH3;
			lopinfo[i].e = monomers[i+2].CH2;
			lopinfo[i].f = monomers[i+5].CH2;

			vx1 = dumpinfo[lopinfo[i].a - 1].x_dump_unwrapped - dumpinfo[lopinfo[i].b - 1].x_dump_unwrapped;
			vy1 = dumpinfo[lopinfo[i].a - 1].y_dump_unwrapped - dumpinfo[lopinfo[i].b - 1].y_dump_unwrapped;
			vz1 = dumpinfo[lopinfo[i].a - 1].z_dump_unwrapped - dumpinfo[lopinfo[i].b - 1].z_dump_unwrapped;
			vx2 = dumpinfo[lopinfo[i].e - 1].x_dump_unwrapped - dumpinfo[lopinfo[i].f - 1].x_dump_unwrapped;
			vy2 = dumpinfo[lopinfo[i].e - 1].y_dump_unwrapped - dumpinfo[lopinfo[i].f - 1].y_dump_unwrapped;
			vz2 = dumpinfo[lopinfo[i].e - 1].z_dump_unwrapped - dumpinfo[lopinfo[i].f - 1].z_dump_unwrapped;

			dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
			magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
			magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
			costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
			theta = acosf (costheta);

			lopinfo[i].LOP2 = ((3 * costheta * costheta) - 1) / 2;
			nLOPs_local++;
			sum_lops += lopinfo[i].LOP2;

			fprintf(LOPdump, "%d,%.4f,%d,%d,%d,%d\n", nLOPs_local, lopinfo[i].LOP2, lopinfo[i].a, lopinfo[i].b, lopinfo[i].e, lopinfo[i].f);
		}
	}
	*AVERAGE_LOCAL_ORDER_PARAMETER = (sum_lops / nLOPs_local);
	*nLOPs = nLOPs_local;

	fprintf(LOPavg, "%.4f\n", (sum_lops / nLOPs_local));
	fclose (LOPavg);
	fclose (LOPdump);
	return lopinfo;
}

struct isotactic_polypropylene *classifyiPP (int natoms, struct lammpsdump *dumpinfo, int *nMonomers, int *nChains, int nTimeframes)
{
	int currentMonomer = 0, currentChain = 0;

	struct isotactic_polypropylene *monomers;
	monomers = (struct isotactic_polypropylene *) calloc ((natoms / 3), sizeof (struct isotactic_polypropylene));

	for (int i = 0; i < natoms; ++i)
	{
		if ((i+2) < natoms && dumpinfo[i].molType == 1 && dumpinfo[i].atomType == 1 && dumpinfo[i+1].atomType == 3 && dumpinfo[i+2].atomType == 2)
		{
			monomers[currentMonomer].CH = dumpinfo[i].id;
			monomers[currentMonomer].CH3 = dumpinfo[i+1].id;
			monomers[currentMonomer].CH2 = dumpinfo[i+2].id;
			monomers[currentMonomer].chain_id_number = currentChain;

			/* Checking */
			// if (nTimeframes > 508)
			// {
			// 	fprintf(stdout, "%s: %d; %s: %d; %s: %d;\n", "CH", dumpinfo[i].id, "CH3", dumpinfo[i+1].id, "CH2", dumpinfo[i+2].id);
			// 	fprintf(stdout, "%s: %.4f; %s: %.4f; %s: %.4f;\n", "x", dumpinfo[i].x_dump_unwrapped, "y", dumpinfo[i].y_dump_unwrapped, "z", dumpinfo[i].z_dump_unwrapped);
			// 	fprintf(stdout, "%s: %.4f; %s: %.4f; %s: %.4f;\n", "x", dumpinfo[i+1].x_dump_unwrapped, "y", dumpinfo[i+1].y_dump_unwrapped, "z", dumpinfo[i+1].z_dump_unwrapped);
			// 	fprintf(stdout, "%s: %.4f; %s: %.4f; %s: %.4f;\n", "x", dumpinfo[i+2].x_dump_unwrapped, "y", dumpinfo[i+2].y_dump_unwrapped, "z", dumpinfo[i+2].z_dump_unwrapped);
			// 	fflush(stdout);
			// 	if (i % 100 == 0)
			// 		sleep(1);
			// }

			currentMonomer++;
		}
		else if (dumpinfo[i].atomType == 1 && dumpinfo[i+1].atomType == 3 && dumpinfo[i+2].atomType == 3)
		{
			currentChain++;
		}
	}

	*nMonomers = currentMonomer;
	*nChains = currentChain;
	return monomers;
}

struct lammpsdump *lscanf (FILE *read, int natoms)
{
	char lineString[3000];
	struct lammpsdump *dumpinfo;
	dumpinfo = (struct lammpsdump *) malloc (natoms * sizeof (struct lammpsdump));

	int lineNumber = 0;
	float xlo, xhi, ylo, yhi, zlo, zhi;

	if (fgets (lineString, 3000, read) == NULL)
	{
		printf("NULL detected in input dump\n");
		// dumpinfo = (struct lammpsdump *) malloc (1 * sizeof (struct lammpsdump));
		dumpinfo[0].status = 0;
		return dumpinfo;
	}
	else
	{
		while (fgets (lineString, 3000, read) != NULL)
		{
			lineNumber++;

			// if (lineNumber == 1)
			// 	fprintf(stdout, "%s\n", lineString);

			if (lineNumber == 3)
			{
				sscanf (lineString, "%d", &natoms);
				// dumpinfo = (struct lammpsdump *) realloc (dumpinfo, (natoms * 2) * sizeof (struct lammpsdump));
			}

			if (lineNumber == 5)
				sscanf (lineString, "%f %f\n", &xlo, &xhi);
			if (lineNumber == 6)
				sscanf (lineString, "%f %f\n", &ylo, &yhi);
			if (lineNumber == 7)
				sscanf (lineString, "%f %f\n", &zlo, &zhi);

			if (lineNumber > 8 && lineNumber <= (natoms + 8))
			{
				sscanf (lineString, "%d %d %d %f %f %f %f %f %f %d %d %d\n", 
					&dumpinfo[lineNumber-9].id, 
					&dumpinfo[lineNumber-9].molType, 
					&dumpinfo[lineNumber-9].atomType, 
					&dumpinfo[lineNumber-9].x, 
					&dumpinfo[lineNumber-9].y, 
					&dumpinfo[lineNumber-9].z, 
					&dumpinfo[lineNumber-9].xs, 
					&dumpinfo[lineNumber-9].ys, 
					&dumpinfo[lineNumber-9].zs, 
					&dumpinfo[lineNumber-9].ix, 
					&dumpinfo[lineNumber-9].iy, 
					&dumpinfo[lineNumber-9].iz);

				/* Unwrapping coordinates */
				dumpinfo[lineNumber-9].x_dump_unwrapped = unwrapCoordinates (dumpinfo[lineNumber-9].x, dumpinfo[lineNumber-9].ix, xlo, xhi);
				dumpinfo[lineNumber-9].y_dump_unwrapped = unwrapCoordinates (dumpinfo[lineNumber-9].y, dumpinfo[lineNumber-9].iy, ylo, yhi);
				dumpinfo[lineNumber-9].z_dump_unwrapped = unwrapCoordinates (dumpinfo[lineNumber-9].z, dumpinfo[lineNumber-9].iz, zlo, zhi);
				dumpinfo[lineNumber-9].status = 1;

				// /* Checking */
				// fprintf(stdout, "%d ==> %.4f; %.4f; %.4f;\n", lineNumber-9, dumpinfo[lineNumber-9].x_dump_unwrapped, dumpinfo[lineNumber-9].y_dump_unwrapped, dumpinfo[lineNumber-9].z_dump_unwrapped);
				// fflush(stdout);
				// sleep(1);
			}
			if (lineNumber >= (natoms + 8))
			{
				lineNumber = 0;
				return dumpinfo;
			}
		}
	}

	return dumpinfo;
}

void computeCOM (int a, int b, float *x_coords, float *y_coords, float *z_coords, float *x_com, float *y_com, float *z_com)
{
	assert (x_com);
	assert (y_com);
	assert (z_com);

	*x_com = (x_coords[a] + x_coords[b]) / 2;
	*y_com = (y_coords[a] + y_coords[b]) / 2;
	*z_com = (z_coords[a] + z_coords[b]) / 2;
}

float arrayAverage2d (float **array, int x, int y)
{
	float array_sum = 0, array_average, array_denominator = 0;
	for (int i = 0; i < x; ++i)
	{
		for (int j = 0; j < y; ++j)
		{
			array_sum += array[i][j];
			array_denominator++;
		}
	}
	array_average = array_sum / array_denominator;
	return array_average;
}

float computeOrderParameter (int i, int j, int k, int l, float x_dump_unwrapped[], float y_dump_unwrapped[], float z_dump_unwrapped[])
{
	float vx1, vy1, vz1, vx2, vy2, vz2, dotproduct, magnitude1, magnitude2, costheta, theta, LOCAL_ORDER_PARAMETER;

	vx1 = x_dump_unwrapped[i-1] - x_dump_unwrapped[j-1];
	vy1 = y_dump_unwrapped[i-1] - y_dump_unwrapped[j-1];
	vz1 = z_dump_unwrapped[i-1] - z_dump_unwrapped[j-1];
	vx2 = x_dump_unwrapped[k-1] - x_dump_unwrapped[l-1];
	vy2 = y_dump_unwrapped[k-1] - y_dump_unwrapped[l-1];
	vz2 = z_dump_unwrapped[k-1] - z_dump_unwrapped[l-1];

	dotproduct = (vx1 * vx2) + (vy1 * vy2) + (vz1 * vz2);
	magnitude1 = (vx1 * vx1) + (vy1 * vy1) + (vz1 * vz1);
	magnitude2 = (vx2 * vx2) + (vy2 * vy2) + (vz2 * vz2);
	costheta = dotproduct / (sqrt (magnitude1) * sqrt (magnitude2));
	theta = acosf (costheta);
	LOCAL_ORDER_PARAMETER = ((3 * costheta * costheta) - 1) / 2;

	return LOCAL_ORDER_PARAMETER;
}

void *reinitialize_int (int *input)
{
	assert (input);
	*input = 0;
	return 0;
}

void *reinitialize_2int (int *input1, int *input2)
{
	assert (input1);
	assert (input2);
	*input1 = 0;
	*input2 = 0;
	return 0;
}

void *reinitialize_3int (int *input1, int *input2, int *input3)
{
	assert (input1);
	assert (input2);
	assert (input3);
	*input1 = 0;
	*input2 = 0;
	*input3 = 0;
	return 0;
}

void *reinitialize_float (float *input)
{
	assert (input);
	*input = 0;
	return 0;
}

void *reinitialize_2float (float *input1, float *input2)
{
	assert (input1);
	assert (input2);
	*input1 = 0;
	*input2 = 0;
	return 0;
}

void *reinitialize_3float (float *input1, float *input2, float *input3)
{
	assert (input1);
	assert (input2);
	assert (input3);
	*input1 = 0;
	*input2 = 0;
	*input3 = 0;
	return 0;
}

float degreeToRadians(float angle_degrees)
{
	float angle_radians;
	angle_radians = angle_degrees/57.2958;
	return angle_radians;
}

void getDimension_lammpsdump (char *inputDumpfile, float *xlo, float *xhi, float *ylo, float *yhi, float *zlo, float *zhi)
{
	assert (xlo);
	assert (xhi);
	assert (ylo);
	assert (yhi);
	assert (zlo);
	assert (zhi);

	FILE *read = fopen (inputDumpfile, "r");
	char lineString[3000];

	/* Skipping initial lines */
	for (int i = 0; i < 5; ++i)
	{
		fgets (lineString, 3000, read);
	}
	/* Scanning for lo and hi dimensions */
	fscanf (read, "%f %f\n", &(*xlo), &(*xhi));
	fscanf (read, "%f %f\n", &(*ylo), &(*yhi));
	fscanf (read, "%f %f\n", &(*zlo), &(*zhi));

	fclose (read);
}

float unwrapCoordinates (float coordinate, float image, float lo, float hi)
{
	if (image > 0)
	{
		coordinate = hi + ((abs(image) - 1) * (hi - lo)) + (coordinate - lo);
	}
	else if (image < 0)
	{
		coordinate = lo - ((abs(image) - 1) * (hi - lo)) - (hi - coordinate);
	}
	return coordinate;
}

int checkNextTimeframe_lammpsdump (FILE *read)
{
	char lineString[3000];
	int lineNumber = 0, natoms, sino_real, sino_expected = 0;

	if (fgets(lineString, 3000, read) == NULL)
	{
		return 0;
	}
	else
	{
		lineNumber++;
		while (fgets (lineString, 3000, read) != NULL)
		{
			lineNumber++;
			if (lineNumber == 4)
			{
				sscanf(lineString, "%d", &natoms);
			}
			else if (lineNumber > 9 && lineNumber <= (natoms + 9))
			{
				if (strstr(lineString, "ITEM: TIMESTEP"))
				{
					return -1;
				}
				else
				{
					sscanf (lineString, "%d\n", &sino_real);
					sino_expected++;
					if (sino_real != sino_expected)
					{
						return -2;
					}
				}
				if (lineNumber == (natoms + 9))
				{
					lineNumber = 0;
					return 1;
				}
			}
		}
	}
	return 0;
}

char *parseNextTimeframe_lammpsdump (FILE *read, int **atom_sino, int **atom_id, int **atom_type, float **atom_x, float **atom_y, float **atom_z, float **scaled_x, float **scaled_y, float **scaled_z, int **image_x, int **image_y, int **image_z)
{
	/* Repacing the old values in array */
	*atom_sino = NULL;
	*atom_id = NULL;
	*atom_type = NULL;
	*atom_x = NULL;
	*atom_y = NULL;
	*atom_z = NULL;
	*scaled_x = NULL;
	*scaled_y = NULL;
	*scaled_z = NULL;
	*image_x = NULL;
	*image_y = NULL;
	*image_z = NULL;

	char lineString[3000], *successString;
	successString = (char *) malloc (10 * sizeof (char));

	int lineNumber = 0, natoms;

	if (fgets (lineString, 3000, read) == NULL)
	{
		// lineNumber++;
		printf("NULL detected in input dump\n");
		return NULL;
	}
	else
	{
		while (fgets (lineString, 3000, read) != NULL)
		{
			lineNumber++;

			if (lineNumber == 3)
			{
				sscanf (lineString, "%d", &natoms);

				/* Memory allocation */
				*atom_sino = (int *) malloc (natoms * sizeof (int));
				*atom_id = (int *) malloc (natoms * sizeof (int));
				*atom_type = (int *) malloc (natoms * sizeof (int));
				*atom_x = (float *) malloc (natoms * sizeof (float));
				*atom_y = (float *) malloc (natoms * sizeof (float));
				*atom_z = (float *) malloc (natoms * sizeof (float));
				*scaled_x = (float *) malloc (natoms * sizeof (float));
				*scaled_y = (float *) malloc (natoms * sizeof (float));
				*scaled_z = (float *) malloc (natoms * sizeof (float));
				*image_x = (int *) malloc (natoms * sizeof (int));
				*image_y = (int *) malloc (natoms * sizeof (int));
				*image_z = (int *) malloc (natoms * sizeof (int));
			}

			if (lineNumber > 8 && lineNumber <= (natoms + 8))
			{
				sscanf (lineString, "%d %d %d %f %f %f %f %f %f %d %d %d\n", 
					&(*atom_sino)[lineNumber-9],
					&(*atom_id)[lineNumber-9],
					&(*atom_type)[lineNumber-9],
					&(*atom_x)[lineNumber-9],
					&(*atom_y)[lineNumber-9],
					&(*atom_z)[lineNumber-9],
					&(*scaled_x)[lineNumber-9],
					&(*scaled_y)[lineNumber-9],
					&(*scaled_z)[lineNumber-9],
					&(*image_x)[lineNumber-9],
					&(*image_y)[lineNumber-9],
					&(*image_z)[lineNumber-9]);
			}
			if (lineNumber == (natoms + 8))
			{
				lineNumber = 0;
				return successString;
			}
		}
	}
	return NULL;
}

void *parseLastTimeframe_lammmpsdump (const char *inputFileName, int **atom_sino, int **atom_id, int **atom_type, float **atom_x, float **atom_y, float **atom_z, float **scaled_x, float **scaled_y, float **scaled_z, int **image_x, int **image_y, int **image_z)
{
	/* Repacing the old values in array */
	*atom_sino = NULL;
	*atom_id = NULL;
	*atom_type = NULL;
	*atom_x = NULL;
	*atom_y = NULL;
	*atom_z = NULL;
	*scaled_x = NULL;
	*scaled_y = NULL;
	*scaled_z = NULL;
	*image_x = NULL;
	*image_y = NULL;
	*image_z = NULL;

	FILE *read = fopen (inputFileName, "r");
	char lineString[3000];
	int lineNumber = 0, natoms;

	while (fgets (lineString, 3000, read) != NULL)
	{
		lineNumber++;

		if (lineNumber == 4)
		{
			sscanf (lineString, "%d", &natoms);

			/* Memory allocation */
				*atom_sino = (int *) malloc (natoms * sizeof (int));
				*atom_id = (int *) malloc (natoms * sizeof (int));
				*atom_type = (int *) malloc (natoms * sizeof (int));
				*atom_x = (float *) malloc (natoms * sizeof (float));
				*atom_y = (float *) malloc (natoms * sizeof (float));
				*atom_z = (float *) malloc (natoms * sizeof (float));
				*scaled_x = (float *) malloc (natoms * sizeof (float));
				*scaled_y = (float *) malloc (natoms * sizeof (float));
				*scaled_z = (float *) malloc (natoms * sizeof (float));
				*image_x = (int *) malloc (natoms * sizeof (int));
				*image_y = (int *) malloc (natoms * sizeof (int));
				*image_z = (int *) malloc (natoms * sizeof (int));
		}

		if (lineNumber > 9 && lineNumber <= (natoms + 9))
		{
				sscanf (lineString, "%d %d %d %f %f %f %f %f %f %d %d %d\n", 
					&(*atom_sino)[lineNumber-10],
					&(*atom_id)[lineNumber-10],
					&(*atom_type)[lineNumber-10],
					&(*atom_x)[lineNumber-10],
					&(*atom_y)[lineNumber-10],
					&(*atom_z)[lineNumber-10],
					&(*scaled_x)[lineNumber-10],
					&(*scaled_y)[lineNumber-10],
					&(*scaled_z)[lineNumber-10],
					&(*image_x)[lineNumber-10],
					&(*image_y)[lineNumber-10],
					&(*image_z)[lineNumber-10]);
		}
		if (lineNumber == (natoms + 10))
		{
			lineNumber = 0;
		}
	}
	return 0;
}

int isFile(const char *name)
{
  DIR *directory = opendir (name);
  if (directory!=NULL)
  {
    closedir(directory);
    return 0;
  }
  if(errno==ENOTDIR)
  {
    return 1;
  }
  return -1;
}

char *getInputFileName()
{
	int nFiles = 0;
	char *inputFileName, fileExtension[200], terminalString[200];
	int fileRequired;
	inputFileName = (char *) malloc (200 * sizeof (char));

	printf("Enter the file extension or a match string to search in current directory...\n"); 
	fgets (terminalString, sizeof (terminalString), stdin);
	sscanf (terminalString, "%s", fileExtension); 
	printf("\n");
	nFiles = displayFiles(fileExtension);

	if (nFiles > 0)
	{
		printf("\nWhich file would you like to input? Enter a number between (1 to %d): ", nFiles); 
		fgets(terminalString, sizeof (terminalString), stdin);
		sscanf (terminalString, "%d", &fileRequired); 
	}
	else
	{
		printf("No files found with the match string\n"); exit(1);
	}

	nFiles = 0;
	DIR *parentDirectory;
	parentDirectory = opendir ("./");

	struct dirent *filePointer;

	/* Scan all the files using filePointer */
	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileExtension))
		{
			nFiles++;
			if (fileRequired == nFiles)
			{
				strcpy (inputFileName, filePointer -> d_name);
			}
		}
	}
	return inputFileName;
}

int displayFiles(const char *fileExtension)
{
	int nFiles = 0;
	DIR *parentDirectory;
	parentDirectory = opendir ("./");

	struct dirent *filePointer;
	/* Scan all the files using filePointer */
	while ((filePointer = readdir (parentDirectory)))
	{
		if (isFile(filePointer -> d_name) && strstr(filePointer -> d_name, fileExtension))
		{
			nFiles++;
			printf("%d --> %s\n", nFiles, filePointer -> d_name);
		}
	}
	return nFiles;
}

int getNatoms (const char *inputFileName)
{
	FILE *read = fopen (inputFileName, "r");
	char lineString[2000];
	int lineNumber = 0, natoms;

	while (fgets (lineString, 2000, read) != NULL)
	{
		lineNumber++;

		if (lineNumber == 4)
		{
			sscanf (lineString, "%d", &natoms);
			return natoms;
		}
	}
	return 0;
}