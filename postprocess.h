#ifndef POSTPROCESS_H
#define POSTPROCESS_H

/* UAs in iPP chains are classified into CH3, CH2 and CH, then chain number is stored to differentiate UAs of different chains */
typedef struct isotactic_polypropylene
{
	int CH3, CH2, CH, chain_id_number;
} IPP_DATA;

typedef struct IOPCoordinates_iPP
{
	float x_ab, y_ab, z_ab, x_cd, y_cd, z_cd, x_iop, y_iop, z_iop;
} IOP_COORDS;

typedef struct lammpsdump
{
	int id, molType, atomType, ix, iy, iz;
	float x, y, z, xs, ys, zs, x_dump_unwrapped, y_dump_unwrapped, z_dump_unwrapped, mass;
	int status;
} LAMMPSDUMP;

typedef struct iPP_LOP
{
	int a, b, c, d, e, f;
	float LOP1, LOP2;
} LOP_DATA;

typedef struct iPP_IOP
{
	int a, b, c, d;
	float x1, y1, z1, x2, y2, z2, x_center, y_center, z_center;
} IOP_DATA;

typedef struct iPP_RIOP
{
	float r_min, r_max, RIOP;
} RIOP_DATA;

typedef struct iPP_RIOP_TAVG
{
	float r_min, r_max, riop_sum, RIOP_TAVG;
	int nTimeframes;
} RIOP_TAVG;

typedef struct iPP_OOP
{
	float x1, y1, z1, x2, y2, z2, OOP;
} OOP_DATA;

typedef struct iPP_OOP_TAVG
{
	float oop_sum, OOP_TAVG;
} OOP_TAVG;

// Curve parameters for trignometric functions
typedef struct curve_parameter
{
	int localMinima_Xpos, localMaxima_Xpos;
	double localMinima_Yval, localMaxima_Yval, center;
} CURVE;

typedef struct curve_parameter_highPrecision
{
	float localMinima_Xpos, localMaxima_Xpos;
	float localMinima_Yval, localMaxima_Yval, center;
} CURVE_HIGHPRECISION;

typedef struct cosine
{
	// C3 and Phaseshift are the same. phaseShift is used to write output in user friendly fashion with values from 0 to 90, whereas C3 contains actual values from calculations
	float C1, C2, C3, C4, errorPercent, covariance, phaseShift;
} COSINE;

// Update the cosine parameters in OOP calculations. It adds the next timeframes' data to the same array, which can later be used to calculate average and standard deviation
void updateCurveParameter (COSINE **cosParameters, int nTimeframes, CURVE_HIGHPRECISION curveParameters, float precision, OOP_DATA *oopinfo);

// Check the input OOP and fits the curve with cosine function, returns struct curve_parameter. Check the struct for more details
CURVE findCurveParameters (float *OOP);
CURVE_HIGHPRECISION findCurveParameters_highPrecision (struct iPP_OOP *oopinfo, float precision);

// Pass constants to get back the cosine value
float computeCosineCurve (float c1, float c2, float c3, float c4, float theta);

/* Pass the folder name to create a folder in current working directory */
int makeDirectory (const char *folderName);

/* Prints the final oop average values from function void iPP_computeOOP_TAVG () */
void iPP_printOOP_TAVG (struct iPP_OOP_TAVG *oopinfo_tavg, const char *outputFolder);

/* Computes the running time average of oop values. Address of (struct iPP_OOP_TAVG **oopinfo_tavg) should be passed. */
void iPP_computeOOP_TAVG (struct iPP_OOP *oopinfo, struct iPP_OOP_TAVG **oopinfo_tavg, int nTimeframes);

/* Rotate a line (with two points) along X axis */
float rotateX (float angle_radians);

/* Rotate a line (with two points) along Z axis */
float rotateZ (float angle_radians);

// Compute errorPercent and covariance in cosine curve
float computeError (OOP_DATA *oopinfo, OOP_DATA *modelCurve, float precision);
float computeCovariance (OOP_DATA *oopinfo, OOP_DATA *modelCurve, float precision);

// Find average and standard deviation of cosine values from all timesteps at the end of scan
COSINE findAverage_Cos (COSINE *cosParameters, float precision);
COSINE findSD_Cos (COSINE *cosParameters, float precision, COSINE averageCosine);

/* Find coords of two points of 359 lines. Each lines are rotated 1 degree from each other. Returns the line as type (struct iPP_OOP) */
struct iPP_OOP *iPP_computeOOP_findXYZ ();
struct iPP_OOP *iPP_computeOOP_findXYZ_highPrecision (float precision);

/* Calculates OOP and returns the OOP values as type (struct iPP_OOP). The function (struct iPP_OOP *iPP_computeOOP_findXYZ ()) is required for this function to run. */
struct iPP_OOP *iPP_computeOOP (struct iPP_IOP *iopinfo, int nIOPs);
struct iPP_OOP *iPP_computeOOP_highPrecision (struct iPP_IOP *iopinfo, int nIOPs, float precision, struct iPP_OOP *oopinfo);
struct iPP_OOP *iPP_computeOOP_highPrecision_2D (struct iPP_IOP *iopinfo, int nIOPs, float precision, struct iPP_OOP *oopinfo);
struct iPP_OOP *iPP_computeOOP_highPrecision_2D_bounds (struct iPP_IOP *iopinfo, int nIOPs, float precision, struct iPP_OOP *oopinfo);


/* Compute the average value of riop values. There is no return value, the function changes the values using the address pointer passed. The average values are present in (struct iPP_RIOP_TAVG *riopinfo_tavg) */
void *iPP_computeRIOP_TAVG (struct iPP_RIOP *riopinfo, struct iPP_RIOP_TAVG **riopinfo_tavg, int nTimeframes, float simBoxSize, float dr);

/* Average riop values computed by the above function can be printed using this one. Pass number of RIOPs along with the average struct. FILE pointer is declared within this void function */
void iPP_printRIOP_TAVG (struct iPP_RIOP_TAVG *riopinfo_tavg, int nRIOPs, const char *outputFolder);

/* Compute intermolecular order parameter */
struct iPP_IOP *iPP_computeIOP (struct lammpsdump *dumpinfo, struct isotactic_polypropylene *monomers, int nMonomers, int *nIOPs);

/* Compute radial intermolecluar order parameter */
struct iPP_RIOP *iPP_computeRIOP (struct iPP_IOP *iopinfo, struct lammpsdump *dumpinfo, struct isotactic_polypropylene *monomers, int nIOPs, float dr, float simBoxSize);

/* Scans lammps dump and the iPP chain to return the struct iPP_LOP */
struct iPP_LOP *iPP_computeLOP (struct lammpsdump *dumpinfo, struct isotactic_polypropylene *monomers, int nMonomers, int *nLOPs, float *AVERAGE_LOCAL_ORDER_PARAMETER, const char *outputFolder);

/* Scans through the stored arrays (atomType and molType) to classify UAs based on struct isotactic_polypropylene */
struct isotactic_polypropylene *classifyiPP (int natoms, struct lammpsdump *dumpinfo, int *nMonomers, int *nChains, int nTimeframes);

/* Scans through the lammps dump file and returns type (struct lammpsdump) */
struct lammpsdump *lscanf (FILE *read, int natoms);

/* To display progress bar */
void displayProgressbar (int dx, int x, int progressWidth, const char *processName);

/* Pass X, Y and Z coords to get back center of mass */
void computeCOM (int a, int b, float *x_coords, float *y_coords, float *z_coords, float *x_com, float *y_com, float *z_com);

/* Pass X/Y/Z coordinates, then 4 points (i, j, k and l). The function will return order parameter as type float. */
float computeOrderParameter (int i, int j, int k, int l, float x_dump_unwrapped[], float y_dump_unwrapped[], float z_dump_unwrapped[]);

/* Calculate the average of all elements in a 2d array. Pass the 2d array, rows and cols to get the float of average in return. */
float arrayAverage2d (float **array, int x, int y);

/* Re-initialize multiple variables of type float/int to zero */
void *reinitialize_int (int *input);
void *reinitialize_2int (int *input1, int *input2);
void *reinitialize_3int (int *input1, int *input2, int *input3);
void *reinitialize_float (float *input);
void *reinitialize_2float (float *input1, float *input2);
void *reinitialize_3float (float *input1, float *input2, float *input3);

/* Convert degrees to radians */
float degreeToRadians(float angle_degrees);

/* Get the simulation box size from lammps dump file. Supply a reference to the variables xlo/xhi/ylo/yhi/zlo/zhi along with dump file name */
void getDimension_lammpsdump (char *inputDumpfile, float *xlo, float *xhi, float *ylo, float *yhi, float *zlo, float *zhi);

/* Use X/Y/Z coordinate, image information to unwrap the coordinates in a periodic simulation box */
float unwrapCoordinates (float coordinate, float image, float lo, float hi);

/* Check timeframes for the exact number of natoms */
int checkNextTimeframe_lammpsdump (FILE *read);

/* Parse the next timeframe and load them into arrays passed */
char *parseNextTimeframe_lammpsdump (FILE *read, int **atom_sino, int **atom_id, int **atom_type, float **atom_x, float **atom_y, float **atom_z, float **scaled_x, float **scaled_y, float **scaled_z, int **image_x, int **image_y, int **image_z);

/* Take the last timeframe and modify the arrays for storing dump variables */
void *parseLastTimeframe_lammmpsdump (const char *inputFileName, int **atom_sino, int **atom_id, int **atom_type, float **atom_x, float **atom_y, float **atom_z, float **scaled_x, float **scaled_y, float **scaled_z, int **image_x, int **image_y, int **image_z);

/* Checks if the object is a file or a folder. Returns 1 is true, 0 is false and -1 if not found. */
int isFile(const char *name);

/* Prompts the user to enter the file extension, uses displayFiles() function to display matching files, then return the file name the user wants as input */
char *getInputFileName();

/* Checks all files/folders in the current directory, checks for matching string. This function prints all the matching files and returns the number of matching files in the present directory */
/* Direct use of this function is not recommended, since this is part of *getInputFileName() function */
int displayFiles(const char *fileExtension);

/* This function takes file name as input and returns the number of atoms in the dump file. NOTE that this does not check the entire file, just the header lines to get the number of atoms present. Sometimes the number of atoms mentioned in header lines can be very different from the number of actual atoms present, in which case, caution must be taken */
int getNatoms (const char *inputFileName);

#endif