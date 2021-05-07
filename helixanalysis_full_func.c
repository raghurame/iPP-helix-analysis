#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "postprocess.h"

int main(int argc, char const *argv[])
{
	FILE *readdump, *readdihedral, *write, *writeRandL, *chiralorderparameterlog;
	char *inputDump, *inputDihedral;

	printf("Searching for input dump file...\n");
	inputDump = (char *) malloc (200 * sizeof (char));
	inputDump = getInputFileName();
	printf("\n\nSearching for input dihedral file...\n");
	inputDihedral = (char *) malloc (200 * sizeof (char));
	inputDihedral = getInputFileName();

	readdump = fopen (inputDump, "r");
	readdihedral = fopen (inputDihedral, "r");
	chiralorderparameterlog = fopen ("chiral_order_parameter.log", "w");
	write = fopen ("VisualizationCombined.lammpstrj", "w");
	writeRandL = fopen ("RandL", "w");

	/*

	VARABLE INITIALIZATION STARTS

	*/

	int MAX_TIMEFRAMES_TO_SCAN = 3000;

	int natoms, ndihedrals, currentLine_overall = 0, currentLine_dump = 0, currentLine_dihedral = 0, currentTimeframe = 0;
	char lineString_dump[4000], lineString_dihedral[4000];

	// Initializing variables to store dihedral parameters
	int *sino_dih, *atom1_dih, *atom2_dih, *atom3_dih, *atom4_dih;
	float *angle_dih;
	// Initializing variables to store sorted dihedral parameters
	int *sino_dih_sorted, *atom1_dih_sorted, *atom2_dih_sorted, *atom3_dih_sorted, *atom4_dih_sorted;
	float *angle_dih_sorted;
	int sortArrayID = 0;
	// Binary parameters to store T/G information
	int *isTrans, *isGauchePlus, *isGaucheMinus;
	int *isTGPlus, *isTGMinus;
	int index_dih = 0, index_dump = 0;

	// Initializing variables to store dump parameters
	int *sino_dump, *mol_dump, *type_dump, *ix_dump, *iy_dump, *iz_dump;
	float *x_dump, *y_dump, *z_dump, *xs_dump, *ys_dump, *zs_dump;
	// Binary variable to check if the chain belongs to main backbone or not
	int *isBackbone;

	// Array to store atom information about all backbone atoms
	int *sino_dump_bb, *mol_dump_bb, *type_dump_bb, *ix_dump_bb, *iy_dump_bb, *iz_dump_bb;
	float *x_dump_bb, *y_dump_bb, *z_dump_bb, *xs_dump_bb, *ys_dump_bb, *zs_dump_bb;
	int natoms_bb;
	// Array to store information about all backbone dihedrals
	int *sino_dih_sorted_bb, *atom1_dih_sorted_bb, *atom2_dih_sorted_bb, *atom3_dih_sorted_bb, *atom4_dih_sorted_bb;
	float *angle_dih_sorted_bb;
	int ndihedrals_bb;

	// Variables to check dihedral occurances as a function of Y-axis
	// frequency_within and frequency_outside count the overall R- and L-handed helices
	// nRight and nLeft counts the R- and L-handed helices in each timeframe alone
	int *frequency_within, *frequency_outside;
	frequency_within = (int *) calloc (3, sizeof (int));
	frequency_outside = (int *) calloc (3, sizeof (int));
	// To store running values
	int frequency_within_none_currentvalue = 0, frequency_within_left_currentvalue = 0, frequency_within_right_currentvalue = 0;
	int frequency_outside_none_currentvalue = 0, frequency_outside_left_currentvalue = 0, frequency_outside_right_currentvalue = 0;

	float ymin = 0, ymax = 0;

	if (strstr(argv[1], "compressed"))
	{
		ymin = -10000; ymax = -10000;
		printf("ymin: %.0f; ymax: %.0f;\n", ymin, ymax);
	}

	if (strstr(argv[1], "uncompressed"))
	{
		ymin = -30; ymax = -13;
		printf("ymin: %.0f; ymax: %.0f;\n", ymin, ymax);
	}

	if (ymin == 0 || ymax == 0)
	{
		printf("No args passed\n");
		exit(1);
	}

	fprintf(chiralorderparameterlog, "chiral_order_parameter_overall, chiral_order_parameter_inside, chiral_order_parameter_outside\n");

	// float ymin = -10000, ymax = 10000;
	// float ymin = -30, ymax = -13;
	int nRight, nLeft;

	// printf("ymin: %.0f; ymax: %.0f;\n", ymin, ymax);

	// Memory allocation for 'frequency_within'
	/*
	 * frequency_within[0] = number of dihedrals within ymin to ymax with atom type: 1
	 * frequency_within[1] = number of dihedrals within ymin to ymax with atom type: 2
	 * frequency_within[2] = number of dihedrals within ymin to ymax with atom type: 3
	 */

	// Ratios and fractions
	float *ratio_overall, *ratio_right, *ratio_left, *fraction_dihedrals_overall, *fraction_dihedrals_right_overall, *fraction_dihedrals_left_overall, *fraction_dihedrals_right_inside, *fraction_dihedrals_left_inside, *fraction_dihedrals_right_outside, *fraction_dihedrals_left_outside;
	// Chiral order parameter
	float *chiral_order_parameter_inside, *chiral_order_parameter_outside, *chiral_order_parameter_overall;
	float chiral_order_parameter_overall_running, chiral_order_parameter_inside_running, chiral_order_parameter_outside_running;

	// Memory declaration for ratio and fractions above, based on MAX_TIMEFRAMES_TO_SCAN
	ratio_overall = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	ratio_right = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	ratio_left = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	fraction_dihedrals_overall = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	fraction_dihedrals_right_overall = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	fraction_dihedrals_left_overall = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	fraction_dihedrals_right_inside = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	fraction_dihedrals_left_inside = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	fraction_dihedrals_right_outside = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	fraction_dihedrals_left_outside = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	chiral_order_parameter_inside = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	chiral_order_parameter_outside = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));
	chiral_order_parameter_overall = (float *) calloc (MAX_TIMEFRAMES_TO_SCAN, sizeof (float));

	// Average and sum
	float ratio_overall_sum = 0, ratio_right_sum = 0, ratio_left_sum = 0, fraction_dihedrals_overall_sum = 0, fraction_dihedrals_right_overall_sum = 0, fraction_dihedrals_left_overall_sum = 0, fraction_dihedrals_right_inside_sum = 0, fraction_dihedrals_left_inside_sum = 0, fraction_dihedrals_right_outside_sum = 0, fraction_dihedrals_left_outside_sum = 0;
	float ratio_overall_average = 0, ratio_right_average = 0, ratio_left_average = 0, fraction_dihedrals_overall_average = 0, fraction_dihedrals_right_overall_average = 0, fraction_dihedrals_left_overall_average = 0, fraction_dihedrals_right_inside_average = 0, fraction_dihedrals_left_inside_average = 0, fraction_dihedrals_right_outside_average = 0, fraction_dihedrals_left_outside_average = 0;

	float chiral_order_parameter_inside_sum = 0, chiral_order_parameter_outside_sum = 0, chiral_order_parameter_overall_sum = 0;
	float chiral_order_parameter_inside_average = 0, chiral_order_parameter_outside_average = 0, chiral_order_parameter_overall_average = 0;

	// For standard deviation - main variables to store
	float ratio_overall_sd = 0, ratio_right_sd = 0, ratio_left_sd = 0, fraction_dihedrals_overall_sd = 0, fraction_dihedrals_right_overall_sd = 0, fraction_dihedrals_left_overall_sd = 0, fraction_dihedrals_right_inside_sd = 0, fraction_dihedrals_left_inside_sd = 0, fraction_dihedrals_right_outside_sd = 0, fraction_dihedrals_left_outside_sd = 0;
	float chiral_order_parameter_inside_sd = 0, chiral_order_parameter_outside_sd = 0, chiral_order_parameter_overall_sd = 0;

	// Declaring sdv1 - variable for standard deviation calculation
	float ratio_overall_sdv1 = 0, ratio_right_sdv1 = 0, ratio_left_sdv1 = 0, fraction_dihedrals_overall_sdv1 = 0, fraction_dihedrals_right_overall_sdv1 = 0, fraction_dihedrals_left_overall_sdv1 = 0, fraction_dihedrals_right_inside_sdv1 = 0, fraction_dihedrals_left_inside_sdv1 = 0, fraction_dihedrals_right_outside_sdv1 = 0, fraction_dihedrals_left_outside_sdv1 = 0;
	float chiral_order_parameter_inside_sdv1 = 0, chiral_order_parameter_outside_sdv1 = 0, chiral_order_parameter_overall_sdv1 = 0;

	// Declaring sdv2 - variable for standard deviation calculation; sdv2 is the summation variable
	float ratio_overall_sdv2 = 0, ratio_right_sdv2 = 0, ratio_left_sdv2 = 0, fraction_dihedrals_overall_sdv2 = 0, fraction_dihedrals_right_overall_sdv2 = 0, fraction_dihedrals_left_overall_sdv2 = 0, fraction_dihedrals_right_inside_sdv2 = 0, fraction_dihedrals_left_inside_sdv2 = 0, fraction_dihedrals_right_outside_sdv2 = 0, fraction_dihedrals_left_outside_sdv2 = 0;
	float chiral_order_parameter_inside_sdv2 = 0, chiral_order_parameter_outside_sdv2 = 0, chiral_order_parameter_overall_sdv2 = 0;

	// Helix length and frequency
	int index_helixlength = 0, *frequency_helixlength;

	// Extra variables for ease of calculations and for code readability
	float numerator1, numerator2, denominator1, denominator2;

	// Boolean
	int isFirstTimeframe = 1;

	printf("\n");

	fprintf(writeRandL, "Time R-handed L-handed\n");

	/*

	VARABLE INITIALIZATION ENDS

	*/

	
	/* 

	MEMORY ALLOCATION STARTS 

	*/
	natoms = getNatoms(inputDump);

	sino_dump_bb = (int *) malloc (natoms * sizeof (int));
	mol_dump_bb = (int *) malloc (natoms * sizeof (int));
	type_dump_bb = (int *) malloc (natoms * sizeof (int));
	x_dump_bb = (float *) malloc (natoms * sizeof (float));
	y_dump_bb = (float *) malloc (natoms * sizeof (float));
	z_dump_bb = (float *) malloc (natoms * sizeof (float));
	xs_dump_bb = (float *) malloc (natoms * sizeof (float));
	ys_dump_bb = (float *) malloc (natoms * sizeof (float));
	zs_dump_bb = (float *) malloc (natoms * sizeof (float));
	ix_dump_bb = (int *) malloc (natoms * sizeof (int));
	iy_dump_bb = (int *) malloc (natoms * sizeof (int));
	iz_dump_bb = (int *) malloc (natoms * sizeof (int));

	// MEMORY ALLOCATION FOR DIHEDRAL VARIABLES
	sino_dih = malloc (ndihedrals * sizeof (int));
	atom1_dih = malloc (ndihedrals * sizeof (int));
	atom2_dih = malloc (ndihedrals * sizeof (int));
	atom3_dih = malloc (ndihedrals * sizeof (int));
	atom4_dih = malloc (ndihedrals * sizeof (int));
	angle_dih = malloc (ndihedrals * sizeof (float));

	sino_dih_sorted = (int *) malloc (ndihedrals * sizeof (int));
	atom1_dih_sorted = (int *) malloc (ndihedrals * sizeof (int));
	atom2_dih_sorted = (int *) malloc (ndihedrals * sizeof (int));
	atom3_dih_sorted = (int *) malloc (ndihedrals * sizeof (int));
	atom4_dih_sorted = (int *) malloc (ndihedrals * sizeof (int));
	angle_dih_sorted = (float *) malloc (ndihedrals * sizeof (float));

	sino_dih_sorted_bb = (int *) malloc (ndihedrals * sizeof (int));
	atom1_dih_sorted_bb = (int *) malloc (ndihedrals * sizeof (int));
	atom2_dih_sorted_bb = (int *) malloc (ndihedrals * sizeof (int));
	atom3_dih_sorted_bb = (int *) malloc (ndihedrals * sizeof (int));
	atom4_dih_sorted_bb = (int *) malloc (ndihedrals * sizeof (int));
	angle_dih_sorted_bb = (float *) malloc (ndihedrals * sizeof (float));

	// Trans and Gauche are used to mark dihedrals
	isTrans = (int *) calloc (ndihedrals, sizeof (int));
	isGauchePlus = (int *) calloc (ndihedrals, sizeof (int));
	isGaucheMinus = (int *) calloc (ndihedrals, sizeof (int));

	// TG+ and TG- are markers for monomers (atoms)
	isTGPlus = (int *) calloc (natoms, sizeof (int));
	isTGMinus = (int *) calloc (natoms, sizeof (int));

	// Memory allocation for helix length calculation
	frequency_helixlength = (int *) calloc (natoms, sizeof (int));

	/* 

	MEMORY ALLOCATION ENDS 

	*/


	int nTimeframes = 0;

	while (parseNextTimeframe(readdump, &sino_dump, &mol_dump, &type_dump, &x_dump, &y_dump, &z_dump, &xs_dump,  &ys_dump, &zs_dump, &ix_dump, &iy_dump, &iz_dump) != NULL)
	{
		/* Variable initialization / re-initialization */
		frequency_within_none_currentvalue = 0;
		frequency_within_right_currentvalue = 0;
		frequency_within_left_currentvalue = 0;
		frequency_outside_none_currentvalue = 0;
		frequency_outside_right_currentvalue = 0;
		frequency_outside_left_currentvalue = 0;
		numerator1 = 0;
		numerator2 = 0;
		denominator1 = 0;
		denominator2 = 0;
		nRight = 0;
		nLeft = 0;
		index_dump = 0;
		index_dih = 0;
		currentLine_overall++;
	}

	fprintf(stdout, "nTimeframes: %d\n", nTimeframes);
	fflush (stdout);

	return 0;
}