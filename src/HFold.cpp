// a simple driver for the HFold
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <algorithm>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
#include "h_common.h"
//#include "h_externs.h"
#include "structs.h"
#include "constants.h"
#include "params.h"

// Hosna June 20th, 2007
//#include "W_final.h"
#include "hfold.h"
#include "dinoknot.h"

/*
As requested, any blocks of changes have "//AP" near them.

Throughout the program there will be multiple lines that look like the following:

	for (auto &energy_model : *energy_models) {
		energy_model.energy_value = H->compute_energy_restricted_emodel (i, j, fres, &energy_model);
	}
	min_en[0] = emodel_energy_function (i, j, energy_models);

The for loop is used to go iterate through the energy model vector and use each energy model to
perform a certain calculation, such as computing the energy for H. This value is stored inside the
energy model to be used in emodel_energy_function. Once a value has been calculated for each energy
model, the energy model vector will be passed into emodel_energy_function which can be found in
common.cpp. This function is used to calculate the average value depending on where in the sequence the
i and j value are. This will be changed once the PMO energy table is implemented. This can be changed
if checks are done before the calculation occurs to determine whether or not certain energy models will
just get a value of zero. Doing this will save time on a running extra iterations of the code that will
only give you zero anyways.

	if (sequence[i] == 4 || sequence[j] == 4 || sequence[ip] == 4 || sequence[jp] == 4)
		return 0;
	if (sequence[i+1] == 4 || sequence[j+1] == 4 || sequence[ip+1] == 4 || sequence[jp+1] == 4 || sequence[i-1] == 4 || sequence[j-1] == 4 || sequence[ip-1] == 4 || sequence[jp-1] == 4)
		return INF;

These if statements are used in the s_..._loop.cpp files to check for the nucleotide 'X'. This has to be done manually because the program will crash if you try to change the NUCL variable. The variable is set currently set to four but if you try to change it to five to tell the program that there is an extra 'X' nucleotide then it crashes due to the various hard coded functions throughout simfold. Since the 'X' nucleotide is set to four in the sequence, we have to manually check for instances of 'X' and act accordingly.
*/

int main (int argc, char *argv[]) {
    char sequence[MAXSLEN];
    char structure[MAXSLEN];
    char restricted[MAXSLEN];
    double energy;
    double energies[MAXSUBSTR];

	// A vector is used to store the energy models in order to keep generality throughout the rest of the program. This means that you can add extra energy models withour having to change the code of how many energy models to loop through.
	std::vector<energy_model> energy_models;
	energy_model *model;

    if (argc != 3)
    {
		//"GCAACGAUGACAUACAUCGCUAGUCGACGC" "(____________________________)"
        printf ("Usage: %s <sequence> <restricted_structure>\n", argv[0]);
        printf ("Example: %s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" \"(____________________________)\"\n", argv[0]);
        printf ("\tRestricted structure symbols:\n");
        printf ("\t\t() restricted base pair\n");
        printf ("\t\t_ no restriction\n");
        return 0;
    }

	strcpy (sequence, argv[1]);
    strcpy (restricted, argv[2]);

	// Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
    // and to call the function init_data, which loads the thermodynamic parameters into memory

	// To add more energy models, just copy this block of code and change the parameters as needed.
	model = new energy_model();
	init_energy_model(model); // Initializes the data structures in the energy model.
	model->config_file = "./simfold/params/multirnafold.conf"; // configuration file, the path should be relative to the location of this executable
	model->dna_or_rna = RNA; // what to fold: RNA or DNA
	model->temperature = 37.0; // temperature: any integer or real number between 0 and 100 Celsius
	energy_models.push_back(*model);

	for (auto &energy_model : energy_models) {
		// initialize the thermodynamic parameters
		// call init_data only once for the same dna_or_rna and same temperature
		// if one of them changes, call init_data again
		init_data_emodel (argv[0], energy_model.config_file.c_str(), energy_model.dna_or_rna, energy_model.temperature, &energy_model);

		// Hosna, July 18, 2012
		// In simfold we have the following for RNA && temp=37
		fill_data_structures_with_new_parameters_emodel ("./simfold/params/turner_parameters_fm363_constrdangles.txt", &energy_model);

		// Hosna, July 25, 2012
		// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
		// so we add it here
		fill_data_structures_with_new_parameters_emodel ("./simfold/params/parameters_DP09.txt", &energy_model);
	}

	energy = hfold_emodel(sequence, restricted, structure, &energy_models);

    //delete min_fold;

    // check if restricted is included in structure
	// Hosna March 7, 2012
	// for optimality we get the value once, use it many times
	int seqLen = strlen (sequence);
    for (int i=0; i < seqLen; i++)
    {
        if ((restricted[i] == '(' || restricted[i] == ')' || restricted[i] == '.') &&
            (restricted[i] != structure[i]))
        {
            fprintf (stderr, "There is something wrong with the structure, doesn't match restricted\n");
			fprintf (stderr, "  %s\n  %s\n  %s\t%.2lf\n", sequence, restricted, structure, energy);
			exit(1);
        }
    }

	printf ("Seq: %s\n", sequence);
    printf ("RES: %s  %.2lf\n", structure, energy);

	// Call the destructor for each energy model.
	for (auto &energy_model : energy_models) {
		destruct_energy_model(&energy_model);
	}

	// Clean up the energy model vector that contains N number of energy models.
	energy_models.erase(std::remove_if(energy_models.begin(), energy_models.end(), [&](energy_model const & emodel) {return &emodel!=NULL; }), energy_models.end());

    return 0;
}
