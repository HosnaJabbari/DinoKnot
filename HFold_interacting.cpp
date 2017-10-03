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
//#include "h_externs.h"
#include "constants.h"
#include "params.h"
#include "common.h"

// Hosna June 20th, 2007
//#include "W_final.h"
//#include "hfold_interacting.h"
#include "hfold.h"

//kevin 23 June 2017
#include <getopt.h>
#include "hfold_validation.h"
#include <unistd.h>
#include "h_common.h"

void printUsage();

int validateModelType(char* type);

#define KEVIN_DEBUG 1
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

	if (sequence[i] == X || sequence[j] == X || sequence[ip] == X || sequence[jp] == X)
		return 0;
	if (sequence[i+1] == X || sequence[j+1] == X || sequence[ip+1] == X || sequence[jp+1] == X || sequence[i-1] == X || sequence[j-1] == X || sequence[ip-1] == X || sequence[jp-1] == X)
		return INF;

These if statements are used in the s_..._loop.cpp files to check for the nucleotide 'X'. This has to be done manually because the program will crash if you try to change the NUCL variable. The variable is set currently set to four but if you try to change it to five to tell the program that there is an extra 'X' nucleotide then it crashes due to the various hard coded functions throughout simfold. Since the 'X' nucleotide is set to four in the sequence, we have to manually check for instances of 'X' and act accordingly.
*/

int main (int argc, char *argv[]) {
	char sequence[MAXSLEN];
	char structure[MAXSLEN];
	char restricted[MAXSLEN];
	double energy;

	// A vector is used to store the energy models in order to keep generality throughout the rest of the program. This means that you can add extra energy models withour having to change the code of how many energy models to loop through.
	std::vector<energy_model> energy_models;
	energy_model *model_1;
	energy_model *model_2;

/*
	if (argc != 5) {
		//printf ("Usage: %s <oligo_sequence> <oligo_structure> <gene_sequence> <gene_structure>\n", argv[0]);
		printf ("Usage: %s <sequence_one> <structure_one> <sequence_two> <structure_two>\n", argv[0]);
		printf ("Example: %s \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" \"(____________________________)\" \"GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC\" \"(__________________________________________________________)\"\n", argv[0]);
		printf ("  Restricted structure symbols:\n");
		printf ("    () restricted base pair\n");
		printf ("    _ no restriction\n");
		return 0;
	}
*/

	//kevin: june 23 2017
	//validation for command line argument
	char inputSequence1[MAXSLEN];
	char inputSequence2[MAXSLEN];
	char inputStructure1[MAXSLEN];
	char inputStructure2[MAXSLEN];
	char* inputPath;
	inputPath = (char*) malloc(sizeof(char) * 1000);

	char* outputPath;
	outputPath = (char*) malloc(sizeof(char) * 1000);

	int model_1_Type = -1;
	int model_2_Type = -1;

	bool sequence1Found = false;
	bool structure1Found = false;
	bool sequence2Found = false;
	bool structure2Found = false;
	bool inputPathFound = false;
	bool outputPathFound = false;
	bool errorFound = false;
	bool type1Found = false;
	bool type2Found = false;

	int option;

	START_HYBRID_PENALTY = 310.0; //default value for hybrid penalty, can be changed with --pen

	//kevin: june 23 2017 https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html
	while (1){
		static struct option long_options[] =
			{
				{"s1", required_argument, 0, 'a'}, 	//sequence1
				{"r1", required_argument, 0, 'b'},	//structure for sequence1
				{"s2", required_argument, 0, 'c'},	//sequence2
				{"r2", required_argument, 0, 'd'},	//structure for sequence2
				{"t1", required_argument, 0, 'e'},	//type for sequence1
				{"t2", required_argument, 0, 'f'},	//type for sequence2
				{"pen",required_argument, 0, 'g'},  //start_hybrid_penalty
				{0, 0, 0, 0}
			};
		// getopt_long stores the option index here.
		int option_index = 0;

		option = getopt_long (argc, argv, "a:b:c:d:i:o:e:f:",
						long_options, &option_index);

		// Detect the end of the options.
		if (option == -1)
			break;

		switch (option)
			{
			case 'a': //--s1 (sequence1)
				if(sequence1Found){
					printf("--s1 is duplicated\n");
					errorFound = true;
					break;
				}
				if(inputPathFound){
					printf("Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				strcpy(inputSequence1,optarg);
				if(!validateSequence(inputSequence1)){
					printf("--s1 is invalid\n");
					errorFound = true;
					break;
				}
				sequence1Found = true;
				break;
			case 'b': //--r1 (structure1)
				if(structure1Found){
					printf("--r1 is duplicated\n");
					errorFound = true;
					break;
				}
				if(inputPathFound){
					printf("Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				strcpy(inputStructure1,optarg);
				structure1Found = true;
				break;
			case 'c': //--s2 (sequence2)
				if(sequence2Found){
					printf("--s2 is duplicated\n");
					errorFound = true;
					break;
				}
				if(inputPathFound){
					printf("Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				strcpy(inputSequence2,optarg);
				if(!validateSequence(inputSequence2)){
					printf("--s2 is invalid\n");
					errorFound = true;
					break;
				}
				sequence2Found = true;
				break;
			case 'd': //--r2 (structure2)
				if(structure2Found){
					printf("--r2 is duplicated\n");
					errorFound = true;
					break;
				}
				if(inputPathFound){
					printf("Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				strcpy(inputStructure2,optarg);
				structure2Found = true;
				break;
			case 'i':
				if(sequence1Found || structure1Found || sequence2Found || structure2Found){
					printf("Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				if(!validateInteractingInputFile(optarg, inputSequence1, inputStructure1, inputSequence2, inputStructure2)){ //store sequences, structures in corresponding variables if valid
					printf("Input file is invalid\n");
					errorFound = true;
					break;
				}
				strcpy(inputPath, optarg);
				inputPathFound = true;
				break;
			case 'o':
				strcpy(outputPath, optarg);
				//printf("access: %d\n",access(outputPath, F_OK));
				if(access(outputPath, F_OK) != -1) { //if file already exist
					addTimestamp(&outputPath);
				}
				outputPathFound = true;
				break;
			case 'e':
				if(type1Found){
						printf("--t1 is duplicated\n");
						errorFound = true;
						break;
				}
				model_1_Type = validateModelType(optarg);
				if(model_1_Type < 0){
					printf("Model Type is invalid\n");
					errorFound = true;
					break;
				}
				type1Found = true;
				break;
			case 'f':
				if(type2Found){
					printf("--t2 is duplicated\n");
					errorFound = true;
					break;
				}
				model_2_Type = validateModelType(optarg);
				if(model_2_Type < 0){
					printf("Model Type is invalid\n");
					errorFound = true;
					break;
				}
				type2Found = true;
				break;
			case 'g':
				START_HYBRID_PENALTY = atof(optarg);
				break;
			default:
				errorFound = true;
				break;
			}
			//clean up when error
			if(errorFound){
				free(inputPath);
				free(outputPath);
				printUsage();
				exit(1);
			}
	}

	if(!inputPathFound){
		//if sequence1 or sequence2 or structure1 or structure2 is missing when input file is not present
		if(!(sequence1Found && structure1Found && sequence2Found && structure2Found)){
			fprintf(stderr, "--s1/--r1/--s2/--r2 is missing\n");
			free(inputPath);free(outputPath);
			printUsage();
			exit(1);
		}else{
			//validate both structures
			if(!(validateStructure(inputStructure1,inputSequence1))){
				fprintf(stderr, "--r1 is invalid\n");
				free(inputPath);free(outputPath);
				printUsage();
				exit(1);
			}else{
				replaceBrackets(inputStructure1);
			}
			if(!(validateStructure(inputStructure2,inputSequence2))){
				fprintf(stderr, "--r2 is invalid\n");
				free(inputPath);free(outputPath);
				printUsage();
				exit(1);
			}else{
				replaceBrackets(inputStructure2);
			}
		}
	}

	if(!type1Found){
		fprintf(stderr, "--t1 is missing\n");
		free(inputPath);free(outputPath);
		printUsage();
		exit(1);
	}

	if(!type2Found){
		fprintf(stderr, "--t2 is missing\n");
		free(inputPath);free(outputPath);
		printUsage();
		exit(1);
	}

	//kevin: if we have output path and input path, try to combine both
	if(outputPathFound && inputPathFound){
		addPath(&outputPath, inputPath);
		//printf("out path: %s\n",outputPath);

	}

	//kevin: june 23 2017
	//end of validation for command line arguments

	// Hosna: November 16, 2015
	// changed this part accordingly
	linker_pos = strlen(inputSequence1);

	strcpy(sequence, inputSequence1);
	if(KEVIN_DEBUG){
		strcat(sequence, "XXXXX");
		//linker_length = 0;
	}else{
		strcat(sequence, linker);
	}
	strcat(sequence, inputSequence2);

	strcpy(restricted, inputStructure1);
	strcat(restricted, ".....");
	strcat(restricted, inputStructure2);

	// Before calling any function in the library, you have to initialize config_file, dna_or_rna, temperature
	// and to call the function init_data, which loads the thermodynamic parameters into memory

	// To add more energy models, just copy this block of code and change the parameters as needed.
	model_1 = new energy_model();
	init_energy_model(model_1); // Initializes the data structures in the energy model.
	model_1->config_file = SIMFOLD_HOME "/params/multirnafold.conf"; // configuration file, the path should be relative to the location of this executable
	model_1->dna_or_rna = model_1_Type; // what to fold: RNA or DNA
	model_1->temperature = 37.0; // temperature: any integer or real number between 0 and 100 Celsius
	energy_models.push_back(*model_1);

	model_2 = new energy_model();
	init_energy_model(model_2); // Initializes the data structures in the energy model.
	model_2->config_file = SIMFOLD_HOME "/params/multirnafold.conf"; // configuration file, the path should be relative to the location of this executable
	model_2->dna_or_rna = model_2_Type; // what to fold: RNA or DNA
	model_2->temperature = 37.0; // temperature: any integer or real number between 0 and 100 Celsius
	energy_models.push_back(*model_2);

	for (auto &energy_model : energy_models) {
		// initialize the thermodynamic parameters
		// call init_data only once for the same dna_or_rna and same temperature
		// if one of them changes, call init_data again
		init_data_emodel (argv[0], energy_model.config_file.c_str(), energy_model.dna_or_rna, energy_model.temperature, &energy_model);

		// Hosna, July 18, 2012
		// In simfold we have the following for RNA && temp=37
		fill_data_structures_with_new_parameters_emodel (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt", &energy_model);

		// Hosna, July 25, 2012
		// in HotKnots and ComputeEnergy package the most up-to-date parameters set is DP09.txt
		// so we add it here
		fill_data_structures_with_new_parameters_emodel (SIMFOLD_HOME "/params/parameters_DP09.txt", &energy_model);
	}

    int method_used = -1;
	//energy = hfold_emodel(sequence, restricted, structure, &energy_models);
	//kevin july 13 changed to call hfold_interacting_emodel instead of hfold_emodel
	energy = hfold_interacting_emodel(sequence, restricted, structure, &energy_models, method_used);

	if(outputPathFound){
		FILE* fp;
		fp = fopen(outputPath,"w");
		if(fp){
			fprintf(fp,"Sequence: %s\n",sequence);
			fprintf(fp,"Input_structure: %s\n",restricted);
			fprintf(fp,"Output_structure: %s\n",structure);
			fprintf(fp,"Energy: %.2lf\n",energy);
			fprintf(fp,"Method Used: %d",method_used);
			fclose(fp);
		}
	}else{
		printf ("Seq: %s\n", sequence);
    	printf ("RES: %s  %.2lf\n", structure, energy);
	}

    // clean up
    free(inputPath);
	free(outputPath);

	destruct_energy_model(model_1);
	destruct_energy_model(model_2);
	delete model_1;
	delete model_2;

	return 0;
}

void printUsage(){
	printf("\nUsage ./HFold_interacting_multimodel --s1 <sequence1> --r2 <restricted_structure1> --s1 <sequence_2> --r2 <restricted_structure2> --t1 <type_for_sequence1> --t2 <type_for_sequence1> [-o </path/to/file>]\n");
	printf("or\n");
	printf("Usage ./HFold -i </path/to/file> --t1 <type_for_sequence1> --t2 <type_for_sequence1> [-o </path/to/file>]\n");
	printf ("  Restricted structure symbols:\n");
	printf ("    () restricted base pair\n");
	printf ("    _ no restriction\n");
	printf("Example:\n");
	printf("./HFold_interacting_multimodel --s1 \"GCAACGAUGACAUACAUCGCUAGUCGACGC\" --r1 \"(____________________________)\" --s2 \"GCAACGAUGACAUACAUCGCUAGUCGACGCGCAACGAUGACAUACAUCGCUAGUCGACGC\" --r2 \"(__________________________________________________________)\" --t1 RNA --t2 DNA\n");
	printf("./HFold_interacting_multimodel -i \"/home/username/Desktop/myinputfile.txt\" -o \"/home/username/Desktop/some_folder/outputfile.txt\"\n");
	printf("Please read README for more details\n");

}

//return code for model type if valid
//return -1 if not valid
int validateModelType(char* type){
	if(strcmp(type, "RNA") == 0){
		return RNA;
	}
	if(strcmp(type, "DNA") == 0){
		return DNA;
	}
	//kevin 30 June 2017
	//added PMO
	if(strcmp(type, "PMO") == 0){
		return PMO;
	}
	return -1;
}
