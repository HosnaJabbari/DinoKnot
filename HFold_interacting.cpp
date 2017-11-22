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

//kevin 26 Sept 2017
#include "s_specific_functions.h"
#include "Hotspot.h"
#include "h_common.h"

#include "hfold_interacting.h"

#include <sys/types.h>
#include <sys/stat.h>



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

	char* outputDir;
	outputDir = (char*) malloc(sizeof(char) * 1000);

	int model_1_Type = -1;
	int model_2_Type = -1;

	bool sequence1Found = false;
	bool structure1Found = false;
	bool sequence2Found = false;
	bool structure2Found = false;
	bool inputPathFound = false;
	bool outputPathFound = false;
	bool outputDirFound = false;
	bool errorFound = false;
	bool type1Found = false;
	bool type2Found = false;

	int number_of_suboptimal_structure = 0;

	int option;

	int max_hotspot = 20; //default value for number of hotspot

	START_HYBRID_PENALTY = -1;

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
				{"n"  ,required_argument, 0, 'h'}, 	//number of suboptimal structure
				{"o_dir", required_argument, 0, 'j'},	//output every file to directory
				{"hotspot_num", required_argument, 0, 'k'}, //max number of hotspot
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
					fprintf(stderr, "--s1 is duplicated\n");
					errorFound = true;
					break;
				}
				if(inputPathFound){
					fprintf(stderr, "Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				strcpy(inputSequence1,optarg);
				if(!validateSequence(inputSequence1,true)){
					fprintf(stderr, "--s1 is invalid\n");
					errorFound = true;
					break;
				}
				sequence1Found = true;
				break;
			case 'b': //--r1 (structure1)
				if(structure1Found){
					fprintf(stderr, "--r1 is duplicated\n");
					errorFound = true;
					break;
				}
				if(inputPathFound){
					fprintf(stderr, "Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				strcpy(inputStructure1,optarg);
				structure1Found = true;
				break;
			case 'c': //--s2 (sequence2)
				if(sequence2Found){
					fprintf(stderr, "--s2 is duplicated\n");
					errorFound = true;
					break;
				}
				if(inputPathFound){
					fprintf(stderr, "Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				strcpy(inputSequence2,optarg);
				if(!validateSequence(inputSequence2,true)){
					fprintf(stderr, "--s2 is invalid\n");
					errorFound = true;
					break;
				}
				sequence2Found = true;
				break;
			case 'd': //--r2 (structure2)
				if(structure2Found){
					fprintf(stderr, "--r2 is duplicated\n");
					errorFound = true;
					break;
				}
				if(inputPathFound){
					fprintf(stderr, "Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				strcpy(inputStructure2,optarg);
				structure2Found = true;
				break;
			case 'i':
				if(sequence1Found || structure1Found || sequence2Found || structure2Found){
					fprintf(stderr, "Cannot combine -i with --s1/--r1/--s2/--r2 \n");
					errorFound = true;
					break;
				}
				/*
				if(!validateInteractingInputFile(optarg, inputSequence1, inputStructure1, inputSequence2, inputStructure2)){ //store sequences, structures in corresponding variables if valid
					printf("Input file is invalid\n");
					errorFound = true;
					break;
				}
				*/
				if(!validateInteractingInputFile2(optarg, inputSequence1, inputStructure1, inputSequence2, inputStructure2,&sequence1Found, &structure1Found , &sequence2Found , &structure2Found)){
					fprintf(stderr, "Input file is invalid\n");
					errorFound = true;
					break;
				}
				strcpy(inputPath, optarg);

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
						fprintf(stderr, "--t1 is duplicated\n");
						errorFound = true;
						break;
				}
				model_1_Type = validateModelType(optarg);
				if(model_1_Type < 0){
					fprintf(stderr, "Model Type is invalid\n");
					errorFound = true;
					break;
				}
				type1Found = true;
				break;
			case 'f':
				if(type2Found){
					fprintf(stderr, "--t2 is duplicated\n");
					errorFound = true;
					break;
				}
				model_2_Type = validateModelType(optarg);
				if(model_2_Type < 0){
					fprintf(stderr, "Model Type is invalid\n");
					errorFound = true;
					break;
				}
				type2Found = true;
				break;
			case 'g':
				START_HYBRID_PENALTY = atof(optarg);
				break;
			case 'h':
				number_of_suboptimal_structure = atoi(optarg);
				if(number_of_suboptimal_structure <= 0){
					fprintf(stderr, "number must be > 0\n");
					errorFound = true;
					break;
				}
				break;
			case 'j':
				if(!isDirectory(optarg)){
					fprintf(stderr, "argument is not a directory or directory does not exist\n");
					errorFound = true;
					break;
				}
				strcpy(outputDir,optarg);
				//printf("outdir: %s\n",outputDir);
				//printf("optarg: %s\n",optarg);
				outputDirFound = true;
				break;
			case 'k':
				if(atoi(optarg) <= 0){
					fprintf(stderr, "number of hotspot must be > 0\n");
					errorFound = true;
					break;
				}
				max_hotspot = atoi(optarg);
				break;
			default:
				errorFound = true;
				break;
			}
			//clean up when error
			if(errorFound){
				free(inputPath);
				free(outputPath);
				free(outputDir);
				printUsage();
				exit(1);
			}
	}


	if(!sequence1Found || !sequence2Found){
		fprintf(stderr, "sequence1 or sequence2 is missing\n");
		free(inputPath);free(outputPath);free(outputDir);
		printUsage();
		exit(1);
	}

	if(!type1Found){
		fprintf(stderr, "--t1 is missing\n");
		free(inputPath);free(outputPath);free(outputDir);
		printUsage();
		exit(1);
	}

	if(!type2Found){
		fprintf(stderr, "--t2 is missing\n");
		free(inputPath);free(outputPath);free(outputDir);
		printUsage();
		exit(1);
	}

	//kevin: if we have output path and input path, try to combine both
	if(outputPathFound && inputPathFound){
		addPath(&outputPath, inputPath);
		//printf("out path: %s\n",outputPath);

	}

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

	if (START_HYBRID_PENALTY == -1)
	    START_HYBRID_PENALTY = get_START_HYBRID_PENALTY(model_1_Type,model_2_Type);
	//printf("penalty: %lf\n",START_HYBRID_PENALTY);
	//exit(999);

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

	//kevin: june 23 2017
	//end of validation for command line arguments

	// Hosna: November 16, 2015
	// changed this part accordingly
	linker_pos = strlen(inputSequence1);

	strcpy(sequence, inputSequence1);
	if(KEVIN_DEBUG){
		strcat(sequence, "XXXXX");
	}else{
		strcat(sequence, linker);
	}
	strcat(sequence, inputSequence2);

	//set up for RNA so we can use this for building hotspot
	if(!structure1Found || !structure2Found){
		char config_file[200];
		strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");
		int dna_or_rna;
		dna_or_rna = RNA;
		double temperature = 37.0;
		init_data ("./HFold", config_file, dna_or_rna, temperature);
		fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");
		fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/parameters_DP09.txt");
	}

	std::vector<Hotspot*> hotspot_list1;
	std::vector<Hotspot*> hotspot_list2;
	Hotspot* hotspot;
	if(structure1Found){
		if(!(validateStructure(inputStructure1,inputSequence1,true))){
			fprintf(stderr, "--r1 is invalid\n");
			free(inputPath);free(outputPath);
			printUsage();
			exit(1);
		}else{
			hotspot = new Hotspot(0,strlen(inputStructure1)-1,strlen(inputStructure1));
			hotspot->set_structure(inputStructure1);
			hotspot_list1.push_back(hotspot);
		}
	}else{
		//printf("start hotspot for struct1\n");
		get_hotspots(inputSequence1, &hotspot_list1,max_hotspot);
	}

	if(structure2Found){
		if(!(validateStructure(inputStructure2,inputSequence2,true))){
			fprintf(stderr, "--r2 is invalid\n");
			free(inputPath);free(outputPath);free(outputDir);
			printUsage();
			exit(1);
		}else{
			hotspot = new Hotspot(0,strlen(inputStructure2)-1,strlen(inputStructure2));
			hotspot->set_structure(inputStructure2);
			hotspot_list2.push_back(hotspot);
		}
	}else{
		//printf("start hotspot for struct2\n");
		get_hotspots(inputSequence2, &hotspot_list2,max_hotspot);
	}


	Result* result;
	std::vector<Result*> result_list;
	for(int i =0; i < hotspot_list1.size(); i++){
		for(int j = 0; j < hotspot_list2.size(); j++){
			//printf("i: %d j:%d\n",i,j);
			//printf("hot1: %s\n",hotspot_list1[i]->get_structure());
			//printf("hot2: %s\n",hotspot_list2[j]->get_structure());
			strcpy(restricted, hotspot_list1[i]->get_structure());
			strcat(restricted, ".....");
			strcat(restricted, hotspot_list2[j]->get_structure());

			int method_used = -1;
			energy = hfold_interacting_emodel(sequence, restricted, structure, &energy_models, method_used);

			result = new Result(sequence,restricted, structure,energy, method_used);
			result_list.push_back(result);
		}
	}
	std::sort(result_list.begin(), result_list.end(),compare_result_ptr);

	//kevin 5 oct 2017
	int number_of_output;
	//printf("number_of_suboptimal_structure: %d\n",number_of_suboptimal_structure);
	if(number_of_suboptimal_structure != 0){
		number_of_output = MIN(result_list.size(),number_of_suboptimal_structure);
	}else{
		number_of_output = 1;
	}

	//kevin: june 22 2017
	//output to file
	if(outputPathFound){
		bool write_success = write_output_file(outputPath, number_of_output, result_list);
		if(!write_success){
			exit(4);
		}
	}else if(outputDirFound){
		bool write_success = write_directory(outputDir, number_of_output, result_list);
		if(!write_success){
			exit(4);
		}
	}else{
		//kevin 5 oct 2017
		printf("Seq: %s\n",sequence);
		for (int i=0; i < number_of_output; i++) {
			printf("Restricted_%d: %s\n",i, result_list[i]->get_restricted());
			printf("Result_%d: %s \nEnergy_%d: %lf\n",i, result_list[i]->get_final_structure(),i,result_list[i]->get_final_energy());
		}
	}


    // clean up
    free(inputPath);
	free(outputPath);
	free(outputDir);

	for (int i=0; i < result_list.size(); i++) {
		delete result_list[i];
	}

	for(int i =0; i<hotspot_list1.size(); i++){
		delete hotspot_list1[i];
	}

	for(int i =0; i<hotspot_list2.size(); i++){
		delete hotspot_list2[i];
	}

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

bool write_output_file(char* path_to_file, int num_of_output, std::vector<Result*> result_list){
	FILE* fp = fopen(path_to_file,"w");
	if (fp == NULL) {
        perror("Write to file error:");
		return false;
	}
	fprintf(fp,"Seq: %s\n",result_list[0]->get_sequence());
	for (int i=0; i < num_of_output; i++) {
		fprintf(fp,"Restricted_%d: %s\n",i, result_list[i]->get_restricted());
		fprintf(fp,"Result_%d: %s \nEnergy_%d: %lf\n",i, result_list[i]->get_final_structure(),i,result_list[i]->get_final_energy());
	}
	fclose(fp);
	return true;
}


bool write_directory(char* path_to_dir, int num_of_output, std::vector<Result*> result_list){
	for (int i=0; i < num_of_output; i++) {
		char* path_to_file = (char*)malloc(sizeof(int)*1000);
		sprintf(path_to_file, "%s/output_%d", path_to_dir, i);
		//printf("path_to_file: %s\n",path_to_file);
		FILE* fp = fopen(path_to_file,"w");
		if (fp == NULL) {
			perror("Write to directory file error:");
			return false;
		}
		fprintf(fp,"Seq: %s\n",result_list[0]->get_sequence());
		fprintf(fp,"Restricted_%d: %s\n",i, result_list[i]->get_restricted());
		fprintf(fp,"Result_%d: %s \nEnergy_%d: %lf\n",i, result_list[i]->get_final_structure(),i,result_list[i]->get_final_energy());
		free(path_to_file);
		fclose(fp);
	}

	return true;
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

double get_START_HYBRID_PENALTY(int type1, int type2){
	if(type1 == type2){ //if both model are the same
		if(type1 == RNA){ //if both are RNA
			return 3.0;
		}else if(type1 == DNA){ //if both are DNA
			return 4.4785825;
		}else if(type1 == PMO){
			fprintf(stderr, "ERROR: model cannot be both PMO\n");
			exit(1);
		}
	}else{
		return 58.4511432; //when 2 different model
	}
}

int isDirectory(const char *path) {
   struct stat statbuf;
   if (stat(path, &statbuf) != 0)
       return 0;
   return S_ISDIR(statbuf.st_mode);
}
