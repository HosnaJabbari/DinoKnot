// a simple driver for the HFold
#include "dinoknot.hh"
#include "cmdline.hh"
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <algorithm>

// include the simfold header files
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
#include "constants.h"
#include "params.h"
#include "common.h"

#include "hfold.h"

//kevin 23 June 2017
#include "hfold_validation.h"
#include "h_common.h"

//kevin 26 Sept 2017
#include "s_specific_functions.h"
#include "Hotspot.h"
#include "h_common.h"



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
bool validateStructure(std::string sequence, std::string structure){
	if(structure.length() != sequence.length()){
		std::cout << " The length of the sequence and corresponding structure must have the same length" << std::endl;
		return false;
	}

	//check if any characters are not ._()
	for(char c : structure) {
		if (!(c == '.' || c == '_' || c == '(' || c == ')')){
			std::cout << "Structure must only contain ._(): " << c << std::endl;
			return false;
		}
	}

	return true;
}

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
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
			return 22.96551130344778;
		}else if(type1 == DNA){ //if both are DNA
			return 34.14979695798525;
		}else if(type1 == PMO){
			fprintf(stderr, "ERROR: model cannot be both PMO\n");
			exit(1);
		}
	}
	
	return 166.0; //when 2 different model old: 58.4511432
	
}

int isDirectory(const char *path) {
   struct stat statbuf;
   if (stat(path, &statbuf) != 0)
       return 0;
   return S_ISDIR(statbuf.st_mode);
}

int main (int argc, char *argv[]) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string inputSequence1 = args_info.sequence1_given ? sequence_1 : "";
	std::string inputSequence2 = args_info.sequence2_given ? sequence_2 : "";

	if(!args_info.sequence1_given || !args_info.sequence2_given){
		std::cout << "sequence1 or sequence2 is missing" << std::endl;
		exit(EXIT_FAILURE);
	}
	std::string seq = inputSequence1 + "XXXXX" + inputSequence2;

	std::string inputStructure1 = args_info.structure1_given ? structure_1 : "";
	if(args_info.structure1_given) 
		if(!validateStructure(inputSequence1,inputStructure1)){
			std::cout << "--r1 is invalid" << std::endl;
			exit(EXIT_FAILURE);
		}
		
	std::string inputStructure2 = args_info.structure2_given ? structure_2 : "";
	if(args_info.structure2_given) 
		if(!validateStructure(inputSequence2,inputStructure2)){
			std::cout << "--r2 is invalid" << std::endl;
			exit(EXIT_FAILURE);
		}
		

	std::string outputDir = args_info.dir_given ? output_dir : "";
	std::string outputFile = args_info.dir_given ? output_file : "";
	std::string hotspotDir = args_info.h_only_given ? hotspot_dir : "";

	int model_1_Type = args_info.type1_given;
	int model_2_Type = args_info.type2_given;

	int max_hotspot = args_info.h_num_given ? hotspot_num : 20;
	int number_of_suboptimal_structure = args_info.subopt_given ? subopt : 1;

	bool hotspot_only = args_info.h_only_given;

	START_HYBRID_PENALTY = args_info.pen_given ? hybrid_pen : get_START_HYBRID_PENALTY(model_1_Type,model_2_Type);

	cmdline_parser_free(&args_info);

	linker_pos = inputSequence1.length();

	// A vector is used to store the energy models in order to keep generality throughout the rest of the program. This means that you can add extra energy models withour having to change the code of how many energy models to loop through.
	std::vector<energy_model> energy_models;

	energy_model model_1;
	init_energy_model(&model_1); // Initializes the data structures in the energy model.
	model_1.config_file = SIMFOLD_HOME "/params/multirnafold.conf"; // configuration file, the path should be relative to the location of this executable
	model_1.dna_or_rna = model_1_Type; // what to fold: RNA or DNA
	model_1.temperature = 37.0; // temperature: any integer or real number between 0 and 100 Celsius
	energy_models.push_back(model_1);

	energy_model model_2;
	init_energy_model(&model_2); // Initializes the data structures in the energy model.
	model_2.config_file = SIMFOLD_HOME "/params/multirnafold.conf"; // configuration file, the path should be relative to the location of this executable
	model_2.dna_or_rna = model_2_Type; // what to fold: RNA or DNA
	model_2.temperature = 37.0; // temperature: any integer or real number between 0 and 100 Celsius
	energy_models.push_back(model_2);
	

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

	//set up for RNA so we can use this for building hotspot
	if(!args_info.structure1_given || !args_info.structure2_given){
		char config_file[200];
		strcpy (config_file, SIMFOLD_HOME "/params/multirnafold.conf");
		int dna_or_rna = RNA;
		double temperature = 37.0;
		init_data ("./HFold", config_file, dna_or_rna, temperature);
		fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt");
		fill_data_structures_with_new_parameters ( SIMFOLD_HOME "/params/parameters_DP09.txt");
	}

	std::vector<Hotspot*> hotspot_list1;
	std::vector<Hotspot*> hotspot_list2;
	Hotspot* hotspot;
	if(inputStructure1 != ""){
		char structure1[inputStructure1.length()+1];
		strcpy(structure1,inputStructure1.c_str());
		hotspot = new Hotspot(0,strlen(structure1)-1,strlen(structure1));
		hotspot->set_structure(structure1);
		hotspot_list1.push_back(hotspot);
	}else{
		int length = inputSequence1.length();
		char sequence1[length+1];
		strcpy(sequence1,inputSequence1.c_str());
		get_hotspots(sequence1, &hotspot_list1,max_hotspot);
	}

	if(inputStructure2 != ""){
		
		char structure2[inputStructure2.length()];
		strcpy(structure2,inputStructure2.c_str());
		hotspot = new Hotspot(0,strlen(structure2)-1,strlen(structure2));
		hotspot->set_structure(structure2);
		hotspot_list2.push_back(hotspot);
	}else{
		int length = inputSequence2.length();
		char sequence2[length+1];
		strcpy(sequence2,inputSequence2.c_str());
		get_hotspots(sequence2, &hotspot_list2,max_hotspot);
	}

	if(hotspot_only){
		FILE* fp = fopen(hotspotDir.c_str(),"w");
		if(fp == NULL){
			fprintf(stderr, "cannot write to hotspot file\n");
			exit(4);
		}
		for(int i =0; i < hotspot_list1.size(); i++){
			fprintf(fp,"Seq1_hotspot_%d: %s\n",i,hotspot_list1[i]->get_structure());
		}
		fprintf(fp,"-----\n");
		for(int j = 0; j < hotspot_list2.size(); j++){
			fprintf(fp,"Seq2_hotspot_%d: %s\n",j,hotspot_list2[j]->get_structure());
		}
		fclose(fp);
		exit(0);
	}

	Result* result;
	std::vector<Result*> result_list;
	
	int length = seq.length();
	char sequence[length+1];
	strcpy(sequence,seq.c_str());
	for(int i =0; i < hotspot_list1.size(); i++){
		for(int j = 0; j < hotspot_list2.size(); j++){
			
			char structure[length+1];
			char restricted[length+1];
			strcpy(restricted, hotspot_list1[i]->get_structure());
			strcat(restricted, ".....");
			strcat(restricted, hotspot_list2[j]->get_structure());

			double energy = hfold_interacting_emodel(sequence, restricted, structure, energy_models);

			result = new Result(sequence,restricted, structure,energy);
			result_list.push_back(result);
		}
	}
	std::sort(result_list.begin(), result_list.end(),compare_result_ptr);

	//kevin 5 oct 2017
	int number_of_output = 1;
	//printf("number_of_suboptimal_structure: %d\n",number_of_suboptimal_structure);
	if(number_of_suboptimal_structure != 1){
		number_of_output = MIN(result_list.size(),number_of_suboptimal_structure);
	}

	//Mateo 7/19/2023
	//output to file
	if(outputFile != ""){
		std::ofstream out(output_file);

		out << "Seq:          " << seq << std::endl;
		for (int i=0; i < number_of_output; i++) {
			out << "Restricted_" << i << ": " << result_list[i]->get_restricted() << std::endl;;
			out << "Result_" << i << ":     " << result_list[i]->get_final_structure() << " (" << result_list[i]->get_final_energy() << ")" << std::endl;
		
		}
		out.close();
	}
	// if(args_info.odir_given){
	// 	for (int i=0; i < number_of_output; ++i) {
	// 		char* path_to_file = (char*)malloc(sizeof(int)*1000);
	// 		sprintf(path_to_file, "%s/output_%d", outputDir.c_str(), i);
	// 		//printf("path_to_file: %s\n",path_to_file);
	// 		FILE* fp = fopen(path_to_file,"w");
	// 		if (fp == NULL) {
	// 			perror("Write to directory file error:");
	// 			return false;
	// 		}
	// 		fprintf(fp,"Seq: %s\n",result_list[0]->get_sequence());
	// 		fprintf(fp,"Restricted_%d: %s\n",i, result_list[i]->get_restricted());
	// 		fprintf(fp,"Result_%d: %s \nEnergy_%d: %lf\n",i, result_list[i]->get_final_structure(),i,result_list[i]->get_final_energy());
	// 		free(path_to_file);
	// 		fclose(fp);
	// 	} }
	else{
	// 	//kevin 5 oct 2017
			std::cout << "Seq:          " << seq << std::endl;
		for (int i=0; i < number_of_output; i++) {
			std::cout << "Restricted_" << i << ": " << result_list[i]->get_restricted() << std::endl;;
			std::cout << "Result_" << i << ":     " << result_list[i]->get_final_structure() << " (" << result_list[i]->get_final_energy() << ")" << std::endl;
		}
	}


    // clean up

	for (int i=0; i < result_list.size(); i++) {
		delete result_list[i];
	}

	for(int i =0; i<hotspot_list1.size(); i++){
		delete hotspot_list1[i];
	}

	for(int i =0; i<hotspot_list2.size(); i++){
		delete hotspot_list2[i];
	}

	destruct_energy_model(&model_1);
	destruct_energy_model(&model_2);

	return 0;
}
