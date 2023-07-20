// Dinoknot files
#include "dinoknot.hh"
#include "cmdline.hh"
#include "Result.hh"

// Simfold files
#include "simfold.h"
#include "externs.h"
#include "h_globals.h"
#include "params.h"
#include "s_specific_functions.h"
#include "Hotspot.h"
#include "h_common.h"

// Non user files
#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <vector>
#include <sys/stat.h>



#define KEVIN_DEBUG 1
/*
As requested, any blocks of changes have "//AP" near them.

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
void validateStructure(std::string sequence, std::string structure){
	if(structure.length() != sequence.length()){
		std::cout << " The length of the sequence and corresponding structure must have the same length" << std::endl;
		exit(EXIT_FAILURE);
	}

	//check if any characters are not ._()
	for(char c : structure) {
		if (!(c == '.' || c == '_' || c == '(' || c == ')')){
			std::cout << "Structure must only contain ._(): " << c << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}

//check if sequence is valid with regular expression
//check length and if any characters other than GCAUT
void validateSequence(std::string sequence){

	if(sequence.length() == 0){
		std::cout << "sequence1 or sequence2 is missing" << std::endl;
		exit(EXIT_FAILURE);
	}
  // return false if any characters other than GCAUT -- future implement check based on type
  for(char c : sequence) {
    if (!(c == 'G' || c == 'C' || c == 'A' || c == 'U' || c == 'T')) {
		std::cout << "Sequence contains character " << c << " that is not G,C,A,U, or T." << std::endl;
		exit(EXIT_FAILURE);
    }
  }
}

bool exists (const std::string path) {
  struct stat buffer;   
  return (stat (path.c_str(), &buffer) == 0); 
}

void get_input(std::string file, std::string &sequence1, std::string &sequence2, std::string &structure1, std::string &structure2 ){
	if(!exists(file)){
		std::cout << "Input file does not exist" << std::endl;
		exit(EXIT_FAILURE);
	}
	// Partial IUPAC notation
	std::string bases = "ACGTUWSMKRY";
	std::string s = "._()";
	std::ifstream in(file.c_str());
	bool secondseq = false;
	bool secondstruct = false;
	std::string str;
	while(getline(in,str)){
		if(bases.find(str[0])+1){
			if(secondseq) sequence2 = str;
			else {
				sequence1 = str; 
				secondseq = true;
			}
		}
		if(s.find(str[0])+1){
			
			if(secondstruct) structure2 = str;
			else {
				structure1 = str; 
				secondstruct = true;
			}
		}
	}

	in.close();
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

int main (int argc, char *argv[]) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	int model_1_Type = args_info.type1_given;
	int model_2_Type = args_info.type2_given;
	std::string inputSequence1 = "";
	std::string inputSequence2 = "";
	std::string inputStructure1 = "";
	std::string inputStructure2 = "";

	std::string inputFile = args_info.input_given ? input_file : "";
	if(args_info.input_given) get_input(inputFile,inputSequence1,inputSequence2,inputStructure1,inputStructure2);
	std::cout << inputSequence1 << "\n" << inputSequence2 << "\n" << inputStructure1 << "\n" << inputStructure2 << std::endl;

	if(args_info.sequence1_given) inputSequence1 =  sequence_1;
	if(args_info.sequence2_given) inputSequence2 =  sequence_2;
	if(args_info.structure1_given) inputStructure1 = structure_1;
	if(args_info.structure2_given) inputStructure2 = structure_2;
	validateSequence(inputSequence1);
	validateSequence(inputSequence2);
	if(inputStructure1 != "") validateStructure(inputSequence1,inputStructure1);
	if(inputStructure2 != "") validateStructure(inputSequence2,inputStructure2);
	
	std::string seq = inputSequence1 + "XXXXX" + inputSequence2;

				
	std::string outputDir = args_info.dir_given ? output_dir : "";
	std::string outputFile = args_info.output_given ? output_file : "";
	std::string hotspotDir = args_info.h_only_given ? hotspot_dir : "";

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
		// Hosna, July 18, 2012
		// initialize the thermodynamic parameters
		init_data_emodel (argv[0], energy_model.config_file.c_str(), energy_model.dna_or_rna, energy_model.temperature, &energy_model);

		// In simfold we have the following for RNA && temp=37
		fill_data_structures_with_new_parameters_emodel (SIMFOLD_HOME "/params/turner_parameters_fm363_constrdangles.txt", &energy_model);

		// Hosna, July 25, 2012
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

	std::vector<Hotspot> hotspot_list1;
	std::vector<Hotspot> hotspot_list2;
	// Hotspot* hotspot;
	if(inputStructure1 != ""){
		Hotspot hotspot1(0,inputStructure1.length()-1,inputStructure1.length());
		hotspot1.set_structure(inputStructure1);
		hotspot_list1.push_back(hotspot1);
	}else{
		int length = inputSequence1.length();
		char sequence1[length+1];
		strcpy(sequence1,inputSequence1.c_str());
		get_hotspots(sequence1, hotspot_list1,max_hotspot);
	}

	if(inputStructure2 != ""){
		
		Hotspot hotspot2(0,inputStructure2.length()-1,inputStructure2.length());
		hotspot2.set_structure(inputStructure2);
		hotspot_list2.push_back(hotspot2);
		
	}else{
		int length = inputSequence2.length();
		char sequence2[length+1];
		strcpy(sequence2,inputSequence2.c_str());
		get_hotspots(sequence2, hotspot_list2,max_hotspot);
	}

	if(hotspot_only){
		std::ofstream out(hotspotDir.c_str());
		if(!exists(hotspot_dir)){
			std::cout << "Input File does not exist!" << std::endl;
			exit (EXIT_FAILURE);
    	}
		for(int i =0; i < hotspot_list1.size(); i++){
			out << "Seq1_hotspot_" << i << ": " << hotspot_list1[i].get_structure() << "(" << hotspot_list1[i].get_energy() << ")" << std::endl;
		}
		out << "---------------" << std::endl;
		for(int j = 0; j < hotspot_list2.size(); j++){
			out << "Seq2_hotspot_" << j << ": " << hotspot_list2[j].get_structure() << "(" << hotspot_list2[j].get_energy() << ")" << std::endl;
		}
		out.close();
	}

	std::vector<Result> result_list;
	
	int length = seq.length();
	char sequence[length+1];
	strcpy(sequence,seq.c_str());
	for(int i =0; i < hotspot_list1.size(); i++){
		for(int j = 0; j < hotspot_list2.size(); j++){
			
			char structure[length+1];
			char restricted[length+1];

			std::string struc = hotspot_list1[i].get_structure() + "....." + hotspot_list2[j].get_structure();
			strcpy(restricted, struc.c_str());
			

			double energy = hfold_interacting_emodel(sequence, restricted, structure, energy_models);
			
			std::string res(restricted);
			std::string final(structure);
			Result result(seq,res,final,energy);
			result_list.push_back(result);
		}
	}
	Result::Result_comp result_comp;
	std::sort(result_list.begin(), result_list.end(),result_comp );

	destruct_energy_model(&model_1);
	destruct_energy_model(&model_2);

	//kevin 5 oct 2017
	int number_of_output = 1;
	// //printf("number_of_suboptimal_structure: %d\n",number_of_suboptimal_structure);
	if(number_of_suboptimal_structure != 1){
		number_of_output = std::min( (int) result_list.size(),number_of_suboptimal_structure);
	}

	//Mateo 7/19/2023
	//output to file
	if(outputFile != ""){
		std::ofstream out(output_file);
		if(!exists(output_file)){
			std::cout << "file is not valid" << std::endl;
			exit(EXIT_FAILURE);
		}

		out << "Seq:          " << seq << std::endl;
		for (int i=0; i < number_of_output; i++) {
			out << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
			out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
		
		}
		out.close();
	}
	else if(outputDir != ""){
		// Mateo 2023
		if(exists(outputDir)){
			if(outputDir[outputDir.length()] != '/') outputDir += '/';
			for (int i=0; i < number_of_output; ++i) {
				std::string path_to_file = outputDir + "output_" + std::to_string(i) + ".txt";
				std::ofstream out(path_to_file);
				out << "Seq:          " << seq << std::endl;
				out << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
				out << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;  
				out.close();
			}
		}
		else{
			std::cout << "Not a valid output directory" << std::endl;
			exit(EXIT_FAILURE);
		}
	} else{
		// Mateo 2023
		std::cout << "Seq:          " << seq << std::endl;
		for (int i=0; i < number_of_output; i++) {
			std::cout << "Restricted_" << i << ": " << result_list[i].get_restricted() << std::endl;;
			std::cout << "Result_" << i << ":     " << result_list[i].get_final_structure() << " (" << result_list[i].get_final_energy() << ")" << std::endl;
		}
	}

	return 0;
}
