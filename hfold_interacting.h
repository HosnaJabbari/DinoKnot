#ifndef HFOLD_INTERACTING_H_
#define HFOLD_INTERACTING_H_


double hfold_interacting(char *sequence, char *restricted, char *structure);
bool write_output_file(char* path_to_file, int num_of_output, std::vector<Result*> result_list);
void printUsage();
int validateModelType(char* type);  
double get_START_HYBRID_PENALTY(int type1, int type2);
int isDirectory(const char *path);
bool write_directory(char* path_to_dir, int num_of_output, std::vector<Result*> result_list);

#endif /*HFOLD_INTERACTING_H_*/
