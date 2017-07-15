#ifndef HFOLD_H_
#define HFOLD_H_

double hfold(char *sequence, char *restricted, char *structure);
double hfold_emodel(char *sequence, char *restricted, char *structure, std::vector<energy_model> *energy_models);
double hfold_pkonly(char *sequence, char *restricted, char *structure); // April 3, 2012

#endif /*HFOLD_H_*/
