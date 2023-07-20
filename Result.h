#ifndef RESULT_HEADER
#define RESULT_HEADER

#include <stdio.h>  
#include <assert.h> 
#include <stdlib.h>
#include <string>

class Result{
    public:
        //constructor
        Result(std::string sequence,std::string restricted, std::string final_structure, double final_energy);
        //destructor
        ~Result();

        //getter
        std::string get_sequence();
        std::string get_restricted();
        std::string get_final_structure();
        double get_final_energy();

        struct Result_comp{
		bool operator ()(Result &x, Result &y) const {
			return x.get_final_energy() < y.get_final_energy();
		}
		} result_comp;

        
    private:
        std::string sequence;
        std::string restricted;
        std::string final_structure;
        double final_energy;
};



#endif
