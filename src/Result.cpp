#include "Result.hh"

//constructor
Result::Result(std::string sequence,std::string restricted, std::string final_structure, double final_energy){
    this->sequence = sequence;
    this->restricted = restricted;
    this->final_structure = final_structure;
    this->final_energy = final_energy;

}

//destructor
Result::~Result(){
   
}

std::string Result::get_sequence(){
    return this->sequence;
}
std::string Result::get_restricted(){
    return this->restricted;
}
std::string Result::get_final_structure(){
    return this->final_structure;
}
double Result::get_final_energy(){
    return this->final_energy;
}